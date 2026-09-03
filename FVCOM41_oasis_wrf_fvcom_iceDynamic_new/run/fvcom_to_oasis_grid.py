#!/usr/bin/env python
"""
fvcom_to_oasis.py
==================

Build OASIS3-MCT coupling grid files (grids.nc, areas.nc, masks.nc) from an
FVCOM grid/output NetCDF file.

Background
----------
OASIS3-MCT (and the oasisgrids.py / esmgrids tooling used for structured
models like MOM/NEMO) expects three files:

    grids.nc   <name>.lat / <name>.lon / <name>.cla / <name>.clo
    areas.nc   <name>.srf
    masks.nc   <name>.msk

For structured grids these arrays are 2-D (ny, nx). FVCOM is unstructured,
so it is exposed to OASIS as a "reduced"/1-D grid: ny = 1, nx = number of
points (nodes or elements). This is the same convention used for other
unstructured-mesh couplings (FESOM, ICON, MPAS, WW3 unstructured, etc.).

FVCOM has two natural meshes:
  * NODE mesh   - scalars live here (zeta, temp, salinity). Dimension "node".
  * ELEMENT mesh - vectors live here (u, v, wind stress). Dimension "nele",
                    triangles defined by the "nv" connectivity array.

This script builds an OASIS grid for each:
  * <grid_name[:3]>n  -> node grid   (point + reconstructed dual-cell corners)
  * <grid_name[:3]>e  -> element grid (triangle, corners = its 3 nodes)

Assumptions / things to check against your actual file
--------------------------------------------------------
1. Coordinate variables: script looks for lon/lat first (spherical case),
   falling back to x/y (Cartesian/UTM case -- OASIS needs degrees, so if
   your FVCOM run is Cartesian-only you must supply/derive lon/lat separately
   before using this for real coupling).
2. Connectivity: "nv" with shape (nele, 3) or (3, nele), values either
   1-based (typical FVCOM/Fortran output) or 0-based -- auto-detected.
3. Areas: uses "art1" (node control-volume area) and "art" (element/triangle
   area) if present in the file; otherwise computes planar/spherical polygon
   area from the reconstructed corners.
4. Mask: defaults to all-unmasked (0 everywhere, OASIS convention:
   1 = masked/land, 0 = unmasked/valid). If a wet/dry variable
   (wet_nodes / wet_cells, 1=wet) is present, mask = 1 - wet_*.
5. Node dual-cell corners are approximated as the centroids of the
   elements surrounding each node, ordered by polar angle. This is a
   standard approximation for unstructured coupling meshes but is not
   exact at open/land boundary nodes (the fan of surrounding elements is
   incomplete there); refine if you need strictly conservative remapping
   at the domain boundary.

Usage
-----
    python fvcom_to_oasis.py fvcom_input.nc \
        --grid_name fvo \
        --grids grids.nc --areas areas.nc --masks masks.nc \
        [--which both|node|element] [--mask_var wet_nodes]

Run it once per FVCOM domain; like oasisgrids.py, it *appends* new grid
variables into existing grids.nc/areas.nc/masks.nc so you can add more
model grids afterwards.
"""

import argparse
import os

import numpy as np
import netCDF4 as nc


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def pick_coords(ds):
    """Return (lon, lat, lonc, latc) 1-D node/element coordinate arrays."""
    node_names = [("lon", "lat"), ("x", "y")]
    elem_names = [("lonc", "latc"), ("xc", "yc")]

    lon = lat = lonc = latc = None
    for a, b in node_names:
        if a in ds.variables and b in ds.variables:
            lon = ds.variables[a][:].astype("f8")
            lat = ds.variables[b][:].astype("f8")
            break
    for a, b in elem_names:
        if a in ds.variables and b in ds.variables:
            lonc = ds.variables[a][:].astype("f8")
            latc = ds.variables[b][:].astype("f8")
            break

    if lon is None:
        raise ValueError("Could not find node coordinates (lon/lat or x/y) in file.")

    return lon, lat, lonc, latc


def get_nv(ds, nele):
    """Return 0-based (nele, 3) node-connectivity array."""
    nv = ds.variables["nv"][:]
    nv = np.asarray(nv)
    if nv.shape[0] == 3 and nv.shape[1] == nele:
        nv = nv.T
    elif nv.shape == (nele, 3):
        pass
    else:
        raise ValueError("Unexpected nv shape {}".format(nv.shape))

    # Auto-detect 1-based vs 0-based indexing.
    if nv.min() >= 1:
        nv = nv - 1
    return nv.astype("int64")


def element_centroids(nv, lon, lat):
    return lon[nv].mean(axis=1), lat[nv].mean(axis=1)


def build_node_polygons(nv, n_nodes, lon, lat, lonc, latc):
    """
    Reconstruct a polygon of surrounding-element centroids for every node
    (the FVCOM "dual" cell). Returns corner_lon, corner_lat of shape
    (n_nodes, max_corners), padded by repeating the last valid corner.
    """
    node_to_elems = [[] for _ in range(n_nodes)]
    for e, tri in enumerate(nv):
        for n in tri:
            node_to_elems[n].append(e)

    max_corners = max(len(v) for v in node_to_elems)
    max_corners = max(max_corners, 3)

    clon = np.full((n_nodes, max_corners), np.nan)
    clat = np.full((n_nodes, max_corners), np.nan)

    for i in range(n_nodes):
        elems = node_to_elems[i]
        if len(elems) == 0:
            # Isolated node (shouldn't happen) - degenerate polygon at point.
            clon[i, :] = lon[i]
            clat[i, :] = lat[i]
            continue

        ex = lonc[elems]
        ey = latc[elems]

        # Order surrounding centroids by polar angle around the node so the
        # polygon does not self-intersect.
        ang = np.arctan2(ey - lat[i], ex - lon[i])
        order = np.argsort(ang)
        ex, ey = ex[order], ey[order]

        clon[i, : len(ex)] = ex
        clat[i, : len(ey)] = ey
        # pad by repeating the last real corner
        clon[i, len(ex):] = ex[-1]
        clat[i, len(ey):] = ey[-1]

    return clon, clat, max_corners


def build_elem_polygons(nv, lon, lat):
    """Triangle corners straight from connectivity -> shape (nele, 3)."""
    return lon[nv], lat[nv]


def polygon_area_deg(clon, clat):
    """
    Approximate spherical polygon area (m^2) via the shoelace formula on an
    equirectangular projection scaled by cos(lat) -- adequate for coupling
    area weights, not a substitute for exact spherical excess if you need
    high precision at very coarse/large cells.
    """
    R = 6371000.0
    lat0 = np.nanmean(clat, axis=1, keepdims=True)
    x = np.radians(clon) * R * np.cos(np.radians(lat0))
    y = np.radians(clat) * R

    x2 = np.roll(x, -1, axis=1)
    y2 = np.roll(y, -1, axis=1)
    area = 0.5 * np.abs(np.nansum(x * y2 - x2 * y, axis=1))
    return area


# ---------------------------------------------------------------------
# OASIS file writers (mirrors the grids.nc/areas.nc/masks.nc convention
# used by OASIS3-MCT / oasisgrids.py-esmgrids for structured models)
# ---------------------------------------------------------------------

def _open_append(path):
    return nc.Dataset(path, "a") if os.path.exists(path) else nc.Dataset(path, "w")


def _ensure_dim(f, name, size):
    if name not in f.dimensions:
        f.createDimension(name, size)


def _ensure_var(f, name, dtype, dims):
    if name in f.variables:
        return f.variables[name]
    return f.createVariable(name, dtype, dims)


def write_grid(grids_path, areas_path, masks_path, name, cell,
                lon, lat, clon, clat, area, mask, title):
    n = lon.shape[0]
    ny_dim = "ny{}_{}".format(cell, name)
    nx_dim = "nx{}_{}".format(cell, name)
    nc_dim = "nc{}_{}".format(cell, name)
    ncorners = clon.shape[1]

    lon_var = "{}{}.lon".format(name[:3], cell)
    lat_var = "{}{}.lat".format(name[:3], cell)
    clo_var = "{}{}.clo".format(name[:3], cell)
    cla_var = "{}{}.cla".format(name[:3], cell)
    srf_var = "{}{}.srf".format(name[:3], cell)
    msk_var = "{}{}.msk".format(name[:3], cell)

    with _open_append(grids_path) as f:
        _ensure_dim(f, ny_dim, 1)
        _ensure_dim(f, nx_dim, n)
        _ensure_dim(f, nc_dim, ncorners)

        v = _ensure_var(f, lat_var, "f8", (ny_dim, nx_dim))
        v.units = "degrees_north"
        v.title = "{} {}-point latitude".format(title, cell)
        v[:] = lat.reshape(1, n)

        v = _ensure_var(f, lon_var, "f8", (ny_dim, nx_dim))
        v.units = "degrees_east"
        v.title = "{} {}-point longitude".format(title, cell)
        v[:] = lon.reshape(1, n)

        v = _ensure_var(f, cla_var, "f8", (nc_dim, ny_dim, nx_dim))
        v.units = "degrees_north"
        v.title = "{} {}-point corner latitude".format(title, cell)
        v[:] = clat.T.reshape(ncorners, 1, n)

        v = _ensure_var(f, clo_var, "f8", (nc_dim, ny_dim, nx_dim))
        v.units = "degrees_east"
        v.title = "{} {}-point corner longitude".format(title, cell)
        v[:] = clon.T.reshape(ncorners, 1, n)

    with _open_append(areas_path) as f:
        _ensure_dim(f, ny_dim, 1)
        _ensure_dim(f, nx_dim, n)
        v = _ensure_var(f, srf_var, "f8", (ny_dim, nx_dim))
        v.units = "m^2"
        v.title = "{} {}-point area".format(title, cell)
        v[:] = area.reshape(1, n)

    with _open_append(masks_path) as f:
        _ensure_dim(f, ny_dim, 1)
        _ensure_dim(f, nx_dim, n)
        v = _ensure_var(f, msk_var, "i4", (ny_dim, nx_dim))
        v.units = "0/1:o/l"
        v.title = "{} {}-point land-sea mask (1=masked)".format(title, cell)
        v[:] = mask.reshape(1, n)


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("fvcom_file", help="FVCOM grid/output NetCDF file")
    p.add_argument("--grid_name", default="fvo",
                    help="OASIS grid id prefix (script uses first 3 chars). Default: fvo")
    p.add_argument("--grids", default="grids.nc")
    p.add_argument("--areas", default="areas.nc")
    p.add_argument("--masks", default="masks.nc")
    p.add_argument("--which", choices=["both", "node", "element"], default="both")
    p.add_argument("--mask_var", default=None,
                    help="Optional wet/dry variable (1=wet) to derive the mask from, "
                         "e.g. wet_nodes / wet_cells. Default: no masking.")
    args = p.parse_args()

    ds = nc.Dataset(args.fvcom_file)

    lon, lat, lonc, latc = pick_coords(ds)
    n_nodes = lon.shape[0]
    nele = ds.dimensions["nele"].size if "nele" in ds.dimensions else None

    nv = None
    if nele is not None:
        nv = get_nv(ds, nele)
        if lonc is None or latc is None:
            lonc, latc = element_centroids(nv, lon, lat)

    if args.which in ("node", "both"):
        if nv is None:
            raise ValueError("Need 'nv' connectivity (and nele dim) to build node dual-cell corners.")
        clon, clat, _ = build_node_polygons(nv, n_nodes, lon, lat, lonc, latc)

        if "art1" in ds.variables:
            area = ds.variables["art1"][:].astype("f8")
        else:
            area = polygon_area_deg(clon, clat)

        if args.mask_var and args.mask_var in ds.variables:
            wet = ds.variables[args.mask_var][:].astype("i4")
            mask = 1 - wet
        else:
            mask = np.zeros(n_nodes, dtype="i4")

        write_grid(args.grids, args.areas, args.masks, args.grid_name, "n",
                   lon, lat, clon, clat, area, mask, "FVCOM node grid")
        print("Wrote node grid '{}n' ({} points)".format(args.grid_name[:3], n_nodes))

    if args.which in ("element", "both"):
        if nv is None:
            raise ValueError("Need 'nv' connectivity (and nele dim) to build element grid.")
        eclon, eclat = build_elem_polygons(nv, lon, lat)

        if "art" in ds.variables:
            earea = ds.variables["art"][:].astype("f8")
        else:
            earea = polygon_area_deg(eclon, eclat)

        mask_var_e = None
        if args.mask_var and args.mask_var.replace("nodes", "cells") in ds.variables:
            mask_var_e = args.mask_var.replace("nodes", "cells")
        if mask_var_e:
            wet = ds.variables[mask_var_e][:].astype("i4")
            emask = 1 - wet
        else:
            emask = np.zeros(nele, dtype="i4")

        write_grid(args.grids, args.areas, args.masks, args.grid_name, "e",
                   lonc, latc, eclon, eclat, earea, emask, "FVCOM element grid")
        print("Wrote element grid '{}e' ({} points)".format(args.grid_name[:3], nele))

    ds.close()


if __name__ == "__main__":
    raise SystemExit(main())
