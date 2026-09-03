import netCDF4, numpy as np

f   = "gl_restart_0001.nc"
new = "1997-01-01T09:00:00.000000"

ds  = netCDF4.Dataset(f, "r+")

T      = ds.variables["Times"]
strlen = T.shape[-1]                       # DateStrLen (e.g. 26)
T[0, :] = " "
for i, c in enumerate(new):
    T[0, i] = c

ds.variables["time"][0]   = 50449.375
ds.variables["Itime"][0]  = 50449
ds.variables["Itime2"][0] = 32400000
ds.close()
