import netCDF4 as nc
import numpy as np

d = nc.Dataset('gl_0001.nc')
last = -1  # last time record

for name, lo, hi in [('zeta', -30, 30), ('u', -30, 30), ('v', -30, 30),
                      ('temp', -10, 50), ('salinity', -2, 50)]:
    if name in d.variables:
        arr = d.variables[name][last]
        bad = np.where((arr < lo) | (arr > hi))
        if bad[0].size > 0:
            print(name, "extreme values at indices:", bad, "values:", arr[bad])
