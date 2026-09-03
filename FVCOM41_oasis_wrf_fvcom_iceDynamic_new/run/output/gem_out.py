import numpy as np
import netCDF4 as nc

# 1. Load the FVCOM netCDF file
file_path = "gl_0001.nc"
dataset = nc.Dataset(file_path, "r")

# 2. Extract your target variable (e.g., sea surface height 'zeta' or velocity 'u')
for var in ["zeta", "u", "v", "temp", "salinity"]:
# 'var' loads as a masked array or raw numpy array
    var_data = dataset.variables[var][:]

    # 3. Clean the data by removing NaNs and fill values
    # This flattens the array and drops NaNs entirely before calculation
    clean_data = var_data[~np.isnan(var_data)]
    print(f"({var}): non-nan shape = {clean_data.shape} of {var_data.flatten().shape}")
    # 4. Find the extreme values safely
    if clean_data.size > 0:
        max_val = np.max(clean_data)
        min_val = np.min(clean_data)
        print(f"({var}): Maximum Value (Before NaN): {max_val}")
        print(f"({var}): Minimum Value (Before NaN): {min_val}")
    else:
        print("({var}): Dataset contains only NaN values.")

dataset.close()

