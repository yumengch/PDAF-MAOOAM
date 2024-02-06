import xarray as xr
import numpy as np

f = xr.open_dataset('var-129.nc')
i = np.arange(129)
i = i[::8]
f = f.isel(nx=i)
f = f.isel(ny=i)
f.to_netcdf('traj_var.nc')
