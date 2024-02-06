"""This script creates initial condition of 
   free run as prior for the following DA experiments
"""
import xarray as xr

it = 15
f = xr.open_dataset('free_for.nc', decode_times=False)
for i in range(16):
    f.isel(time=[it], ens=i).to_netcdf('initFree/maooam_{:03}.nc'.format(i+1))
f.close()