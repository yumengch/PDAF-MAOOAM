"""This script creates initial condition of 
   truth used for free run, and used for obs. generation.
"""
import xarray as xr

_Freeinit = 100000
f = xr.open_dataset('truth_for.nc', decode_times=False)
f = f.isel(time=[_DAinit]).to_netcdf('free_init.nc')
f.close()

_DAinit = 100000 + 15
f = xr.open_dataset('truth_for.nc', decode_times=False)
f = f.isel(time=[_DAinit]).to_netcdf('DA_init.nc')
f.close()