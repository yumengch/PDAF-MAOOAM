"""selecting truth model trajectory used for covariance matrix generation
"""
import xarray as xr
f = xr.open_dataset('traj.nc',decode_times=False)
f0 = 1.032*1e-4
dt = 1/f0
dy = 360*24*3600/dt
f = f.sel(time=slice(1000*dy, 1100*dy))
nt = f.dims['time']
f = f.isel(time=range(0, nt, 10))
print (f.dims['time'])
f.to_netcdf('trajectory.nc')
