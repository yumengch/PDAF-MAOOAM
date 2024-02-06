"""generation of initial condition of the truth
"""
import numpy as np
import xarray as xr

varnames = ['psi_a_f', 'T_a_f', 'psi_o_f', 'T_o_f']
f = xr.open_dataset('/home/users/yumengch/MAOOAM_EXPs/results/truth/trajectory.nc', decode_times=False).isel(time=[0])
f = f[varnames]
for varname in varnames:
	f[varname][:] = np.random.random(len(f[varname][0]))*0.01
	print (f[varname][:])
f.to_netcdf('init.nc')
