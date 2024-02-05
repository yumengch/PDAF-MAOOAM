from MAOOAM import Model
from MAOOAM.parallelization import parallelization
from MAOOAM.model.params_maooam import f0, dt, natm, noc

import numpy as np
import xarray as xr

import datetime

def getFileAttrs():
    attrs = dict()
    attrs['Conventions'] = 'CF-1.8'
    attrs['title'] = 'NetCDF output from MAOOAM-pyPDAF'
    attrs['institution'] = 'NCEO-AWI-UoR'
    attrs['source'] = 'pyPDAF and MAOOAM'
    attrs['history'] = f'{datetime.datetime.now().isoformat(timespec="seconds")}: Data created'
    attrs['reference'] = 'https://github.com/yumengch/pyPDAF'
    return attrs


def saveData(field):
    data = dict()
    attrs = dict()

    varnames = ['psi_a', 'T_a', 'psi_o', 'T_o']
    dims = ['natm', 'natm', 'noc', 'noc']
    offsets = [0, natm, 2*natm, 2*natm + noc, 2*(natm + noc)]

    for i, (varname, dim) in enumerate(zip(varnames, dims)):
        attrs['standard_name'] = 'atmosphere_' if 'a' in varname else 'ocean_'
        attrs['standard_name'] += 'streamfunction_' if 'psi' in varname else 'temperature_'
        attrs['standard_name'] += 'coefficient'

        s = attrs['standard_name'].split('_')
        attrs['long_name'] = f'coefficient of {s[1]} in {s[0]}'
        data[varname] = xr.DataArray(field[:, offsets[i]:offsets[i+1]],
                                     dims=['time', dim], name=varname, attrs=attrs)
    return data


def genTrajectory(T_spinup, T, save_freq):
    pe = parallelization(dim_ens=1, n_modeltasks=1, screen=2)
    nt_spinup = int(f0*T_spinup/dt)
    nt = int(f0*T/dt)
    n_save_freq = int(save_freq*f0/dt)

    print (f'... The spin-up time steps is {nt_spinup}...')
    print (f'... The runtime steps is {nt}...')
    print (f'... The runtime field is saved every {n_save_freq} steps...')

    model = Model.Model(t0=0, nt=nt_spinup, pe=pe)
    model.init_field()
    for step in range(nt_spinup):
        model.step(pe, step, True)

    field_p = model.field_p
    model = Model.Model(t0=nt_spinup*dt, nt=nt, pe=pe)
    model.field_p = field_p
    field = np.zeros((nt//n_save_freq + 1, model.nx))
    cnt = 0
    for step in range(nt):
        model.step(pe, step, True)
        if step % n_save_freq == 0:
            print(f"...{step/nt*100}%...")
            field[cnt] = model.field_p.copy()
            cnt = cnt + 1

    return field, nt_spinup*dt/f0, nt*dt/f0, n_save_freq*dt/f0


if __name__ == '__main__':
    yearInSecs = 3600*24*365
    monthInSecs = 24*3600*30
    nyear_spin = 0.01
    T_spinup = nyear_spin*yearInSecs
    nyear = 0.1
    T = nyear*yearInSecs
    model_field, start_time, model_time, save_freq = genTrajectory(T_spinup, T, save_freq=monthInSecs/30./30.)
    time = np.arange(start_time, start_time + model_time, save_freq)
    assert len(time) == len(model_field)
    ds = xr.Dataset(saveData(model_field),
                    coords={'time' : ('time', time)},
                    attrs=getFileAttrs())
    ds.to_netcdf('trajectory.nc')