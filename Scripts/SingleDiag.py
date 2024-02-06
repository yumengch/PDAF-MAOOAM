import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgs
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import time
from mpi4py import MPI
import mod_model 


_varnames = ['psi_a', 'T_a', 'psi_o', 'T_o']
_titles = [r'$\Psi_a$', r'$T_a$', r'$\Psi_o$', r'$T_o$']
_DAtime = 100000 + 15
_nx = 129
_ny = 129
_Ne = 16

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['font.size'] = 20


def toPhysical(coeffs, fields, suffix):
    mod_model.mod_model.tophysical_a_multistep(coeffs, fields['psi_a'+suffix], fields['T_a'+suffix])
    mod_model.mod_model.tophysical_o_multistep(coeffs, fields['psi_o'+suffix], fields['T_o'+suffix])
    return fields


def getRMSE(truth, forecast):
    return np.sqrt(((forecast - truth)*(forecast - truth)).mean(axis=(-1)))


def CPUdistribution(nt_total):
    nPE = MPI.COMM_WORLD.Get_size()
    iRank = MPI.COMM_WORLD.Get_rank()
    nt= nt_total//nPE
    if iRank == nPE - 1: nt = nt_total - (nPE-1)*nt

    if nt < 100:
        n_subnt = 1
    else:
        n_subnt = 100
        while True:
            n_subnt -= 1
            if np.isclose(nt//n_subnt, nt/n_subnt):
                break
    print (iRank, nt, n_subnt)
    assert np.isclose(nt//n_subnt, nt/n_subnt), f'try different n_subnt for {nt} steps'
    tslices = [slice(nt*iRank + i*n_subnt, nt*iRank + (i+1)*n_subnt) for i in range(nt//n_subnt)]
    return nPE, iRank, nt, tslices, n_subnt


def initialiseFields(nt, n_subnt):
    # allocate arrays
    error = dict()
    fieldsSpace = dict()
    fields = dict()
    for varname in _varnames:
        fieldsSpace[varname+'_t'] = np.zeros((n_subnt, _nx, _ny), order='F')
        fields[varname+'_t'] = np.zeros(nt)

        for suffix in ['_a', '_f']:
            fieldsSpace[varname+suffix] = np.zeros((n_subnt, _nx, _ny), order='F')
            fields[varname+suffix] = np.zeros(nt)
            error[varname+suffix] = np.zeros(nt)

    return error, fields, fieldsSpace


def saveError(gridpoints, alpha):
    ft = xr.open_dataset('/home/users/yumengch/MAOOAM_EXPs/results/truth/truth_for.nc', decode_times=False)
    f = xr.open_dataset(f'/home/users/yumengch/MAOOAM_EXPs/For/AssmSingle/netCDFFiles/maooam{alpha}.nc', decode_times=False)
    f = f.isel(time=range(f.dims['time']-1))
    ft = ft.isel(time=range(_DAtime, _DAtime+f.dims['time']))
    print (f['time'][0].to_numpy(), ft['time'][0].to_numpy())
    nt_total = ft.dims['time']

    nPE, iRank, nt, tslices, n_subnt = CPUdistribution(nt_total)
    error, fields, fieldsSpace = initialiseFields(nt, n_subnt)

    mask = np.zeros((_nx, _ny), dtype=bool)
    if gridpoints == 'obs':
        mask[::8][:, ::8] = True
    elif gridpoints == 'unobs':
        mask[::8][:, ::8] = True
        mask = ~mask
    else:
        mask[:] = True

    mod_model.mod_model.initialize_model()
    start = time.time()

    for i, tslice in enumerate(tslices):
        fslice = slice(i*n_subnt, (i+1)*n_subnt)
        print ((iRank, tslice, nt))
        coeffs_t = np.concatenate([ft[varname+'_a'][tslice].to_numpy() for varname in _varnames], axis=-1)
        fieldsSpace = toPhysical(coeffs_t, fieldsSpace, '_t')

        for suffix in ['_a', '_f']:
            coeffs = np.concatenate([f[varname+suffix][:, tslice].to_numpy() for varname in _varnames], axis=-1)
            coeffs = coeffs.mean(axis=0)
            print ((iRank, suffix, coeffs.shape))
            fieldsSpace = toPhysical(coeffs, fieldsSpace, suffix)

            for varname in _varnames:
                error[varname+suffix][fslice] = \
                    getRMSE(fieldsSpace[varname+'_t'][:, mask],
                            fieldsSpace[varname+suffix][:, mask])
            for varname in _varnames:
                fields[varname+suffix][fslice] = fieldsSpace[varname+suffix].mean(axis=(-1, -2))
        for varname in _varnames:
            fields[varname+'_t'][fslice] = fieldsSpace[varname+'_t'].mean(axis=(-1, -2))

    print (time.time() - start)
    mod_model.mod_model.finalize_model()

    print ((iRank, 'gathering'))
    error_total = dict()
    fields_total = dict()
    for varname in _varnames:
        fields_total['f_'+varname+'_t'] = None
        for suffix in ['_a', '_f']:
            error_total['err_'+varname+suffix] = None
            fields_total['f_'+varname+suffix] = None
            if iRank == 0:
                error_total['err_'+varname+suffix] = np.zeros((nPE, nt))
                fields_total['f_'+varname+suffix] = np.zeros((nPE, nt))
        if iRank == 0:
            fields_total['f_'+varname+'_t'] = np.zeros((nPE, nt))

    for varname in _varnames:
        MPI.COMM_WORLD.Gather(fields[varname+'_t'], fields_total['f_'+varname+'_t'], root=0)
        for suffix in ['_a', '_f']:
            MPI.COMM_WORLD.Gather(error[varname+suffix], error_total['err_'+varname+suffix], root=0)
            MPI.COMM_WORLD.Gather(fields[varname+suffix], fields_total['f_'+varname+suffix], root=0)

    print ((iRank, 'saving'))
    if iRank == 0:
        np.savez(f'error_{gridpoints}_{alpha}.npz', **error_total, **fields_total) 


def saveSpread(gridpoints, alpha):
    f = xr.open_dataset(f'/home/users/yumengch/MAOOAM_EXPs/For/AssmStrongAlpha/netCDFFiles/maooam{alpha}.nc', decode_times=False)
    f = f.isel(time=range(f.dims['time']-1))
    nt_total = f.dims['time']

    nPE = MPI.COMM_WORLD.Get_size()
    iRank = MPI.COMM_WORLD.Get_rank()
    nt = nt_total//nPE
    if iRank == nPE - 1: nt = nt_total - (nPE-1)*nt
    n_subnt = 100
    assert np.isclose(nt//n_subnt, nt/n_subnt), f'try different n_subnt for {nt} steps'
    it_start = nt*iRank
    it_end = it_start + nt

    # allocate arrays
    spread = dict()
    fields = dict()
    for varname in _varnames:
        fields[varname] = np.zeros((n_subnt*_Ne, _nx, _ny), order='F')
        for suffix in ['_a', '_f']:
            spread[varname+suffix] = np.zeros(nt)

    mask = np.zeros((_nx, _ny), dtype=bool)
    if gridpoints == 'obs':
        mask[::8][:, ::8] = True
    elif gridpoints == 'unobs':
        mask[::8][:, ::8] = True
        mask = ~mask
    else:
        mask[:] = True

    mod_model.mod_model.initialize_model()
    start = time.time()
    for it in range(it_start, it_end, n_subnt):
        print ((iRank, it, it+n_subnt))
        for suffix in ['_a', '_f']:
            coeffs = np.concatenate([f[varname+'_a'][:, it:it+n_subnt].to_numpy()
                                     for varname in _varnames], axis=-1)
            coeffs = coeffs.reshape(n_subnt*_Ne, 36)
            fields = toPhysical(coeffs, fields, '')
            for varname in _varnames:
                spread[varname+suffix][it-it_start:it-it_start+n_subnt] = \
                fields[varname].reshape(_Ne, n_subnt, _nx,_ny).std(axis=0, ddof=1)[:, mask].mean(axis=-1)
    print (time.time() - start)
    mod_model.mod_model.finalize_model()

    spread_total = dict()
    for varname in _varnames:
        for suffix in ['_a', '_f']:
            spread_total[varname+suffix] = None
            if iRank == 0:
                spread_total[varname+suffix] = np.zeros((nPE, nt))

    for varname in _varnames:
        for suffix in ['_a', '_f']:
            MPI.COMM_WORLD.Gather(spread[varname+suffix], spread_total[varname+suffix], root=0)

    if iRank == 0:
        np.savez(f'spread_{gridpoints}_{alpha}.npz', **spread_total)


def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0, axis=0), axis=0) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)


def saveErrorMoving(gridpoints, window, alpha):
    ft = xr.open_dataset('/home/users/yumengch/MAOOAM_EXPs/results/truth/truth_for.nc', decode_times=False)
    f = xr.open_dataset(f'/home/users/yumengch/MAOOAM_EXPs/For/AssmSingle/netCDFFiles/maooam{alpha}.nc', decode_times=False)
    f = f.isel(time=range(f.dims['time']-1))
    ft = ft.isel(time=range(_DAtime, _DAtime+f.dims['time']))
    print (f['time'][0].to_numpy(), ft['time'][0].to_numpy())

    coeffs_t = np.concatenate([ft[varname+'_a'].to_numpy() for varname in _varnames], axis=-1)
    coeffs_t = running_mean(coeffs_t, window)
    coeffs = np.concatenate([f[varname+'_a'].to_numpy() for varname in _varnames], axis=-1)
    coeffs = coeffs.mean(axis=0)
    coeffs = running_mean(coeffs, window)
    print (coeffs.shape)
    nt_total = len(coeffs)
    f.close()
    ft.close()

    nPE, iRank, nt, tslices, n_subnt = CPUdistribution(nt_total)
    coeffs_t = coeffs_t[tslices[0].start:tslices[-1].stop]
    coeffs = coeffs[tslices[0].start:tslices[-1].stop]

    # allocate arrays
    error = dict()
    fields = dict()
    fieldsSpace = dict()
    for varname in _varnames:
        fieldsSpace[varname] = np.zeros((n_subnt, _nx, _ny), order='F')
        fieldsSpace[varname+'_t'] = np.zeros((n_subnt, _nx, _ny), order='F')
        error[varname] = np.zeros(nt)
        fields[varname] = np.zeros(nt)
        fields[varname+'_t'] = np.zeros(nt)

    mask = np.zeros((_nx, _ny), dtype=bool)
    if gridpoints == 'obs':
        mask[::8][:, ::8] = True
    elif gridpoints == 'unobs':
        mask[::8][:, ::8] = True
        mask = ~mask
    else:
        mask[:] = True

    mod_model.mod_model.initialize_model()
    start = time.time()
    for i, tslice in enumerate(tslices):
        print ((iRank, tslice))
        fslice = slice(i*n_subnt, (i+1)*n_subnt)
        coeffs_tmp = coeffs[fslice]
        coeffs_t_tmp = coeffs_t[fslice]
        print (iRank, tslice, nt)
        fieldsSpace = toPhysical(coeffs_t_tmp, fieldsSpace, '_t')
        fieldsSpace = toPhysical(coeffs_tmp, fieldsSpace, '')
        for varname in _varnames:
            error[varname][fslice] = \
                getRMSE(fieldsSpace[varname+'_t'][:, mask],
                        fieldsSpace[varname][:, mask])
            fields[varname][fslice] = fieldsSpace[varname].mean(axis=(-1, -2))
            fields[varname+'_t'][fslice] = fieldsSpace[varname+'_t'].mean(axis=(-1, -2))

    print (time.time() - start)
    mod_model.mod_model.finalize_model()

    error_total = dict()
    fields_total = dict()
    for varname in _varnames:
        fields_total['f_'+varname+'_t'] = None
        error_total['err_'+varname] = None
        fields_total['f_'+varname] = None
        if iRank == 0:
            error_total['err_'+varname] = np.zeros((nPE, nt))
            fields_total['f_'+varname] = np.zeros((nPE, nt))
            fields_total['f_'+varname+'_t'] = np.zeros((nPE, nt))

    for varname in _varnames:
        MPI.COMM_WORLD.Gather(fields[varname+'_t'], fields_total['f_'+varname+'_t'], root=0)
        MPI.COMM_WORLD.Gather(error[varname], error_total['err_'+varname], root=0)
        MPI.COMM_WORLD.Gather(fields[varname], fields_total['f_'+varname], root=0)

    if iRank == 0:
        np.savez(f'movingerror_{gridpoints}_{window}_{alpha}.npz', **error_total, **fields_total) 


if __name__ == '__main__':
    import sys
    saveError('all', 'A1')
    saveError('all', 'O1')
    for nyear in [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
        saveErrorMoving('all', nyear*360, 'A1')
        saveErrorMoving('all', nyear*360, 'O1')
