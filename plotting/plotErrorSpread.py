import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt # type: ignore
import matplotlib.gridspec as mgs # type: ignore
import matplotlib.ticker as mticker # type: ignore
import matplotlib.colors as mcolors # type: ignore
# import mod_model # type: ignore
import time

from mpi4py import MPI


_varnames = ['psi_a', 'T_a', 'psi_o', 'T_o']
_titles = [r'$\Psi_a$', r'$T_a$', r'$\Psi_o$', r'$T_o$']
_obserr = [0.0169, 0.0076, 0.0011, 0.049]
_DAtime = 100000 + 15
_nx = 129
_ny = 129
_Ne = 16

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['font.size'] = 20

def mergeData(dirname):
    E = [xr.open_dataset('{:}/maooam_{:03}.nc'.format(dirname, i), decode_times=False) for i in range(1, 17)]
    E = [f.assign_coords({'ens': (['ens',], [i],
                                  {'long_name': 'ensemble index',
                                   'standard_name': 'ensemble number',
                                   'units': '1'
                                  }
                                 )
                         }) for i, f in enumerate(E)]
    return xr.combine_by_coords(E, combine_attrs='drop_conflicts')


def toPhysical(coeffs, fields):
    mod_model.mod_model.tophysical_a_multistep(coeffs, fields['psi_a'], fields['T_a'])
    mod_model.mod_model.tophysical_o_multistep(coeffs, fields['psi_o'], fields['T_o'])
    return fields


def getRMSE(truth, forecast):
    return np.sqrt(((forecast - truth)*(forecast - truth)).mean(axis=(-1)))


def saveError(is_obs, is_all):
    ft = xr.open_dataset('/storage/research/nceo/ia923171/pyPDAF_MAOOAM/results/truth/truth_for.nc', decode_times=False)
    f = xr.open_dataset('/storage/research/nceo/ia923171/MAOOAM/both/maooamA1O1.nc', decode_times=False)
    f = f.isel(time=range(f.dims['time']-1))
    ft = ft.isel(time=range(_DAtime, _DAtime+f.dims['time']))
    print (f['time'][0].to_numpy(), ft['time'][0].to_numpy())
    nt_total = ft.dims['time']

    nPE = MPI.COMM_WORLD.Get_size()
    # assert nPE == 8, 'set 8 processors used'
    iRank = MPI.COMM_WORLD.Get_rank()
    nt= nt_total//nPE
    if iRank == nPE - 1: nt = nt_total - (nPE-1)*nt
    n_subnt = 125
    assert np.isclose(nt//n_subnt, nt/n_subnt), f'try different n_subnt for {nt} steps'
    it_start = nt*iRank
    it_end = it_start + nt

    # allocate arrays
    error = dict()
    error_f = dict()
    fields = dict()
    fields_truth = dict()
    for varname in _varnames:
        fields[varname] = np.zeros((n_subnt*_Ne, _nx, _ny), order='F')
        fields_truth[varname] = np.zeros((n_subnt, _nx, _ny), order='F')
        error[varname] = np.zeros(nt)
        error_f[varname] = np.zeros(nt)

    mask = np.zeros((_nx, _ny), dtype=bool)
    mask[::8][:, ::8] = True
    if not is_obs: mask = ~mask
    if is_all: mask[:] = True
    mod_model.mod_model.initialize_model()
    start = time.time()
    for it in range(it_start, it_end, n_subnt):
        print (it, it+n_subnt)
        coeffs_t = np.concatenate([ft[varname+'_a'][it:it+n_subnt].to_numpy() for varname in _varnames], axis=-1)
        fields_truth = toPhysical(coeffs_t, fields_truth)

        coeffs = np.concatenate([f[varname+'_a'][:, it:it+n_subnt].to_numpy() for varname in _varnames], axis=-1)
        coeffs = coeffs.reshape(n_subnt*_Ne, 36)
        fields = toPhysical(coeffs, fields)

        for varname in _varnames:
            error[varname][it-it_start:it-it_start+n_subnt] = \
                getRMSE(fields_truth[varname][:, mask],
                        fields[varname].reshape(_Ne,n_subnt,_nx,_ny).mean(axis=0)[:, mask])

        coeffs = np.concatenate([f[varname+'_f'][:, it:it+n_subnt].to_numpy() for varname in _varnames], axis=-1)
        coeffs = coeffs.reshape(n_subnt*_Ne, 36)
        fields = toPhysical(coeffs, fields)

        for varname in _varnames:
            error_f[varname][it-it_start:it-it_start+n_subnt] = \
                getRMSE(fields_truth[varname][:, mask],
                        fields[varname].reshape(_Ne,n_subnt,_nx,_ny).mean(axis=0)[:, mask])

    print (time.time() - start)
    mod_model.mod_model.finalize_model()

    error_total = dict()
    error_total_f = dict()
    for varname in _varnames:
        error_total[varname+'_a'] = None
        error_total_f[varname+'_f'] = None
        if iRank == 0:
            error_total[varname+'_a'] = np.zeros((nPE, nt))
            error_total_f[varname+'_f'] = np.zeros((nPE, nt))

    for varname in _varnames:
        MPI.COMM_WORLD.Gather(error[varname], error_total[varname+'_a'], root=0)
        MPI.COMM_WORLD.Gather(error_f[varname], error_total_f[varname+'_f'], root=0)

    if iRank == 0:
        if is_obs:
            np.savez('error_obs.npz', **error_total, **error_total_f)
        else:
            np.savez('error_noobs.npz', **error_total, **error_total_f)
        if is_all: np.savez('error.npz', **error_total, **error_total_f)


def saveSpread(is_obs, is_all):
    f = xr.open_dataset('/storage/research/nceo/ia923171/MAOOAM/both/maooamA1O1.nc', decode_times=False)
    f = f.isel(time=range(f.dims['time']-1))
    nt_total = f.dims['time']

    nPE = MPI.COMM_WORLD.Get_size()
    # assert nPE == 8, 'set 8 processors used'
    iRank = MPI.COMM_WORLD.Get_rank()
    nt = nt_total//nPE
    if iRank == nPE - 1: nt = nt_total - (nPE-1)*nt
    n_subnt = 125
    assert np.isclose(nt//n_subnt, nt/n_subnt), f'try different n_subnt for {nt} steps'
    it_start = nt*iRank
    it_end = it_start + nt

    print ((iRank, it_start, it_end, n_subnt))
    # allocate arrays
    spread = dict()
    fields = dict()
    for varname in _varnames:
        fields[varname] = np.zeros((n_subnt*_Ne, _nx, _ny), order='F')
        spread[varname] = np.zeros(nt)

    mask = np.zeros((_nx, _ny), dtype=bool)
    mask[::8][:, ::8] = True
    if not is_obs: mask = ~mask
    if is_all: mask[:] = True
    mod_model.mod_model.initialize_model()
    start = time.time()
    for it in range(it_start, it_end, n_subnt):
        print ((iRank, it, it+n_subnt))
        coeffs = np.concatenate([f[varname+'_a'][:, it:it+n_subnt].to_numpy() for varname in _varnames], axis=-1)
        coeffs = coeffs.reshape(n_subnt*_Ne, 36)
        fields = toPhysical(coeffs, fields)
        for varname in _varnames:
            spread[varname][it-it_start:it-it_start+n_subnt] = fields[varname].reshape(_Ne, n_subnt, _nx,_ny).std(axis=0, ddof=1)[:, mask].mean(axis=-1)
    print (time.time() - start)
    mod_model.mod_model.finalize_model()

    spread_total = dict()
    for varname in _varnames:
        spread_total[varname] = None
        if iRank == 0:
            spread_total[varname] = np.zeros((nPE, nt))

    for varname in _varnames:
        MPI.COMM_WORLD.Gather(spread[varname], spread_total[varname], root=0)

    if iRank == 0:
        if is_obs:
            np.savez('spread_obs.npz', **spread_total)
        else:
            np.savez('spread_noobs.npz', **spread_total)
        if is_all: np.savez('spread.npz', **spread_total)


def getCorr(fields, nt, ana):
    E = np.zeros((_Ne, 2*_nx*_ny))
    for it in range(nt):
        fig = plt.figure()
        gs = mgs.GridSpec(2, 2)
        fig.clf()
        for i, varname1 in enumerate(['psi_a', 'T_a']):
            for j, varname2 in enumerate(['psi_o', 'T_o']):
                E[:, :_nx*_ny] = fields[varname1][:, it].reshape(_Ne, _nx*_ny)
                E[:, _nx*_ny:] = fields[varname2][:, it].reshape(_Ne, _nx*_ny)
                A = E - E.mean(axis=0, keepdims=True)
                P = A.T@A/(_Ne-1)
                Ps = np.diag(P[:_nx*_ny, _nx*_ny:])
                P = np.diag(P[_nx*_ny:, :_nx*_ny])
                print (np.max(np.abs(P -Ps)))
                corr = P/np.std(E[:, :_nx*_ny], axis=0, ddof=1)/np.std(E[:, _nx*_ny:], axis=0, ddof=1)
                ax = fig.add_subplot(gs[i, j])
                ax.set_title(f'{varname1} .vs. {varname2}')
                c = ax.pcolormesh(corr.reshape(_nx, _ny), cmap='RdBu', norm=mcolors.CenteredNorm())
                fig.colorbar(c, ax=ax)
        fig.savefig(f'corr_{it}_{ana}.png', dpi=300)
        plt.close(fig)


def plotCorr():
    f = xr.open_dataset('maooam_A1O7_For.nc', decode_times=False)
    f = f.isel(time=range(2000))
    nt = f.dims['time']

    # allocate arrays
    fields = dict()
    for varname in _varnames:
        fields[varname] = np.zeros((nt*_Ne, _nx, _ny), order='F')

    mod_model.mod_model.initialize_model()
    start = time.time()

    coeffs = np.concatenate([f[varname+'_a'].to_numpy() for varname in _varnames], axis=-1)
    coeffs = coeffs.reshape(nt*_Ne, 36)
    fields = toPhysical(coeffs, fields)
    for varname in _varnames:
        fields[varname] = fields[varname].reshape(_Ne, nt, _nx, _ny)
    getCorr(fields, nt, 'a')

    coeffs = np.concatenate([f[varname+'_f'].to_numpy() for varname in _varnames], axis=-1)
    coeffs = coeffs.reshape(nt*_Ne, 36)
    fields = toPhysical(coeffs, fields)
    for varname in _varnames:
        fields[varname] = fields[varname].reshape(_Ne, nt, _nx, _ny)
    getCorr(fields, nt, 'f')

    print (time.time() - start)
    mod_model.mod_model.finalize_model()


def plotIncrementInnovationError():
    ft = xr.open_dataset('/home/users/yumengch/MAOOAM_EXPs/results/truth/truth_for.nc', decode_times=False)
    f = xr.open_dataset('maooam_A1O7_For.nc', decode_times=False)
    # f = xr.open_dataset('../../results/free_deflate/free_for.nc', decode_times=False)
    f = f.isel(time=range(2000))
    ft = ft.isel(time=range(_DAtime, _DAtime+f.dims['time']))
    print (f.dims['time'], ft.dims['time'])
    print (f['time'][0].to_numpy(), ft['time'][0].to_numpy())
    nt_total = ft.dims['time']

    nPE = MPI.COMM_WORLD.Get_size()
    # assert nPE == 8, 'set 8 processors used'
    iRank = MPI.COMM_WORLD.Get_rank()
    nt= nt_total//nPE
    if iRank == nPE - 1: nt = nt_total - (nPE-1)*nt
    n_subnt = 50
    # n_subnt = 1
    assert np.isclose(nt//n_subnt, nt/n_subnt), f'try different n_subnt for {nt} steps'
    it_start = nt*iRank
    it_end = it_start + nt

    # allocate arrays
    error = dict()
    incr = dict()
    fields_f = dict()
    fields_a = dict()
    fields_truth = dict()
    for varname in _varnames:
        fields_f[varname] = np.zeros((n_subnt*_Ne, _nx, _ny), order='F')
        fields_a[varname] = np.zeros((n_subnt*_Ne, _nx, _ny), order='F')
        fields_truth[varname] = np.zeros((n_subnt, _nx, _ny), order='F')
        error[varname] = np.zeros((n_subnt, _nx, _ny))
        incr[varname] = np.zeros((n_subnt, _nx, _ny))

    mod_model.mod_model.initialize_model()
    start = time.time()
    for it in range(it_start, it_end, n_subnt):
        print (it, it+n_subnt)
        coeffs_t = np.concatenate([ft[varname+'_a'][it:it+n_subnt].to_numpy() for varname in _varnames], axis=-1)
        fields_truth = toPhysical(coeffs_t, fields_truth)

        coeffs = np.concatenate([f[varname+'_a'][:, it:it+n_subnt].to_numpy() for varname in _varnames], axis=-1)
        coeffs = coeffs.reshape(n_subnt*_Ne, 36)
        fields_f = toPhysical(coeffs, fields_f)

        for varname in _varnames:
            error[varname] = fields_truth[varname] - fields_f[varname].reshape(_Ne,n_subnt,_nx,_ny).mean(axis=0)


        coeffs = np.concatenate([f[varname+'_f'][:, it:it+n_subnt].to_numpy() for varname in _varnames], axis=-1)
        coeffs = coeffs.reshape(n_subnt*_Ne, 36)
        fields_a = toPhysical(coeffs, fields_a)

        for varname in _varnames:
            incr[varname] = \
                fields_a[varname].reshape(_Ne,n_subnt,_nx,_ny).mean(axis=0) - \
                    fields_f[varname].reshape(_Ne,n_subnt,_nx,_ny).mean(axis=0)

        fig = plt.figure(1)
        gs = mgs.GridSpec(4, 2)
        for i in range(n_subnt):
            fig.clf()
            for j, varname in enumerate(_varnames):
                ax = fig.add_subplot(gs[j, 0])
                c = ax.pcolormesh(error[varname][i], cmap='RdBu', norm=mcolors.CenteredNorm())
                ax.set_title(f'${varname}$' + '_err')
                fig.colorbar(c, ax=ax)
                ax = fig.add_subplot(gs[j, 1])
                c = ax.pcolormesh(incr[varname][i], cmap='RdBu', norm=mcolors.CenteredNorm())
                ax.set_title(f'${varname}$' + '_incr')
                fig.colorbar(c, ax=ax)
            fig.savefig(f'err_{it_start+i}.png', dpi=300)


    print (time.time() - start)
    mod_model.mod_model.finalize_model()


def plotError():
    """Plot the error and spread of the free run and the SCDA run.
    """
    t = np.arange(100000)*90*1000/1.032/3600/24/365
    fig = plt.figure(1)
    fig.clf()
    w, h = fig.get_size_inches()
    fig.set_size_inches(4*w, 2*h)
    gs = mgs.GridSpec(2, 4, left=0.04, bottom=0.14, right=0.995, top=0.92, wspace=0.25, hspace=0.54)

    errorObs = np.load('free/error_obs.npz')
    errorNoObs = np.load('free/error_noobs.npz')
    for varname in _varnames:
        print ('analysis error of observed quantities', varname, np.mean(errorObs[varname+'_a'].ravel()[365:]))
        print ('analysis error of unobserved quantities', varname, np.mean(errorNoObs[varname+'_a'].ravel()[365:]))
    errorAll = np.load('free/error.npz')
    spreadAll = np.load('free/spread.npz')

    for i, (varname, title, obserr) in enumerate(zip(_varnames, _titles, _obserr)):
        ax = fig.add_subplot(gs[0, i])
        LineErrAll, = ax.plot(t[365:], errorAll[varname+'_a'].ravel()[365:], color='k', label='err.')
        LineStdAll, = ax.plot(t[365:], spreadAll[varname].ravel()[365:], color='r', label='std dev.')
        lineObs = ax.axhline(obserr, color='grey', label='obs. std')
        ax.set_xlabel('year')
        ax.set_title(f'{title}')
    fig.text(0.5, 0.97, r"Time series of free run on $129 \times 129$ grid points", ha='center', fontsize=20)

    errorObs = np.load('scda/error_obs.npz')
    errorNoObs = np.load('scda/error_noobs.npz')
    for varname in _varnames:
        print ('analysis error of observed quantities', varname, np.mean(errorObs[varname+'_a'].ravel()[365:]))
        print ('analysis error of unobserved quantities', varname, np.mean(errorNoObs[varname+'_a'].ravel()[365:]))
    errorAll = np.load('scda/error.npz')
    spreadAll = np.load('scda/spread.npz')
    for i, (varname, title, obserr) in enumerate(zip(_varnames, _titles, _obserr)):
        ax = fig.add_subplot(gs[1, i])
        LineErrAll, = ax.plot(t[365:], errorAll[varname+'_a'].ravel()[365:], color='k', label='err.')
        LineStdAll, = ax.plot(t[365:], spreadAll[varname].ravel()[365:], color='r', label='std dev.')
        # lineObs = ax.axhline(obserr, color='grey', label='obs. std')
        ax.set_xlabel('year')
        ax.set_title(f'{title}')
    fig.text(0.5, 0.49, r"Time series of SCDA analysis on $129 \times 129$ grid points", ha='center', fontsize=20)

    fig.legend(loc='outside lower center',
               handles=[LineErrAll, LineStdAll], #, lineObs],
               ncols=3)
    fig.savefig('error.pdf', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    # saveError(is_obs=True, is_all=False)
    # saveError(is_obs=False, is_all=False)
    # saveError(is_obs=False, is_all=True)
    # saveSpread(is_obs=True, is_all=False)
    # saveSpread(is_obs=False, is_all=False)
    # saveSpread(is_obs=False, is_all=True)
    plotError()
