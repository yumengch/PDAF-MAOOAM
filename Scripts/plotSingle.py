import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgs
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors


_varnames = ['psi_a', 'T_a', 'psi_o', 'T_o']
_titles = [r'$\Psi_a$', r'$T_a$', r'$\Psi_o$', r'$T_o$']
_linestyles = ['solid', 'dashed', 'dotted']
_DAtime = 100000 + 15
_nx = 129
_ny = 129
_Ne = 16

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['font.size'] = 20


def plotErrorMovingStrong():
    burnin = 360
    years = [0, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    error = dict()
    for alpha in ['A1', 'O1', 'free']:
        error[alpha] = dict()
        for varname in _varnames:
            error[alpha][varname] = np.zeros(len(years))
        for i, nyear in enumerate(years):
            if nyear == 0:
                if alpha == 'free':
                    f = np.load('/home/users/yumengch/MAOOAM_EXPs/results/free_deflate/error.npz')
                    for varname in _varnames:
                        error[alpha][varname][i] = f[varname+'_a'].ravel()[burnin:].mean()
                else:
                    f = np.load(f'getError/error_all_{alpha}.npz')
                    for varname in _varnames:
                        error[alpha][varname][i] = f['err_'+varname+'_a'].ravel()[burnin:].mean()

            else:
                if alpha == 'free':
                    f = np.load(f'/home/users/yumengch/MAOOAM_EXPs/results/free_deflate/movingerror_all_{nyear*360}.npz')
                    for varname in _varnames:
                        error[alpha][varname][i] = f['err_'+varname].ravel()[burnin:].mean()
                else:
                    f = np.load(f'getError/movingerror_all_{nyear*360}_{alpha}.npz')
                    for varname in _varnames:
                        error[alpha][varname][i] = f['err_'+varname].ravel()[burnin:].mean()


    fig = plt.figure(1)
    fig.clf()
    w, h = fig.get_size_inches()
    fig.set_size_inches(2*w, 2*h)
    gs = mgs.GridSpec(2, 2, left=0.09, bottom=0.14, right=0.996, top=0.96, wspace=0.25, hspace=0.3)

    for i, (varname, title) in enumerate(zip(_varnames, _titles)):
        ax = fig.add_subplot(gs[i])
        lines = []
        line, = ax.plot(years, error['free'][varname], color='r', label='free')
        lines.append(line)
        line, = ax.plot(years, error['A1'][varname], color='k', linestyle='solid', label='ocean obs. only')
        lines.append(line)
        line, = ax.plot(years, error['O1'][varname], color='k', linestyle='dotted', label='atmos. obs. only')
        lines.append(line)
        ax.set_ylabel('Time averaged RMSE')
        ax.set_xlabel('Moving average window (years)')
        ax.set_yscale('log')
        ax.set_title(title)

    fig.legend(loc='lower center', handles=lines, ncol=3,fontsize=16)
    fig.savefig(f'err_moving.eps', dpi=300)


if __name__ == '__main__':
    plotErrorMovingStrong()
