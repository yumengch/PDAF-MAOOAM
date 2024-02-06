import os
import numpy as np
import cmocean
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgs
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import matplotlib.colorbar as mcolorbar


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

def plotBar():
    burnin = 360
    fo = np.load('getError/error_obs_A0O0.npz')
    fuo = np.load('getError/error_unobs_A0O0.npz')
    fall = np.load('getError/error_all_A0O0.npz')

    foFree = np.load('/home/users/yumengch/MAOOAM_EXPs/results/free_deflate/error_obs.npz')
    fuoFree = np.load('/home/users/yumengch/MAOOAM_EXPs/results/free_deflate/error_noobs.npz')
    fallFree = np.load('/home/users/yumengch/MAOOAM_EXPs/results/free_deflate/error.npz')

    analysis = dict()
    free = dict()
    for t in ['o', 'uo', 'all']:
        analysis[t] = np.zeros(4)
        free[t] = np.zeros(4)

    for i, varname in enumerate(_varnames):
        analysis['o'][i] = fo['err_'+varname+'_a'].ravel(
                        )[burnin:].mean(axis=0)
        analysis['uo'][i] = fuo['err_'+varname+'_a'].ravel(
                        )[burnin:].mean(axis=0)
        analysis['all'][i] = fall['err_'+varname+'_a'].ravel(
                        )[burnin:].mean(axis=0)

        free['o'][i] = foFree[varname+'_a'].ravel()[burnin:].mean(axis=0)
        free['uo'][i] = fuoFree[varname+'_a'].ravel()[burnin:].mean(axis=0)
        free['all'][i] = fallFree[varname+'_a'].ravel()[burnin:].mean(axis=0)

    fig = plt.figure(1)
    fig.clf()
    w, h = fig.get_size_inches()
    fig.set_size_inches(2*w, 1*h)
    gs = mgs.GridSpec(1, 2, left=0.07, bottom=0.25, right=0.98, top=0.99, wspace=0.35, hspace=0.2)

    x = np.arange(0, 4*2, 2)  # the label locations
    width = 0.25  # the width of the bars
    ax = fig.add_subplot(gs[0])
    rects = []
    for i, (t, hatch, label) in enumerate(zip(['o', 'uo'], ['', 'xxx'], ['obs.', 'no obs.'])):
        offset = 2*i*width
        rect = ax.bar(x + offset, analysis[t], width, color='white',
                      edgecolor='k', hatch=hatch, label=f'analysis RMSE ({label})')
        rects.append(rect)
        offset = (2*i+1)*width
        rect = ax.bar(x + offset, free[t], width, color='gray',
                      edgecolor='k', hatch=hatch, label=f'freerun RMSE ({label})')
        rects.append(rect)
    ax.set_xticks(x + width, _titles)
    plt.yscale('log')
    fig.legend(loc=(0.001, 0.01), handles=rects, ncol=2,fontsize=16)

    fsall = np.load('getError/error_all_A1O1.npz')
    analysisS = np.zeros(4)

    for i, varname in enumerate(_varnames):
        analysisS[i] = fsall['err_'+varname+'_a'].ravel(
                        )[burnin:].mean(axis=0)

    ax = fig.add_subplot(gs[1])
    rects = []
    rect = ax.bar(x, analysis['all'], width, color='white',
                  edgecolor='k', hatch='', label=f'WCDA RMSE')
    rects.append(rect)
    rect = ax.bar(x + width, analysisS, width, color='gray',
                  edgecolor='k', hatch='', label=f'SCDA RMSE')
    rects.append(rect)
    ax.set_xticks(x + width, _titles)
    plt.yscale('log')
    fig.legend(loc=(0.6, 0.01), handles=rects, ncol=2,fontsize=16)

    fig.savefig('WS_bar.eps', dpi=300)
    plt.close(fig)


def plotCorr():
    f = np.load('getError/corr_weak_mid.npz')

    x = np.linspace(0, np.pi*2/1.5, 129)
    y = np.linspace(0, np.pi, 129)
    locator_x = mticker.FixedLocator([0, np.pi/1.5, np.pi*2/1.5])
    formatter_x = mticker.FixedFormatter(['0', r'$\frac{\pi}{n}$', r'$\frac{2\pi}{n}$'])
    locator_y = mticker.FixedLocator([0, np.pi*.5, np.pi])
    formatter_y = mticker.FixedFormatter(['0', r'$\frac{\pi}{2}$', r'$\pi$'])

    fig = plt.figure(1)
    fig.clf()
    w, h = fig.get_size_inches()
    fig.set_size_inches(2*w, 2*h)
    gs = mgs.GridSpec(2, 2, left=0.025, bottom=0.12, right=0.98, top=0.96, wspace=0.13, hspace=0.22)
    for i, (varname, title) in enumerate(zip(_varnames, _titles)):
        ax = fig.add_subplot(gs[i])
        m = ax.pcolormesh(x, y, f[varname], cmap='cmo.balance', vmin=-0.1, vmax=0.1)

        ax.xaxis.set_major_locator(locator_x)
        ax.xaxis.set_major_formatter(formatter_x)
        ax.yaxis.set_major_locator(locator_y)
        ax.yaxis.set_major_formatter(formatter_y)

        ax.set_title('Corr('+title+')')

    ax = fig.add_axes((0.25, 0.04, 0.5, 0.02))
    mcolorbar.Colorbar(ax, m, orientation='horizontal')
    fig.savefig('corr_weak.eps', dpi=300)
    plt.close(fig)

if __name__ == '__main__':
    plotBar()
    plotCorr()