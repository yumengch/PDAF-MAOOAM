import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgs
import matplotlib.ticker as mticker
import cmocean


_varnames = ['psi_a', 'T_a', 'psi_o', 'T_o']
_titles = [r'$\Psi_a$', r'$T_a$', r'$\Psi_o$', r'$T_o$']
_percentage = [0.5, 0.5, 0.7, 0.7]
_nt = 100000
_nx = 129
_ny = 129
_Ne = 16

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['font.size'] = 20


def getRMSE(error):
    return np.sqrt((error*error).mean(axis=(1, 2)))

def plotError():
    f = xr.open_dataset('traj_var.nc', decode_times=False)
    fig = plt.figure(1)
    fig.clf()
    w, h = fig.get_size_inches()
    fig.set_size_inches(2*w, 2*h)
    gs = mgs.GridSpec(2, 2, left=0.025, bottom=0.06, right=0.948, top=0.96, wspace=0.13, hspace=0.2)
    x = np.linspace(0, np.pi*2/1.5, 17)
    y = np.linspace(0, np.pi, 17)
    locator_x = mticker.FixedLocator([0, np.pi/1.5, np.pi*2/1.5])
    formatter_x = mticker.FixedFormatter(['0', r'$\frac{\pi}{n}$', r'$\frac{2\pi}{n}$'])
    locator_y = mticker.FixedLocator([0, np.pi*.5, np.pi])
    formatter_y = mticker.FixedFormatter(['0', r'$\frac{\pi}{2}$', r'$\pi$'])
    for i, (varname, title, p) in enumerate(zip(_varnames, _titles, _percentage)):
        ax = fig.add_subplot(gs[i])
        mappable = ax.pcolormesh(x, y, p*np.sqrt(f[varname+'_var']), cmap='cmo.amp')
        ax.xaxis.set_major_locator(locator_x)
        ax.xaxis.set_major_formatter(formatter_x)
        ax.yaxis.set_major_locator(locator_y)
        ax.yaxis.set_major_formatter(formatter_y)
        ax.set_title(f'{title} ({np.round(np.mean(p*np.sqrt(f[varname+"_var"])).to_numpy(), 4)})')
        plt.colorbar(mappable, ax=ax)

    fig.savefig(f'R.eps', dpi=300)
    plt.close(fig)

if __name__ == '__main__':
    plotError()