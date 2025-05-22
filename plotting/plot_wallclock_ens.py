"""plot the wallclock time of the different
Python and Fortran implementations under different resolutions
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgs


_funcnames = ['PDAF-internal operations','prepoststep_pdaf', 'distribute_state_pdaf',
              'collect_state_pdaf',
              'MPI communication in PDAF', 'obs_op_pdafomi', 'OMI-internal routines',
              'init_dim_obs_pdafomi']
_labels = ['internal','pre-post', 'distribute state',
              'collect state',
              'MPI', 'obs. operator', 'OMI-internal',
              'OMI setup', 'total']
_overallnames = ['Initialize PDAF', 'Ensemble forecast', 'ETKF analysis', 'prepoststep_pdaf']
ngrids : np.ndarray = 2**np.arange(7, 12) + 1
_filename_fortran : str = '/home/users/ia923171/maooam_exps_add_revision/strong_{mem}mem/myout.txt'
_filename_py : str = '/home/users/ia923171/maooam_exps_add_revision/strong_py_{mem}mem/myout.txt'


plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['font.size'] = 16

def read_time(filename: str) -> dict[str, float]:
    """Read the wallclock time of the different functions from the output file

    Parameters
    ----------
    filename : str
        The name of the output file
    """
    time : dict[str, float] = dict()
    for funcname in _funcnames:
        time[funcname] = 0.0
    time['overall'] = 0.0

    with open(filename, 'r') as f:
        for line in f:
            for funcname in _funcnames:
                if funcname in line and 'PDAF' in line:
                    time[funcname] += float(line.split()[-2])/1000.0
            for funcname in _overallnames:
                if funcname in line and 'PDAF' in line:
                    time['overall'] += float(line.split()[-2])/1000.0

    return time


def plot_bar() -> None:
    """Plot the wallclock time of the different functions for the different resolutions"""
    time_py : list[dict[str, float]] = []
    time_fortran : list[dict[str, float]] = []

    # retrieve wall clock time for python and fortran system
    for mem in ['16', '64', '128']:
        time_py.append(read_time(_filename_py.format(mem=mem)))
        time_fortran.append(read_time(_filename_fortran.format(mem=mem)))
    # print (time_py)
    print (time_fortran)
    for name in time_py[1]:
        print (name, time_py[1][name], time_fortran[1][name], time_py[1][name]/time_fortran[1][name])
    for name in time_py[2]:
        print (name, time_py[2][name], time_fortran[2][name], time_py[2][name]/time_fortran[2][name])
    fig = plt.figure()
    gs = mgs.GridSpec(1, 1, figure=fig, left=0.095, right=1., top=0.99, bottom=0.28, wspace=0.3)
    w, h = fig.get_size_inches()
    fig.set_size_inches(w*2, h)

    colors = ['#377eb8', '#ff7f00', '#4daf4a',
              '#f781bf', '#a65628', '#984ea3',
              '#999999', '#e41a1c', '#dede00']
    width = 0.2
    ind = np.arange(len(_funcnames) + 1)*2.5 - 5*width

    ax = fig.add_subplot(gs[0])
    for i, (mem, c) in enumerate(zip(['16', '64', '128'], colors)):
        ax.bar(ind - 5*width + 2*i*width, np.array(list(time_fortran[i].values())),
               width, color=c, edgecolor='k', alpha=0.3,
               label=f'{mem} members (fort)')
        ax.bar(ind - 5*width + (2*i+1)*width, np.array(list(time_py[i].values())),
               width, color=c, label=f'{mem} members (py)')
        print (time_fortran[i]['init_dim_obs_pdafomi'], time_py[i]['init_dim_obs_pdafomi'])
        ax.set_xticks(ind)
        ax.set_xticklabels([label for label in _labels], rotation=15)
        ax.set_yscale('log')
        ax.set_ylabel('Time per analysis step (s) \n in log-scale')
        ax.legend(loc=(-0.05, -0.38), ncols=5, fontsize=11)

    fig.savefig('wallclock_time_ens.pdf', dpi=300)


if __name__ == '__main__':
    plot_bar()
