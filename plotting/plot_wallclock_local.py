"""plot the wallclock time of the different Python and Fortran implementations under different resolutions"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as mgs


_funcnames = ['PDAF-internal operations','prepoststep_pdaf', 'distribute_state_pdaf',
              'collect_state_pdaf',
              'MPI communication in PDAF', 'obs_op_pdafomi', 'OMI-internal routines',
              'init_dim_obs_pdafomi', 'init_n_domains_pdaf', 'init_dim_l_pdaf',
              'g2l_state_pdaf', 'l2g_state_pdaf']
_labels = ['internal','pre-post', 'distribute state',
              'collect state',
              'MPI', 'obs. operator', 'OMI-internal',
              'OMI setup', 'no. domains', 'init local domain',
              'g2l state',
              'l2g state', 'total']
_overallnames = ['Initialize PDAF', 'Ensemble forecast', 'LETKF analysis:', 'prepoststep_pdaf']

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

    with open(filename) as f:
        for line in f:
            for funcname in _funcnames:
                if funcname in line and 'PDAF' in line:
                    time[funcname] += float(line.split()[-2])/1000.0
            for funcname in _overallnames:
                if funcname in line and 'PDAF' in line:
                    time['overall'] += float(line.split()[-2])/1000.0

    return time

def plot_bar_local() -> None:
    """Plot the wallclock time of the different functions for the different resolutions"""
    time_py : list[dict[str, float]] = []
    time_fortran : list[dict[str, float]] = []

    # retrieve wall clock time for python and fortran system
    legends = ['every 8 gp', 'every 4 gp', 'every 8 gp (PDAFlocal)', 'every 4 gp (PDAFlocal)']
    _filepath_fortran: list[str]=['/home/users/ia923171/maooam_exps/strong_local_257/myout.txt',
                                  '/home/users/ia923171/maooam_exps/strong_local_257_4/myout.txt',
                                  '/home/users/ia923171/maooam_exps/strong_local_257/myout.txt',
                                  '/home/users/ia923171/maooam_exps/strong_local_257_4/myout.txt',]
    _filepath_py: list[str]=['/home/users/ia923171/maooam_exps/strong_local_py_257/myout.txt',
                             '/home/users/ia923171/maooam_exps/strong_local_py_257_4/myout.txt',
                             '/home/users/ia923171/maooam_local_exps/no_callback_257/myout.txt',
                             '/home/users/ia923171/maooam_local_exps/no_callback_257_4/myout.txt',]
    for filename in _filepath_py:
        time_py.append(read_time(filename))

    for filename in _filepath_fortran:
        time_fortran.append(read_time(filename))

    fig = plt.figure()
    gs = mgs.GridSpec(1, 1, figure=fig, left=0.08, right=1., top=0.99, bottom=0.33, wspace=0.3)
    w, h = fig.get_size_inches()
    fig.set_size_inches(w*2, h)

    colors = ['#648fff', '#785ef0', '#dc267f',
              '#fe6100', ]
    width = 0.2
    ind = np.arange(len(_funcnames) + 1)*2.5 - 5*width

    ax = fig.add_subplot(gs[0])
    for i, (t_py, t_fort, c, legend) in enumerate(zip(time_py, time_fortran, colors, legends)):
        ax.bar((ind - 5*width + 2*i*width)[5:], np.array(list(t_fort.values()))[5:],
               width, color=c, edgecolor='k', alpha=0.3, label=f'{legend} (fort)')
        ax.bar((ind - 5*width + (2*i+1)*width)[5:], np.array(list(t_py.values()))[5:],
               width, color=c, label=f'{legend} (py)')
        print (legend, t_fort['overall'], t_py['overall'])

        ax.set_xticks(ind[5:])
        ax.set_xticklabels([label for label in _labels[5:]], rotation=25)
        ax.set_yscale('log')
        ax.set_ylabel('Time per analysis step (s)')
        ax.legend(loc=(0.08, -0.5), ncols=4, fontsize=11)

    fig.savefig('wallclock_time_local.pdf', dpi=300)


if __name__ == '__main__':
    plot_bar_local()