"""This file is part of pyPDAF

Copyright (C) 2022 University of Reading and
National Centre for Earth Observation

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import numpy as np
import xarray as xr
import scipy
import sys
sys.path.append('/home/users/yumengch/qgs-0.2.5/')

from qgs.params.params import QgParams
from qgs.integrators.integrator import RungeKuttaIntegrator
from qgs.functions.tendencies import create_tendencies

from ModelWriter import ModelWriter


class Model:

    """Model information in PDAF

    Attributes
    ----------
    field_p : ndarray
        PE-local model field
    nx : ndarray
        integer array for grid size
    nx_p : ndarray
        integer array for PE-local grid size
    total_steps : int
        total number of time steps
    """

    def __init__(self, config, task_id):
        """constructor

        Parameters
        ----------
        config : Config.PDAFConfig
        task_id : int
            the task_id-th ensemble member
        """
        # init model size
        self.init_params(config.getfloat('n', 1.5))
        # init physical grid
        self.nx = config.getint('nx', 9)
        self.ny = config.getint('ny', 9)
        x0 = 0; x1 = 2*np.pi/self.model_parameters.scale_params.n
        y0 = 0; y1 = np.pi
        X = np.linspace(x0, x1, self.nx)
        Y = np.linspace(y0, y1, self.ny)
        self.xc, self.yc = np.meshgrid(X, Y)
        # init model time steps
        self.total_steps = config.getint('total_steps', 0)
        self.tw = config.getint('tw', 0)
        self.dt = 0.1
        # init fields in physical space
        self.varnames = ['psi_a', 'T_a', 'psi_o', 'T_o']
        self.fields = {varname:np.zeros_like(self.xc) for varname in self.varnames}
        # model initial state
        ln_restart = config.getboolean('ln_restart', False)
        restart_it = config.getint('restart_it', 0)
        self.init_field(ln_restart, restart_it, task_id)
        # init writer
        self.writer = ModelWriter('maooam_{:03}.nc'.format(task_id), 
                                  *self.model_parameters.nmod)
        # write the init time step
        self.writer.write(self.t0*self.dt, 'f', self.field_p)
        # # model index
        # self.task_id = task_id


    def init_params(self, n):
        # Setting some model parameters
        # Model parameters instantiation with default specs
        self.model_parameters = QgParams({'n' : n, 'kd': 0.02, 'kdp': 0.02, 'n': n, 'r': 1.e-7,
                                     'h': 165, 'd': 9e-8})
        # Mode truncation at the wavenumber 2 in both x and y spatial coordinate
        self.model_parameters.set_atmospheric_channel_fourier_modes(2, 2)
        # Mode truncation at the wavenumber 2 in the x and at the
        # wavenumber 4 in the y spatial coordinates for the ocean
        self.model_parameters.set_oceanic_basin_fourier_modes(2, 4)

        # Setting MAOOAM parameters according to the publication linked above
        self.model_parameters.set_params({'kd': 0.02, 'kdp': 0.02, 'n': n, 'r': 1.e-7,
                                          'h': 165, 'd': 9e-8})
        self.model_parameters.atemperature_params.set_params({'eps': 0.7, 'T0': 290.2, 'hlambda': 15.06, })
        self.model_parameters.gotemperature_params.set_params({'gamma': 6.6e8, 'T0': 299.35})
        self.model_parameters.atemperature_params.set_insolation(103.3333, 0)
        self.model_parameters.gotemperature_params.set_insolation(310, 0)
        
        # Creating the tendencies functions
        f, _ = create_tendencies(self.model_parameters)
        # Defining an integrator
        self.integrator = RungeKuttaIntegrator(num_threads=1)
        self.integrator.set_func(f)


    def init_field(self, ln_restart, restart_it, task_id):
        if ln_restart:
            f = xr.open_dataset('restart/maooam_{:03}.nc'.format(task_id), decode_times=False)
            self.t0 = f['time'][restart_it].to_numpy()
            self.field_p = np.concatenate([f[varname+'_f'][restart_it].to_numpy() for varname in self.varnames])
            f.close()
        else:
            self.field_p = np.zeros(self.model_parameters.ndim)
            self.t0 = 0


    def step(self, step):
        """step model forward

        Parameters
        ----------
        step : int
            current time step
        """
        dt = self.dt
        self.integrator.integrate((self.t0 + step)*dt, (self.t0 + step + 1)*dt, dt, ic=self.field_p, write_steps=0)
        t, self.field_p = self.integrator.get_trajectories()
        return t


    def printInfo(self, pe):
        """print model info

        Parameters
        ----------
        USE_PDAF : bool
            whether PDAF is used
        pe : `parallelization.parallelization`
            parallelization object
        """
        doPrint = pe.mype_model == 0
        doPrint = do_print or \
            (pe.task_id == 1 and pe.mype_model == 0)

        if doPrint:
            print('MODEL-side: INITIALIZE MAOOAM MODEL')
            print(f'Grid size: {self.nx} x {self.ny}')
            print(f'Time steps {self.total_steps}')
            print(f'-- Domain decomposition over {pe.npes_model} PEs')
            print(f'-- local domain sizes (nx_p): {self.nx*self.ny}')


    def toPhysical_A(self):
        # define spatial domain
        natm, _ = self.model_parameters.nmod
        self.fields['psi_a'][:] = 0.
        self.fields['T_a'][:] = 0.
        # get atmospheric components
        basis = self.model_parameters.atmospheric_basis.num_functions()
        for i, b in enumerate(basis):
            self.fields['psi_a'][:] += self.field_p[i]*b(self.xc, self.yc)
            self.fields['T_a'][:] += self.field_p[i+natm]*b(self.xc, self.yc)


    def toPhysical_O(self):
        # define spatial domain
        natm, noc = self.model_parameters.nmod
        self.fields['psi_o'][:] = 0.
        self.fields['T_o'][:] = 0.
        # get ocean components
        basis = self.model_parameters.oceanic_basis.num_functions()
        for i, b in enumerate(basis):
            self.fields['psi_o'][:] += self.field_p[i+2*natm]*b(self.xc, self.yc)
            self.fields['T_o'][:] += self.field_p[i+2*natm+noc]*b(self.xc, self.yc)


    def toFourier_A(self):
        # define spatial domain
        natm, _ = self.model_parameters.nmod
        dx = self.xc[0, 1] - self.xc[0, 0]
        dy = self.yc[1, 0] - self.yc[0, 0]

        # get atmospheric components
        basis = self.model_parameters.atmospheric_basis.num_functions()
        for i, b in enumerate(basis):
            self.field_p[i] = scipy.integrate.romb(
                            scipy.integrate.romb(
                                self.fields['psi_a']*b(self.xc, self.yc), dx=dx
                                ), dx=dy
                            )
            self.field_p[i+natm] = scipy.integrate.romb(
                                scipy.integrate.romb(
                                    self.fields['T_a']*b(self.xc, self.yc), dx=dx
                                    ), dx=dy
                                    )
        self.field_p[:2*natm] =  self.field_p[:2*natm]*self.model_parameters.scale_params.n/2/np.pi/np.pi


    def toFourier_O(self):
        # define spatial domain
        natm, noc = self.model_parameters.nmod
        dx = self.xc[0, 1] - self.xc[0, 0]
        dy = self.yc[1, 0] - self.yc[0, 0]

        # get ocean components
        basis = self.model_parameters.oceanic_basis.num_functions()
        for i, b in enumerate(basis):
            self.field_p[i + 2*natm] = scipy.integrate.romb(
                            scipy.integrate.romb(
                                self.fields['psi_o']*b(self.xc, self.yc), dx=dx
                                ), dx=dy
                            )
            self.field_p[i + 2*natm + noc] = scipy.integrate.romb(
                                scipy.integrate.romb(
                                    self.fields['T_o']*b(self.xc, self.yc), dx=dx
                                    ), dx=dy
                                )
        self.field_p[2*natm:] =  self.field_p[2*natm:]*self.model_parameters.scale_params.n/2/np.pi/np.pi
