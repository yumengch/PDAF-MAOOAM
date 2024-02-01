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

import mod_model
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
        mod_model.mod_model.initialize_model()
        # init physical grid
        self.nx = config.getint('nx', 9)
        self.ny = config.getint('ny', 9)
        x0 = 0; x1 = 2*np.pi/config.getfloat('n', 1.5)
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
        self.fields = {varname:np.zeros((self.nx, self.ny), order='F') for varname in self.varnames}
        # model initial state
        restart_it = config.getint('restart_it', 0)
        self.init_field(restart_it, task_id)
        self.t = self.t0
        self.field_p_new = self.field_p.copy()
        # init writer
        self.writer = ModelWriter('maooam_{:03}.nc'.format(task_id),
                                  mod_model.mod_model.natm, mod_model.mod_model.noc)
        # write the init time step
        self.writer.write(self.t0*self.dt, 'f', self.field_p)


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


    def init_field(self, restart_it, task_id):
        f = xr.open_dataset('restart/maooam_{:03}.nc'.format(task_id), decode_times=False)
        self.t0 = f['time'][restart_it].to_numpy()
        self.field_p = np.concatenate([f[varname+'_a'][restart_it].to_numpy() for varname in self.varnames])
        f.close()


    def step(self):
        """step model forward
        """
        mod_model.mod_model.step(self.field_p, self.t, self.field_p_new)
        self.field_p[:] = self.field_p_new[:]
        return self.t


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
        self.fields['psi_a'][:] = 0.
        self.fields['T_a'][:] = 0.
        mod_model.mod_model.tophysical_a(self.field_p, self.fields['psi_a'], self.fields['T_a'])


    def toPhysical_O(self):
        self.fields['psi_o'][:] = 0.
        self.fields['T_o'][:] = 0.
        mod_model.mod_model.tophysical_o(self.field_p, self.fields['psi_o'], self.fields['T_o'])


    def toFourier_A(self):
        mod_model.mod_model.tofourier_a(self.nx, self.ny, self.fields['psi_a'], self.fields['T_a'], self.field_p)


    def toFourier_O(self):
        mod_model.mod_model.tofourier_o(self.nx, self.ny, self.fields['psi_o'], self.fields['T_o'], self.field_p)

    def __del__(self):
        mod_model.mod_model.finalize_model()