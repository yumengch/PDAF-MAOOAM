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
import configparser

import numpy as np
import xarray as xr
import mod_model # type: ignore

import parallelisation
from model_writer import model_writer
import log

class model:

    """Model information in PDAF

    Attributes
    ----------
    field_p : ndarray
        PE-local model field in spectral space
    nx : int
        number of grid points in x direction
    ny : int
        number of grid points in y direction
    xc : ndarray
        x-coordinate of the grid points
    yc : ndarray
        y-coordinate of the grid points
    dt : float
        time step size
    t : float
        current time
    t0 : float
        initial time
    tw : int
        time step for writing output
    varnames : list
        variable names
    fields : dict
        model fields in physical space
    writer : model_writer
        writer object
    total_steps : int
        total number of time steps

    Methods
    -------
    init_params(n:int)
        initialise model parameters
    init_field(restart_it:int, task_id:int)
        initialise model field from restart file
    step() -> float
        step model forward
    printInfo(pe:parallelisation.parallelisation)
        print model info
    toPhysical_A()
        convert atmosphere fields to physical space
    toPhysical_O()
        convert ocean fields to physical space
    toFourier_A()
        convert atmosphere fields to spectral space
    toFourier_O()
        convert ocean fields to spectral space
    """

    def __init__(self, config:configparser.SectionProxy, task_id:int) -> None:
        """constructor

        Parameters
        ----------
        config : configparser.SectionProxy
        task_id : int
            the task_id-th ensemble member
        """
        # init model size
        mod_model.mod_model.initialize_model()
        # init physical grid
        self.nx:int = config.getint('nx', 9)
        self.ny:int = config.getint('ny', 9)
        x0:float = 0; x1:float = 2*np.pi/config.getfloat('n', 1.5)
        y0:float = 0; y1:float = np.pi
        X:np.ndarray = np.linspace(x0, x1, self.nx)
        Y:np.ndarray = np.linspace(y0, y1, self.ny)
        self.xc:np.ndarray
        self.yc:np.ndarray
        self.xc, self.yc = np.meshgrid(X, Y)
        # init model time steps
        self.total_steps:int = config.getint('total_steps', 0)
        self.tw:int = config.getint('tw', 0)
        self.dt:float = 0.1
        # init fields in physical space
        self.varnames:list[str] = ['psi_a', 'T_a', 'psi_o', 'T_o']
        self.fields:dict[str, np.ndarray] = {varname:np.zeros((self.nx, self.ny), order='F') for varname in self.varnames}
        # model initial state
        restart_it:int = config.getint('restart_it', 0)
        self.init_field(restart_it, task_id)
        self.t:float = self.t0
        self.field_p_new:np.ndarray = self.field_p.copy()
        # init writer
        self.writer:model_writer = model_writer('maooam_{:03}.nc'.format(task_id),
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


    def init_field(self, restart_it:int, task_id:int) -> None:
        """initialize model field from  restart file"""
        f:xr.Dataset = xr.open_dataset('restart/maooam_{:03}.nc'.format(task_id), decode_times=False)
        self.t0:float = float(f['time'][restart_it].to_numpy())
        self.field_p:np.ndarray = np.concatenate([f[varname+'_a'][restart_it].to_numpy() for varname in self.varnames])
        f.close()


    def step(self) -> float:
        """step model forward
        """
        mod_model.mod_model.step(self.field_p, self.t, self.field_p_new)
        self.field_p[:] = self.field_p_new[:]
        return self.t


    def printInfo(self, pe:parallelisation.parallelisation) -> None:
        """print model info

        Parameters
        ----------
        USE_PDAF : bool
            whether PDAF is used
        pe : `parallelization.parallelization`
            parallelization object
        """
        doPrint:bool = pe.mype_model == 0
        doPrint = pe.task_id == 1 and pe.mype_model == 0

        if doPrint:
            log.logger.info('MODEL-side: INITIALIZE MAOOAM MODEL')
            log.logger.info(f'Grid size: {self.nx} x {self.ny}')
            log.logger.info(f'Time steps {self.total_steps}')
            log.logger.info(f'-- Domain decomposition over {pe.npes_model} PEs')
            log.logger.info(f'-- local domain sizes (nx_p): {self.nx*self.ny}')


    def toPhysical_A(self):
        """Convert atmosphere fields to physical space
        """
        # define spatial domain
        self.fields['psi_a'][:] = 0.
        self.fields['T_a'][:] = 0.
        mod_model.mod_model.tophysical_a(self.field_p, self.fields['psi_a'], self.fields['T_a'])


    def toPhysical_O(self):
        """Convert ocean fields to physical space
        """
        self.fields['psi_o'][:] = 0.
        self.fields['T_o'][:] = 0.
        mod_model.mod_model.tophysical_o(self.field_p, self.fields['psi_o'], self.fields['T_o'])


    def toFourier_A(self):
        """Convert atmosphere fields to spectral space
        """
        mod_model.mod_model.tofourier_a(self.nx, self.ny, self.fields['psi_a'], self.fields['T_a'], self.field_p)


    def toFourier_O(self):
        """Convert ocean fields to spectral space
        """
        mod_model.mod_model.tofourier_o(self.nx, self.ny, self.fields['psi_o'], self.fields['T_o'], self.field_p)

    def finalise(self):
        mod_model.mod_model.finalize_model()
        self.writer.finalise()
