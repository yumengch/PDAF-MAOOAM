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
import sys
sys.path.append('/home/yumengch/NCEO/MAOOAM/qgs-0.2.5')

from qgs.params.params import QgParams
from qgs.integrators.integrator import RungeKuttaIntegrator
from qgs.functions.tendencies import create_tendencies


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

    def __init__(self, config, pe):
        """constructor

        Parameters
        ----------
        nx : ndarray
            integer array for grid size
        nt : int
            total number of time steps
        pe : `parallelization.parallelization`
            parallelization object
        """
        # model size
        self.init_params()
        self.nx = self.model_parameters.ndim
        self.t0 = config.getint('spinup_steps', 0)
        # model size for each CPU
        self.nx_p = self.get_nx_p(pe)
        # model time steps
        self.total_steps = config.getint('total_steps', 0)

    def init_params(self):
        # Time parameters
        self.dt = 0.1
        # Setting some model parameters
        # Model parameters instantiation with default specs
        self.model_parameters = QgParams()
        # Mode truncation at the wavenumber 2 in both x and y spatial coordinate
        self.model_parameters.set_atmospheric_channel_fourier_modes(2, 2)
        # Mode truncation at the wavenumber 2 in the x and at the
        # wavenumber 4 in the y spatial coordinates for the ocean
        self.model_parameters.set_oceanic_basin_fourier_modes(2, 4)

        # Setting MAOOAM parameters according to the publication linked above
        self.model_parameters.set_params({'kd': 0.0290, 'kdp': 0.0290, 'n': 1.5, 'r': 1.e-7,
                                     'h': 136.5, 'd': 1.1e-7})
        self.model_parameters.atemperature_params.set_params({'eps': 0.7, 'T0': 289.3, 'hlambda': 15.06, })
        self.model_parameters.gotemperature_params.set_params({'gamma': 5.6e8, 'T0': 301.46})

        self.model_parameters.atemperature_params.set_insolation(103.3333, 0)
        self.model_parameters.gotemperature_params.set_insolation(310., 0)
        
        # Creating the tendencies functions
        f, _ = create_tendencies(self.model_parameters)
        # Defining an integrator
        self.integrator = RungeKuttaIntegrator(num_threads=1)
        self.integrator.set_func(f)
        self.ic = np.random.rand(self.model_parameters.ndim)*0.01

    def get_nx_p(self, pe):
        """Compute local-PE domain size/domain decomposition

        Parameters
        ----------
        pe : `parallelization.parallelization`
            parallelization object
        """

        try:
            assert pe.npes_model == 1, 'MAOOAM does not support MPI'
            return self.nx
        except AssertionError:
            print((f'...ERROR: MAOOAM does not support MPI...'))
            pe.abort_parallel()


    def init_field(self):
        """initialise PE-local model field

        Parameters
        ----------
        filename : string
            input filename
        mype_model : int
            rank of the process in model communicator
        """
        # model field
        self.field_p = self.ic.copy()


    def step(self, pe, step, usePDAF):
        """step model forward

        Parameters
        ----------
        pe : `parallelization.parallelization`
            parallelization object
        step : int
            current time step
        usePDAF : bool
            whether PDAF is used
        """
        dt = self.dt
        self.integrator.integrate((self.t0 + step - 1)*dt, (self.t0 + step)*dt, dt, ic=self.field_p, write_steps=0)
        t, self.field_p = self.integrator.get_trajectories()
        if usePDAF:
            return

        if pe.task_id == 1:
            np.savetxt(f'true_step{step}.txt', field)


    def printInfo(self, usePDAF, pe):
        """print model info

        Parameters
        ----------
        USE_PDAF : bool
            whether PDAF is used
        pe : `parallelization.parallelization`
            parallelization object
        """
        doPrint = usePDAF and pe.mype_model == 0
        doPrint = do_print or \
            (pe.task_id == 1 and pe.mype_model == 0 and not usePDAF)

        if doPrint:
            print('MODEL-side: INITIALIZE MAOOAM MODEL')
            print(f'Grid size: {self.nx}')
            print(f'Time steps {self.total_steps}')
            print(f'-- Domain decomposition over {pe.npes_model} PEs')
            print(f'-- local domain sizes (nx_p): {self.nx_p}')


    def toPhysical(self):
        # define basis functions
        Fa = lambda M, H, P, x, y: np.sqrt(2)*np.cos(P*y)
        Fk = lambda M, H, P, x, y: 2*np.cos(M*x*params_maooam.n)*np.sin(P*y)
        Fl = lambda M, H, P, x, y: 2*np.sin(H*x*params_maooam.n)*np.sin(P*y)
        phi = lambda M, H, P, x, y: 2*np.sin(0.5*H*x*params_maooam.n)*np.sin(P*y)
        basis = {'A' : Fa, 'K' : Fk, 'L' : Fl}

        # define spatial domain
        natm = params_maooam.natm
        noc = params_maooam.noc
        xx = np.linspace(0, 2*np.pi/params_maooam.n, 100, endpoint=True)
        yy = np.linspace(0, np.pi, 100, endpoint=True)
        Xc, Yc = np.meshgrid(xx, yy)

        # get atmospheric components
        psi_a = np.zeros_like(Xc)
        T_a = np.zeros_like(Xc)
        for i, ti in enumerate(awavenum):
            psi_a += X[i]*basis[ti.typ](ti.M, ti.H, ti.P, Xc, Yc)
            T_a += X[i + natm]*basis[ti.typ](ti.M, ti.H, ti.P, Xc, Yc)

        # get ocean components
        psi_o = np.zeros_like(Xc)
        T_o = np.zeros_like(Xc)
        for i, ti in enumerate(owavenum):
            psi_o += X[i + 2*natm]*basis[ti.typ](ti.M, ti.H, ti.P, Xc, Yc)
            T_o += X[i + 2*natm + noc]*basis[ti.typ](ti.M, ti.H, ti.P, Xc, Yc)

        return psi_a, T_a, psi_o, T_o
