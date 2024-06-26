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

Attributes
----------
firsttime : bool
    global variable for prepoststep_ens_pdaf
"""
import numpy as np


class PDAFUserFuncs:
    """
    Attribute
    ----------
    das : `DAS.DAS`
        an object containing DA system configurations
    firsttime : `bool`
        Whether it is the firsttime userfuncs are called
    """
    def __init__(self, das):
        self.das = das
        self.firsttime = True

    def init_ens_pdaf(self, filtertype, dim_p, dim_ens, state_p,
                      uinv, ens_p, status_pdaf):
        """initialise the ensemble for freerun where the
           initial condition is generated from a model trajectory

        Parameters
        ----------
        filtertype : int
            type of filter
        dim_p : int
            size of state vector (local part in case of parallel decomposed state)
        dim_ens : int
            size of state ensemble
        state_p : ndarray
            1D state vector on local PE
        uinv : ndarray
            2D left eigenvector with shape (dim_ens - 1,dim_ens - 1)
        ens_p : ndarray
            ensemble state vector on local PE (dim_p, dim_ens)
        status_pdaf : int
            status of PDAF

        Returns
        -------
        status_pdaf : int
            status of PDAF
        """
        # convert to physical space
        ens_p[:] = 0.
        return state_p, uinv, ens_p, status_pdaf


    def collect_state_pdaf(self, dim_p, state_p):
        """generate state vector from the model fields

        Parameters
        ----------
        dim_p : int
            size of state vector
        state_p : ndarray
            1D state vector on local PE
        """
        self.das.model.toPhysical_A()
        self.das.model.toPhysical_O()
        state_p[:] = np.concatenate([self.das.model.fields[varname].ravel(order='F') 
                                     for varname in self.das.sv.varnames])
        return state_p


    def distribute_state_pdaf(self, dim_p, state_p):
        """Distribute state vector to model field

        Parameters
        ----------
        dim_p : int
            size of the state vector
        state_p : ndarray
            1D state vector on local PE
        """
        if self.firsttime:
            if (self.das.pe.mype_world == 0):
                print ('distribute_state_pdaf: starting from restart files')
            return state_p

        size = self.das.model.nx*self.das.model.ny
        skipOcean = False
        if self.das.obs['ObsO'].doassim == 0: skipOcean = True

        varnames = self.das.sv.varnames[:].copy()
        if skipOcean: varnames = varnames[:2]
        for i, varname in enumerate(varnames):
            self.das.model.fields[varname][:] = \
                state_p[i*size:(i+1)*size].reshape(self.das.model.nx,
                                                   self.das.model.ny, order='F')

        if self.das.pe.mype_world == 0: print ('distbute atmos')
        self.das.model.toFourier_A()
        if not skipOcean:
            if self.das.pe.mype_world == 0: print ('distbute ocean')
            self.das.model.toFourier_O()
        return state_p


    def next_observation_pdaf(self, stepnow, nsteps, doexit, time):
        """The time for the next observation

        Parameters
        ----------
        stepnow : int
            Current time step
        nsteps : int
            steps between assimilation
        doexit : int
            Whether exit PDAF assimilation
        time : double
            Current model time

        Returns
        -------
        nsteps : int
            steps between assimilation
        doexit : int
            Whether exit PDAF assimilation
        time : double
            Current model time
        """
        nsteps = self.das.obs.delt_obs
        if (stepnow + nsteps <= self.das.model.total_steps):
            doexit = 0
            if (self.das.pe.mype_world == 0):
                print((stepnow, stepnow, 'Next observation at time step',
                       stepnow + nsteps))
        else:
            doexit = 1
            if (self.das.pe.mype_world == 0):
                print((stepnow, 'No more observations - end assimilation'))

        return nsteps, doexit, time


    def prepoststep_ens_pdaf(self, step, dim_p, dim_ens, dim_ens_p,
                             dim_obs_p, state_p, uinv, ens_p, flag):
        """pre- and post-processing of ensemble

        Parameters
        ----------
        step : int
            current time step (negative for call after forecast)
        dim_p : int
            pe-local state dimension
        dim_ens : int
            size of state ensemble
        dim_ens_p : int
            pe-local size of ensemble
        dim_obs_p : int
            pe-local dimension of observation vector
        state_p : ndarray[float]
            pe-local forecast/analysis state
             (the array 'state_p' is not generally not
             initialized in the case of seik.
             it can be used freely here.)
            shape is (dim_p)
        uinv : ndarray[float]
            inverse of matrix u
            shape is (dim_ens-1,dim_ens-1)
        ens_p : ndarray[float]
            pe-local state ensemble
            shape is (dim_p,dim_ens)
        flag : int
            pdaf status flag

        Returns
        -------
        state_p : ndarray[float]
            pe-local forecast/analysis state
             (the array 'state_p' is not generally not
             initialized in the case of seik.
             it can be used freely here.)
        uinv : ndarray[float]
            inverse of matrix u
        ens_p : ndarray[float]
            pe-local state ensemble
        """
        if (self.firsttime):
            print('Analyze initial state ensemble')
        else:
            if (step < 0):
                print('Analyze and write forecasted state ensemble')
            else:
                print('Analyze and write assimilated state ensemble')

        state_p = np.mean(ens_p, axis=1)
        variance_p = np.var(ens_p, axis=-1, ddof=1)

        rmserror_est = np.sqrt(np.sum(
                                variance_p
                                     )/self.das.sv.dim_state)

        if self.das.pe.mype_filter == 0:
            print('RMS error: ', rmserror_est)

        return state_p, uinv, ens_p
