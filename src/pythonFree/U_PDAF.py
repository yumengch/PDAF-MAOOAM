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
import netCDF4

import pyPDAF.PDAF as PDAF

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
        self.das.model.toPhysical_A()
        self.das.model.toPhysical_O()
        state_p[:] = np.concatenate([self.das.model.fields[varname].ravel(order='F') for varname in self.das.sv.varnames])
        if dim_ens > 1:
            f = netCDF4.Dataset('covariance.nc')
            svals = f['sigma'][:dim_ens-1].data
            eofV = np.hstack([f[varname+'_svd'][:dim_ens-1].data.reshape((dim_ens-1,
                                                                          dim_p//4))
                              for varname in ['psi_a', 'T_a', 'psi_o', 'T_o']])
            f.close()
            _, _, ens_p, status_pdaf = PDAF.sampleens(eofV.T, svals, state_p,
                                                      verbose=1, flag=status_pdaf)

            ens_mean = np.mean(ens_p, axis=1, keepdims=True)
            offset = self.das.model.nx*self.das.model.ny
            for i, scale in enumerate(self.das.model.ensscale):
                ens_p[i*offset:(i+1)*offset] = scale*(ens_p[i*offset:(i+1)*offset] - \
                                                ens_mean[i*offset:(i+1)*offset]) + \
                                                    ens_mean[i*offset:(i+1)*offset]
        else:
            ens_p = np.zeros((dim_p, dim_ens))
            ens_p[:, 0] = state_p.copy()

        state_p[:] = np.concatenate([self.das.model.fields[varname].ravel(order='F') for varname in self.das.sv.varnames])

        return state_p, uinv, ens_p, status_pdaf


    def distribute_state_pdaf(self, dim_p, state_p):
        """Distribute state vector to model field

        Parameters
        ----------
        dim_p : int
            size of the state vector
        state_p : ndarray
            1D state vector on local PE
        """
        size = self.das.model.nx*self.das.model.ny
        for i, varname in enumerate(self.das.sv.varnames):
            self.das.model.fields[varname][:] = \
                state_p[i*size:(i+1)*size].reshape(self.das.model.nx,
                                                   self.das.model.ny, order='F')

        if self.das.pe.mype_world == 0: print ('distbute atmos')
        self.das.model.toFourier_A()
        if self.das.pe.mype_world == 0: print ('distbute ocean')
        self.das.model.toFourier_O()
        return state_p


    def next_observation_pdaf(self,
                              stepnow, nsteps, doexit, time):
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
        nsteps = 1
        doexit = 0

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
        print('Analyze initial state ensemble')

        variance = np.zeros(self.das.sv.dim_state)
        state_p = np.mean(ens_p, axis=1)
        variance_p = np.var(ens_p, axis=-1, ddof=1)

        if self.das.pe.mype_filter != 0:
            self.das.pe.COMM_filter.Send(variance_p, 0, self.das.pe.mype_filter)
        else:
            variance[:dim_p] = variance_p
            for i in range(1, self.das.pe.npes_filter):
                self.das.pe.COMM_filter.Recv(
                                  variance[i*dim_p:(i+1)*dim_p], i, i)

        rmserror_est = np.sqrt(np.sum(
                                variance
                                     )/self.das.sv.dim_state)

        if self.das.pe.mype_filter == 0:
            print('RMS error: ', rmserror_est)

        return state_p, uinv, ens_p
