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
import xarray as xr

import U_PDAFomi

import pyPDAF.PDAF as PDAF


def init_ens_pdaf(model, filtertype, dim_p,
                  dim_ens, state_p, 
                  uinv, ens_p, status_pdaf):
    """initialise the ensemble

    Parameters
    ----------
    model : `Model.Model`
        model object
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
        ensemble state vector on local PE (dim_ens, dim_p)
    status_pdaf : int
        status of PDAF

    Returns
    -------
    status_pdaf : int
        status of PDAF
    """
    # convert to physical space
    psi_a, T_a, psi_o, T_o = model.toPhysical()
    state_p[:] = np.concatenate([psi_a.ravel(), T_a.ravel(), psi_o.ravel(), T_o.ravel()])
    if dim_ens > 1:
        f = xr.load_dataset('covariance.nc')
        svals = f['sigma'].to_numpy()[:dim_ens-1]
        eofV = np.hstack([f[varname+'_svd'].to_numpy()[:dim_ens-1].reshape((dim_ens-1, dim_p//4)) for varname in ['psi_a', 'T_a', 'psi_o', 'T_o']])
        _, _, ens_p, status_pdaf = PDAF.sampleens(eofV.T, svals, state_p, verbose=1, flag=status_pdaf)
    else:
        ens_p = np.zeros((dim_p, dim_ens))
        ens_p[:, 0] = state_p.copy()
    state_p[:] = np.concatenate([psi_a.ravel(), T_a.ravel(), psi_o.ravel(), T_o.ravel()])
    return state_p, uinv, ens_p, status_pdaf


def collect_state_pdaf(model, dim_p, state_p):
    """Collect state vector in PDAF from model

    Parameters
    ----------
    model : `Model.Model`
        model object
    dim_p : int
        size of state vector
    state_p : ndarray
        1D state vector on local PE
    """
    psi_a, T_a, psi_o, T_o = model.toPhysical()
    state_p[:] = np.concatenate([psi_a.ravel(), T_a.ravel(), psi_o.ravel(), T_o.ravel()])
    return state_p


def distribute_state_pdaf(model, dim_p, state_p):
    """Distribute state vector to model field

    Parameters
    ----------
    model : `Model.Model`
        model object
    dim_p : int
        size of the state vector
    state_p : ndarray
        1D state vector on local PE
    """
    size = model.nx*model.ny
    psi_a, T_a, psi_o, T_o = (state_p[i*size:(i+1)*size].reshape(model.ny, model.nx) for i in range(4))
    model.toFourier(psi_a, T_a, psi_o, T_o)
    return state_p


def next_observation_pdaf(model, pe, delt_obs,
                          stepnow, nsteps, doexit, time):
    """The time for the next observation

    Parameters
    ----------
    model : `Model.Model`
        model object
    pe : `parallelization.parallelization`
        parallelization object
    delt_obs : int
        frequency of observations
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
    if pe.mype_world == 0:
        print (('next', stepnow))
    if (stepnow + nsteps <= model.total_steps):
        nsteps = delt_obs
        doexit = 0
        if (pe.mype_world == 0):
            print((stepnow, 'Next observation at time step',
                   stepnow + nsteps))
    else:
        nsteps = 0
        doexit = 1
        if (pe.mype_world == 0):
            print((stepnow, 'No more observations - end assimilation'))

    return nsteps, doexit, time


firsttime = True


def prepoststep_ens_pdaf(sv, model, pe, obs, writer,
                         step, dim_p, dim_ens, dim_ens_p,
                         dim_obs_p, state_p, uinv, ens_p, flag):
    """pre- and post-processing of ensemble

    Parameters
    ----------
    model : `Model.Model`
        model object
    pe : `parallelization.parallelization`
        MPI information
    obs: `ObsFactory.ObsFactory`
        Observations information
    writer: `ModelWriter.ModelWriter`
        Model writer
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
    if pe.mype_world == 0:
        print (('prepost', step))
    global firsttime
    if (firsttime):
        print('Analyze initial state ensemble')
    else:
        if (step < 0):
            print('Analyze and write forecasted state ensemble')
        else:
            print('Analyze and write assimilated state ensemble')

    variance = np.zeros(sv.dim_state)
    state_p = np.mean(ens_p, axis=1)
    variance_p = np.var(ens_p, axis=-1, ddof=1)

    if pe.mype_filter != 0:
        pe.COMM_filter.Send(variance_p, 0, pe.mype_filter)
    else:
        variance[:dim_p] = variance_p
        for i in range(1, pe.npes_filter):
            pe.COMM_filter.Recv(
                              variance[i*dim_p:(i+1)*dim_p], i, i)

    rmserror_est = np.sqrt(np.sum(
                            variance
                                 )/sv.dim_state)

    if pe.mype_filter == 0:
        print('RMS error: ', rmserror_est)

    if step < 0:
        print('---writing ensemble forecast---')
        writer.write(-step, 'f', ens_p)
    else:
        print('---writing ensemble analysis---')
        writer.write(step, 'a', ens_p)

    # U_PDAFomi.deallocate_obs_pdafomi(obs)

    firsttime = False

    return state_p, uinv, ens_p
