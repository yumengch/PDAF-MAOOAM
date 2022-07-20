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


def get_obs_f(obsFactory, step, dim_obs_f, observation_f):
    """Save synthetic observations

    Parameters
    ----------
    obsFactory : 'ObsFactory.ObsFactory'
        an object containing all observations
    step : int
        current time step
    dim_obs_f : int
        dimension of observation vector
    observation_f : ndarray
        a vector of synthetic observations

    Returns
    -------
    observation_f : ndarray
        a vector of synthetic observations
    """
    istart = 0
    for obsname, obs in obsFactory.items():
        iend = istart + obs.dim_obs
        obsFactory.writer[obsname].write(step, observation_f[istart:iend])
        istart = iend

    return observation_f


def init_dim_obs_gen_pdafomi(obsFactory, local_range,
                             mype_filter, nx, nx_p, step, dim_obs):
    """initialise observation dimensions

    Parameters
    ----------
    obsFactory : dict
        a dictionary of observations
    nx_p : ndarray
        integer array for PE-local grid size
    step : int
        current time step
    dim_obs : int
        dimension of observation vector

    Returns
    -------
    dim_obs : int
        dimension of observation vector
    """
    dim_obs = 0
    for obsname, obs in obsFactory.items():
        if(obs.doassim):
            obs.init_dim_obs_gen(step, dim_obs, local_range,
                                 mype_filter, nx, nx_p)
            dim_obs += obs.dim_obs
    return dim_obs


def init_dim_obs_pdafomi(obsFactory, local_range,
                         mype_filter, nx, nx_p, step, dim_obs):
    """initialise observation dimensions

    Parameters
    ----------
    obsFactory : dict
        a dictionary of observations
    local_range : float
        range for local observation domain
    mype_filter : int
        rank of the PE in filter communicator
    nx : ndarray
        integer array for grid size
    nx_p : ndarray
        integer array for PE-local grid size
    step : int
        current time step
    dim_obs : int
        dimension of observation vector

    Returns
    -------
    TYPE
        Description
    """
    dim_obs = 0
    for obsname, obs in obsFactory.items():
        if(obs.doassim):
            obs.init_dim_obs(step, dim_obs, local_range,
                             mype_filter, nx, nx_p)
            dim_obs += obs.dim_obs
    return dim_obs


def obs_op_pdafomi(obsFactory, step, dim_p, dim_obs_p, state_p, ostate):
    """turn state vector to observation space

    Parameters
    ----------
    obsFactory : dict
        a dictionary of observations
    step : int
        current time step
    dim_p : int
        size of state vector (local part in case of parallel decomposed state)
    dim_obs_p : int
        size of observation vector
    state_p : ndarray
        local PE state vector
    ostate : ndarray
        state vector in obs space
    """
    for obsname, obs in obsFactory.items():
        ostate = obs.obs_op_gridpoint(step, state_p, ostate)
    return ostate


def init_dim_obs_l_pdafomi(obsFactory, local,
                           domain_p, step, dim_obs, dim_obs_l):
    """initialise local observation dimension

    Parameters
    ----------
    obsFactory : dict
        a dictionary of observations
    local : `Localization.Localization`
        a localization info obejct
    domain_p : int
        index of current local analysis domain
    step : int
        current time step
    dim_obs : int
        dimension of observation vector
    dim_obs_l : int
        dimension of local observation vector
    """
    for obsname, obs in obsFactory.items():
        obs.init_dim_obs_l(localization,
                           domain_p, step, dim_obs, dim_obs_l)


def localize_covar_pdafomi(obsFactory, local,
                           mype_filter, nx_p, HP_p, HPH):
    """localize covariance matrix

    Parameters
    ----------
    obsFactory : dict
        a dictionary of observations
    local : `Localization.Localization`
        a localization info obejct
    mype_filter : int
        rank of the PE in filter communicator
    nx_p : ndarray
        integer array for PE-local grid size
    HP_p : ndarray
        matrix HP
    HPH : ndarray
        matrix HPH
    """
    dim_p = HPH.shape[0]
    coords_p = np.zeros((2, dim_p))
    offset = mype_filter*nx_p

    coords_p[0] = np.arange(nx_p) + offset

    for obsname, obs in obsFactory.items():
        obs.localize_covar(localization, HP_p, HPH, coords_p)


def deallocate_obs_pdafomi(obsFactory):
    """deallocate PDAFomi object

    Parameters
    ----------
    obsFactory : dict
        a dictionary of observations
    step : int
        current time step
    """
    for obsname, obs in obsFactory.items():
        obs.deallocate_obs()
