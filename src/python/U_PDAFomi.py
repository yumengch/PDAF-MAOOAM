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


class PDAFomiUserFuncs:
    """
    Attribute
    ----------
    das : 'DAS.DAS'
        an object containing DA system configurations
    """
    def __init__(self, das):
        self.das = das


    def get_obs_f(self, step, dim_obs_f, observation_f):
        """Save synthetic observations

        Parameters
        ----------
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
        for obsname, obs in self.das.obs.items():
            iend = istart + obs.dim_obs
            self.das.obs.writer[obsname].write(step, observation_f[istart:iend])
            istart = iend

        return observation_f


    def init_dim_obs_gen_pdafomi(self, step, dim_obs):
        """initialise observation dimensions

        Parameters
        ----------
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
        for obsname, obs in self.das.obs.items():
            if(obs.doassim):
                obs.init_dim_obs_gen(step, dim_obs, self.das.model,
                                     self.das.pe.mype_filter,
                                     self.das.sv.dim_state,
                                     self.das.sv.dim_state_p)
                dim_obs += obs.dim_obs
        return dim_obs


    def init_dim_obs_pdafomi(self, step, dim_obs):
        """initialise observation dimensions

        Parameters
        ----------
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
        condition = (not self.das.isStrong) and (self.das.sv.component == 'ao')
        factor = 2 if condition else 1
        shift = -1 if condition else 0
        true_step = (step - shift)//factor
        for obsname, obs in self.das.obs.items():
            obs.init_dim_obs(true_step, dim_obs, self.das.model, 
                             self.das.isStrong,
                             self.das.pe.mype_filter,
                             self.das.sv.dim_state, 
                             self.das.sv.dim_state_p)
            dim_obs += obs.dim_obs
        return dim_obs


    def obs_op_pdafomi(self, step, dim_p, dim_obs_p, state_p, ostate):
        """turn state vector to observation space

        Parameters
        ----------
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
        for obsname, obs in self.das.obs.items():
            if(obs.doassim):
                ostate = obs.obs_op_gridpoint(step, state_p, ostate)
        return ostate


    def init_dim_obs_l_pdafomi(self, domain_p, step, dim_obs, dim_obs_l):
        """initialise local observation dimension

        Parameters
        ----------
        domain_p : int
            index of current local analysis domain
        step : int
            current time step
        dim_obs : int
            dimension of observation vector
        dim_obs_l : int
            dimension of local observation vector
        """
        dim_obs_l = 0
        for obsname, obs in self.das.obs.items():
            if(obs.doassim):
                dim_obs_l += obs.init_dim_obs_l(
                                   domain_p, step, dim_obs, dim_obs_l)
        return dim_obs_l

    def localize_covar_pdafomi(self, dim_state_p, HP_p, HPH):
        """localize covariance matrix

        Parameters
        ----------
        dim_state_p : float
            integer array for PE-local grid size
        HP_p : ndarray
            matrix HP
        HPH : ndarray
            matrix HPH
        """
        dim_p = HPH.shape[0]
        coords_p = np.zeros((2, dim_p))
        coords_p[0] = np.tile(self.das.model.xc.ravel(), self.das.sv.nVar)
        coords_p[1] = np.tile(self.das.model.yc.ravel(), self.das.sv.nVar)

        for obsname, obs in self.das.obs.items():
            if(obs.doassim):
                HP_p, HPH = obs.localize_covar(HP_p, HPH, coords_p)
        return HP_p, HPH

    def deallocate_obs_pdafomi(self):
        """deallocate PDAFomi object

        Parameters
        ----------
        step : int
            current time step
        """
        for obsname, obs in self.das.obs.items():
            obs.deallocate_obs()
