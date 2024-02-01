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
            obs.init_dim_obs_gen(step, self.das.model)
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
        for obsname, obs in self.das.obs.items():
            if obs.doassim == 0: continue
            obs.init_dim_obs(step, self.das.model)
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
            if obs.doassim == 1:
                ostate = obs.obs_op_gridpoint(step, state_p, ostate)
        return ostate
