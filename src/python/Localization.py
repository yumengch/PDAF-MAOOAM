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
import pyPDAF.PDAF.PDAFomi as PDAFomi


class Localization:

    """class for localization information and user-supplied functions

    Attributes
    ----------
    coords_l : ndarray
        coordinates of local analysis domain
    id_lstate_in_pstate : ndarray
        indices of local state vector in PE-local global state vector
    loc_weight : int
        - (0) constant weight of 1
        - (1) exponentially decreasing with SRANGE
        - (2) use 5th-order polynomial
        - (3) regulated localization of R with mean error variance
        - (4) regulated localization of R with single-point error variance
    local_range : float
        range for local observation domain
    srange : float
        support range for 5th order polynomial
        or radius for 1/e for exponential weighting
    """

    def __init__(self, config):
        """class constructor

        Parameters
        ----------

        """
        # - (0) constant weight of 1
        # - (1) exponentially decreasing with SRANGE
        # - (2) use 5th-order polynomial
        # - (3) regulated localization of R with mean error variance
        # - (4) regulated localization of R with single-point error variance
        self.loc_weight = config.getint('loc_weight', 0)
        # range for local observation domain
        self.local_range = config.getint('local_range', 0)
        # support range for 5th order polynomial
        # or radius for 1/e for exponential weighting 
        self.srange = config.getint('srange', 0)

    def init_dim_l_pdaf(self, model, mype_filter, step, domain_p, dim_l):
        """initialise the local dimension of PDAF.

        The function set coordinates of local analysis domain
        and its conversion to global state vector. It also returns
        the dimension of local state vector

        Parameters
        ----------
        model : Model.Model
            model object
        mype_filter : int
            rank of the process in filter communicator
        step : int
            current time step
        domain_p : int
            index of current local analysis domain
        dim_l : int
            dimension of local state vector

        Returns
        -------
        dim_l : int
            dimension of local state vector
        """
        # initialize local state dimension
        dim_l = 1
        # initialize coordinates of local domain
        # we use grid point indices as coordinates,
        #  but could e.g. use meters
        self.coords_l = np.zeros(2)
        self.coords_l[0] = model.xc.ravel()[domain_p]
        self.coords_l[1] = model.yc.ravel()[domain_p]
        # initialize array of indices of the local state
        #  vector elements in the global state vector.

        # allocate array
        self.id_lstate_in_pstate = np.zeros(dim_l)

        # here the local domain is a single grid point
        # and variable given by domain_p
        self.id_lstate_in_pstate[0] = domain_p

        return dim_l

    def init_n_domains_pdaf(self, assim_dim, step):
        """get the number of analysis domains

        Parameters
        ----------
        assim_dim : `AssimilationDimension.AssimilationDimension`
            assim_dim object
        step : int
            current time step

        Returns
        -------
        n_domains_p : int
            PE-local number of analysis domains
        """
        return assim_dim.dim_state_p

    def g2l_state_pdaf(self, step, domain_p, state_p, state_l):
        """convert PE-local global state vector to local state vector

        Parameters
        ----------
        step : int
            current time step
        domain_p : int
            local domain index
        state_p : ndarray
            PE-local global state vector
        state_l : ndarray
            local state vector for local analysis
        """
        # generic initialization
        # using id_lstate_in_pstate set in init_dim_l_pdaf
        state_l = state_p[self.id_lstate_in_pstate]

    def l2g_state_pdaf(self, step, domain_p, state_l, state_p):
        """convert local state vector to PE-local global state vector

        Parameters
        ----------
        step : int
            current time step
        domain_p : int
            local domain index
        state_l : ndarray
            local state vector for local analysis
        state_p : ndarray
            PE-local global state vector
        """
        state_p[self.id_lstate_in_pstate] = state_l
