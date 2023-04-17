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

import pyPDAF.PDAF.PDAFomi as PDAFomi
from Localization import Localization
import configparser

class Obs:

    """observation information and user-supplied routines

    Attributes
    ----------
    delt_obs : int
        time step interval for observations
    dim_obs : int
        dimension size of the observation vector
    dim_obs_p : int
        dimension size of the PE-local observation vector
    disttype : int
        type of distance computation to use for localization
    doassim : int
        whether to assimilate this observation type
    domainsize : ndarray
        size of domain for periodicity (<=0 for no periodicity)
    i_obs : int
        index of the observation type
    icoeff_p : ndarray
        2d array for interpolation coefficients for obs. operator
    id_obs_p : ndarray
        indices of process-local observed field in state vector
    ivar_obs_p : ndarray
        vector of process-local inverse observation error variance
    n_obs : int
        number of observation types
    ncoord : int
        number of coordinate dimension
    nrows : int
        number of rows in ocoord_p
    obs_err_type : int
        type of observation error
    obs_p : ndarray
        vector of process-local observations
    ocoord_p : ndarray
        2d array of process-local observation coordinates
    rms_obs : float
        observation error standard deviation (for constant errors)
    use_global_obs : int
       Whether to use (1) global full obs. or
       (0) obs. restricted to those relevant for a process domain

    local : `Localization.Localization`
        a localization info obejct
    filename : `str`
        observation filename
    filename_var : `str`
        observation variance filename
    varnames : list
        list of observation varnames
    file_timestep : int
        timestep for the observation files
    file_timecount : int
        current time step of the observation file
    
    """

    def __init__(self, i_obs, obsname, mype_filter, model_n):
        """constructor

        Parameters
        ----------
        i_obs : int
            index of the current observation
        obsname : string
            name of the observation type
        mype_filter : int
            rank of the PE in filter communicator
        model_n : float
            model domain aspect ratio parameter
        """

        self.i_obs = i_obs

        if (mype_filter == 0):
            print(('Assimilate observations:', obsname))

        config = configparser.ConfigParser()
        config.read(f'{obsname}.ini')
        
        self.local = Localization(config['local'])

        config = config['init']

        self.filename = config.get('filename', f'MAOOAM_{obsname}.nc')
        self.filename_var = config.get('filename_var', f'MAOOAM_{obsname}.nc')

        self.doassim = config.getint('doassim', 1)
        self.delt_obs = config.getint('delt_obs', 2)
        self.rms_obs = config.getfloat('rms_obs', 1.0)

        # Specify type of distance computation
        # 0=Cartesian 1=Cartesian periodic
        self.disttype = config.getint('disttype', 0)

        # Number of coordinates used for distance computation
        # The distance compution starts from the first row
        self.ncoord = config.getint('ncoord', 2)

        # Allocate process-local index array
        # This array has as many rows as required
        # for the observation operator
        # 1 if observations are at grid points;
        # >1 if interpolation is required
        self.nrows = config.getint('nrows', 1)

        # Size of domain for periodicity for disttype=1
        # (<0 for no periodicity)
        self.domainsize = np.zeros(self.ncoord)
        self.domainsize[0] = 2*np.pi/model_n
        self.domainsize[1] = np.pi

        # Type of observation error: (0) Gauss, (1) Laplace
        self.obs_err_type = config.getint('obs_err_type', 0)

        # Whether to use (1) global full obs.
        # (0) obs. restricted to those relevant for a process domain
        self.use_global_obs = config.getint('use_global_obs', 1)

        self.icoeff_p = None

        self.missing_value = -999

        self.file_timecount = -1
        self.file_timestep = config.getint('file_timestep', None)
        self.obs_den = config.getint('obs_den', None) # 8
        self.varnames = config.get('varnames', None).split(',')


    def init_dim_obs(self, step, dim_obs, model, isStrong,
                     mype_filter, dim_state, dim_state_p):
        """intialise PDAFomi and getting dimension of observation vector

        Parameters
        ----------
        step : int
            current time step
        dim_obs : int
            dimension size of the observation vector
        model : int
            rank of the PE in filter communicator
        dim_state : ndarray
            integer array for grid size
        dim_state_p : ndarray
            integer array for PE-local grid size
        """
        if isStrong:
            if (step % self.delt_obs == 0):
                self.doassim = 1
            else:
                self.doassim = 0
        
        if (self.doassim == 0):
            self.dim_obs =0
            return

        obs_field = self.get_obs_field(step)

        # Count valid observations that
        # lie within the process sub-domain
        pe_start = dim_state_p*mype_filter
        pe_end = dim_state_p*(mype_filter+1)
        obs_field_p = obs_field[pe_start:pe_end]
        cnt_p = np.count_nonzero(obs_field_p > self.missing_value)
        self.dim_obs_p = cnt_p

        # Initialize vector of observations on the process sub-domain
        # Initialize coordinate array of observations
        # on the process sub-domain
        if self.dim_obs_p > 0:
            self.set_obs_p(obs_field_p)
            self.set_id_obs_p(obs_field_p, model, isStrong)
            self.set_ocoord_p(obs_field_p, pe_start, model)
            self.set_ivar_obs_p(obs_field_p)
        else:
            self.obs_p = np.zeros(1)
            self.ivar_obs_p = np.zeros(1)
            self.ocoord_p = np.zeros((self.ncoord, 1))
            self.id_obs_p = np.zeros((self.nrows, 1))

        self.set_PDAFomi()


    def init_dim_obs_gen(self, step, dim_obs, model,
                         mype_filter, dim_state, dim_state_p):
        """intialise PDAFomi and getting dimension of observation vector

        Parameters
        ----------
        step : int
            current time step
        dim_obs : int
            dimension size of the observation vector
        mype_filter : int
            rank of the PE in filter communicator
        dim_state : ndarray
            integer array for grid size
        dim_state_p : ndarray
            integer array for PE-local grid size
        """

        # Count valid observations that
        # lie within the process sub-domain
        obs_field = np.zeros((len(self.varnames), model.ny, model.nx))[:, ::self.obs_den, :: self.obs_den]
        obs_field = obs_field.ravel()
        # Count valid observations that
        # lie within the process sub-domain
        pe_start = dim_state_p*mype_filter
        pe_end = dim_state_p*(mype_filter+1)
        obs_field_p = obs_field[pe_start:pe_end]
        cnt_p = np.count_nonzero(obs_field_p > self.missing_value)
        self.dim_obs_p = cnt_p

        # Initialize vector of observations on the process sub-domain
        # Initialize coordinate array of observations
        # on the process sub-domain
        if self.dim_obs_p > 0:
            self.set_obs_p(obs_field_p)
            self.set_id_obs_p(obs_field_p, model)
            self.set_ocoord_p(obs_field_p, pe_start, model)
            self.set_ivar_obs_p(obs_field_p)
        else:
            self.obs_p = np.zeros(1)
            self.ivar_obs_p = np.zeros(1)
            self.ocoord_p = np.zeros((self.ncoord, 1))
            self.id_obs_p = np.zeros((self.nrows, 1))

        self.set_PDAFomi()


    def set_obs_p(self, obs_field_p):
        """set up PE-local observation vector

        Parameters
        ----------
        obs_field_p : ndarray
            PE-local observation field
        """
        obs_field_tmp = obs_field_p.ravel()
        self.obs_p = np.zeros(self.dim_obs_p)
        self.obs_p[:self.dim_obs_p] = obs_field_tmp[obs_field_tmp > self.missing_value]


    def set_id_obs_p(self, obs_field_p, model, isStrong=True):
        """set id_obs_p

        Parameters
        ----------
        obs_field_p : ndarray
            PE-local observation field
        """
        nx = model.nx
        ny = model.ny
        offset = 0
        nobsvar = len(self.varnames)
        if ('psi_a' not in self.varnames) and isStrong: offset = 2*nx*ny
        self.id_obs_p = np.zeros((self.nrows, self.dim_obs_p))
        obs_field_tmp = obs_field_p.ravel()
        i = np.arange(nx)[::self.obs_den]
        j = np.arange(ny)[::self.obs_den]
        idx = np.arange(nx*ny, dtype=int).reshape(ny, nx)[j][:, i].ravel()
        self.id_obs_p[0, :self.dim_obs_p] = 1 + offset + np.concatenate([idx + i*nx*ny 
                                                                     for i in range(nobsvar)],
                                                                     dtype=int)[obs_field_tmp > self.missing_value]
 

    def set_ocoord_p(self, obs_field_p, offset, model):
        """set ocoord_p

        Parameters
        ----------
        obs_field_p : ndarray
            PE-local observation field
        offset : int
            PE-local offset starting from rank 0
        """
        self.ocoord_p = np.zeros((self.ncoord, self.dim_obs_p))
        nx = model.nx
        ny = model.ny
        dx = model.xc[0, 1] - model.xc[0, 0]
        dy = model.yc[1, 0] - model.yc[0, 0]
        obs_field_tmp = obs_field_p.ravel()
        nobs = len(self.varnames)    
        # assuming first ravel dx then ravel dy
        self.ocoord_p[0, :self.dim_obs_p] = np.tile(np.arange(nx)[::self.obs_den]*dx, 
                                                    (ny//self.obs_den + 1)*nobs
                                                    )[obs_field_tmp > self.missing_value]
        self.ocoord_p[1, :self.dim_obs_p] = np.repeat(np.arange(ny)[::self.obs_den]*dy, 
                                                      (nx//self.obs_den + 1)*nobs 
                                                      )[obs_field_tmp > self.missing_value]


    def set_ivar_obs_p(self, obs_field_p):
        """set ivar_obs_p
        """
        
        f = xr.open_dataset(self.filename_var, decode_times=False)
        var_obs = np.concatenate([f[varname+'_var'].to_numpy().ravel()
                                 for varname in self.varnames])
        f.close()
        obs_field_tmp = obs_field_p.ravel()
        self.ivar_obs_p = np.where(var_obs < 1e-12,
                                   1e14, 
                                   np.ones(
                                          self.dim_obs_p
                                          )/self.rms_obs/var_obs
                                   )


    def get_obs_field(self, step):
        """retrieve observation field

        Parameters
        ----------
        step : int
            current time step

        Returns
        -------
        obs_field : ndarray
            observation field
        """
        f = xr.open_dataset(self.filename, decode_times=False)
        self.file_timecount += self.file_timestep
        obs_field = np.concatenate([f[varname].isel(time=self.file_timecount).to_numpy().ravel() 
                              for varname in self.varnames])
        f.close()
        return obs_field


    def set_PDAFomi(self):
        """set PDAFomi obs_f object
        """
        PDAFomi.set_doassim(self.i_obs, self.doassim)
        PDAFomi.set_disttype(self.i_obs, self.disttype)
        PDAFomi.set_ncoord(self.i_obs, self.ncoord)
        PDAFomi.set_id_obs_p(self.i_obs, self.id_obs_p)
        if self.domainsize is not None:
            PDAFomi.set_domainsize(self.i_obs, self.domainsize)
        if self.obs_err_type is not None:
            PDAFomi.set_obs_err_type(self.i_obs, self.obs_err_type)
        if self.use_global_obs is not None:
            PDAFomi.set_use_global_obs(self.i_obs, self.use_global_obs)
        if self.icoeff_p is not None:
            PDAFomi.set_icoeff_p(self.i_obs, self.icoeff_p)

        self.dim_obs = PDAFomi.gather_obs(self.i_obs,
                                          self.obs_p,
                                          self.ivar_obs_p,
                                          self.ocoord_p,
                                          self.local.local_range)


    def obs_op_gridpoint(self, step, state_p, ostate):
        """convert state vector by observation operator

        Parameters
        ----------
        step : int
            current time step
        state_p : ndarray
            PE-local state vector
        ostate : ndarray
            state vector transformed by identity matrix
        """
        ostate = PDAFomi.obs_op_gridpoint(self.i_obs, state_p, ostate)
        return ostate


    def init_dim_obs_l(self, domain_p, step, dim_obs, dim_obs_l):
        """intialise local observation vector

        Parameters
        ----------
        local : TYPE
            Description
        domain_p : int
            index of current local analysis domain
        step : int
            current time step
        dim_obs : int
            dimension of observation vector
        dim_obs_l : int
            dimension of local observation vector

        Returns
        -------
        dim_obs_l : int
            dimension of local observations
        """
        return PDAFomi.init_dim_obs_l(self.i_obs, self.local.coords_l,
                                      self.local.loc_weight,
                                      self.local.local_range,
                                      self.local.srange, dim_obs_l)


    def localize_covar(self, HP_p, HPH, coords_p):
        """localze covariance matrix

        Parameters
        ----------
        local : `Localization.Localization`
            the localization object
        HP_p : ndarray
            matrix HPH
        HPH : ndarray
            PE local part of matrix HP
        coords_p : ndarray
            coordinates of state vector elements
        """
        return PDAFomi.localize_covar(self.i_obs, self.local.loc_weight,
                               self.local.local_range, self.local.srange,
                               coords_p, HP_p, HPH)


    def deallocate_obs(self):
        """deallocate PDAFomi object

        Parameters
        ----------
        step : int
            current time step
        """
        PDAFomi.deallocate_obs(self.i_obs)
