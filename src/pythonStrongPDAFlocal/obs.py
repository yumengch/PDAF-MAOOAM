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

import pyPDAF_local.PDAF as PDAF # type: ignore
import configparser

import log
import model
import localisation
import parallelisation

class obs:

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

    def __init__(self, i_obs:int, obsname:str, pe:parallelisation.parallelisation,
                 model_t:model.model) -> None:
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

        self.i_obs:int = i_obs
        self.obsname:str = obsname
        if (pe.mype_filter == 0):
            log.logger.info(('Assimilate observations:', obsname))

        config:configparser.ConfigParser = configparser.ConfigParser()
        config.read(f'{obsname}.ini')

        config_init:configparser.SectionProxy = config['init']

        self.filename:str = config_init.get('filename', f'MAOOAM_{obsname}.nc')
        self.filename_var:str = config_init.get('filename_var', f'MAOOAM_{obsname}.nc')

        self.doassim:int = config_init.getint('doassim', 1)
        self.delt_obs:int = config_init.getint('delt_obs', 2)
        self.rms_obs:float = config_init.getfloat('rms_obs', 1.0)

        # Specify type of distance computation
        # 0=Cartesian 1=Cartesian periodic
        self.disttype:int = config_init.getint('disttype', 0)

        # Number of coordinates used for distance computation
        # The distance compution starts from the first row
        self.ncoord:int = config_init.getint('ncoord', 2)

        # Allocate process-local index array
        # This array has as many rows as required
        # for the observation operator
        # 1 if observations are at grid points;
        # >1 if interpolation is required
        self.nrows:int = config_init.getint('nrows', 1)

        # Size of domain for periodicity for disttype=1
        # (<0 for no periodicity)
        self.domainsize:np.ndarray = np.zeros(self.ncoord)
        self.domainsize[0] = model_t.xc[0, -1]
        self.domainsize[1] = model_t.yc[-1, 0]

        # Type of observation error: (0) Gauss, (1) Laplace
        self.obs_err_type:int = config_init.getint('obs_err_type', 0)

        # Whether to use (1) global full obs.
        # (0) obs. restricted to those relevant for a process domain
        self.use_global_obs:int = config_init.getint('use_global_obs', 1)

        self.icoeff_p = None

        self.missing_value = -999

        self.file_timecount:int = -1
        self.file_timestep:int = config_init.getint('file_timestep', 1)
        self.obs_den:int = config_init.getint('obs_den', 1) # 8
        self.varnames:list[str] = config_init.get('varnames', None).split(',')


    def init_dim_obs(self, step:int, model_t:model.model) -> None:
        """intialise PDAFomi and getting dimension of observation vector

        Parameters
        ----------
        step : int
            current time step
        model : int
            rank of the PE in filter communicator
        """
        if (self.doassim == 0):
            self.dim_obs =0
            return

        print ('PDAFomi: Obs varname: ', self.varnames)
        obs_field = self.get_obs_field()

        # Count valid observations that
        # lie within the process sub-domain
        obs_field_p = obs_field
        cnt_p = np.count_nonzero(obs_field_p > self.missing_value)
        self.dim_obs_p = cnt_p

        # Initialize vector of observations on the process sub-domain
        # Initialize coordinate array of observations
        # on the process sub-domain
        if self.dim_obs_p > 0:
            self.set_obs_p(obs_field_p)
            self.set_id_obs_p(obs_field_p, model_t)
            self.set_ocoord_p(obs_field_p, model_t)
            self.set_ivar_obs_p()
        else:
            self.obs_p = np.zeros(1)
            self.ivar_obs_p = np.zeros(1)
            self.ocoord_p = np.zeros((self.ncoord, 1))
            self.id_obs_p = np.zeros((self.nrows, 1))

        self.set_PDAFomi()


    def init_dim_obs_gen(self, model_t:model.model) -> None:
        """intialise PDAFomi and getting dimension of observation vector

        Parameters
        ----------
        model : int
            rank of the PE in filter communicator
        """

        # Count valid observations that
        # lie within the process sub-domain
        obs_field:np.ndarray = np.zeros((len(self.varnames), model_t.ny, model_t.nx))[:, ::self.obs_den, :: self.obs_den]
        obs_field_p:np.ndarray = obs_field.ravel()
        self.dim_obs_p = len(obs_field_p)

        # Initialize vector of observations on the process sub-domain
        # Initialize coordinate array of observations
        # on the process sub-domain
        if self.dim_obs_p > 0:
            self.set_obs_p(obs_field_p)
            self.set_id_obs_p(obs_field_p, model_t)
            self.set_ocoord_p(obs_field_p, model_t)
            self.set_ivar_obs_p()
        else:
            self.obs_p = np.zeros(1)
            self.ivar_obs_p = np.zeros(1)
            self.ocoord_p = np.zeros((self.ncoord, 1))
            self.id_obs_p = np.zeros((self.nrows, 1))

        self.set_PDAFomi()


    def set_obs_p(self, obs_field_p: np.ndarray) -> None:
        """set up PE-local observation vector

        Parameters
        ----------
        obs_field_p : ndarray
            PE-local observation field
        """
        obs_field_tmp:np.ndarray = obs_field_p.ravel()
        self.obs_p = np.zeros(self.dim_obs_p)
        self.obs_p[:self.dim_obs_p] = obs_field_tmp[obs_field_tmp > self.missing_value]

    def set_id_obs_p(self, obs_field_p: np.ndarray, model_t:model.model) -> None:
        """set id_obs_p

        Parameters
        ----------
        obs_field_p : ndarray
            PE-local observation field
        """
        nx:int = model_t.nx
        ny:int = model_t.ny
        nobsvar:int = len(self.varnames)
        offset:int = 0
        if 'psi_a' not in self.varnames: offset = 2*nx*ny
        self.id_obs_p = np.zeros((self.nrows, self.dim_obs_p), dtype=np.intc)
        obs_field_tmp:np.ndarray = obs_field_p.ravel()
        i:np.ndarray = np.arange(nx)[::self.obs_den]
        j:np.ndarray = np.arange(ny)[::self.obs_den]
        idx:np.ndarray = np.arange(nx*ny, dtype=int).reshape(ny, nx)[j][:, i].ravel()
        self.id_obs_p[0, :self.dim_obs_p] = 1 + offset + np.concatenate([idx + i*nx*ny
                                                                     for i in range(nobsvar)],
                                                                     dtype=int)[obs_field_tmp > self.missing_value]

    def set_ocoord_p(self, obs_field_p: np.ndarray, model_t: model.model) -> None:
        """set ocoord_p

        Parameters
        ----------
        obs_field_p : ndarray
            PE-local observation field
        """
        self.ocoord_p = np.zeros((self.ncoord, self.dim_obs_p))
        nx:int = model_t.nx
        ny:int = model_t.ny
        dx:float = model_t.xc[0, 1] - model_t.xc[0, 0]
        dy:float = model_t.yc[1, 0] - model_t.yc[0, 0]
        obs_field_tmp: np.ndarray = obs_field_p.ravel()
        nobs:int = len(self.varnames)
        # assuming first ravel dx then ravel dy
        self.ocoord_p[0, :self.dim_obs_p] = np.tile(np.arange(nx)[::self.obs_den]*dx,
                                                    (ny//self.obs_den + 1)*nobs
                                                    )[obs_field_tmp > self.missing_value]
        self.ocoord_p[1, :self.dim_obs_p] = np.tile(np.repeat(
                                                              np.arange(ny)[::self.obs_den]*dy,
                                                              (nx//self.obs_den + 1)
                                                              ),
                                                    nobs)[obs_field_tmp > self.missing_value]

    def set_ivar_obs_p(self) -> None:
        """set ivar_obs_p
        """

        f:xr.Dataset = xr.open_dataset(self.filename_var, decode_times=False)
        var_obs:np.ndarray = np.concatenate([f[varname+'_var'].to_numpy().ravel()
                                 for varname in self.varnames])
        f.close()

        self.ivar_obs_p = np.where(var_obs < 1e-12,
                                   1e14,
                                   np.ones(
                                          self.dim_obs_p
                                          )/self.rms_obs/var_obs
                                   )

    def get_obs_field(self) -> np.ndarray:
        """retrieve observation field

        Returns
        -------
        obs_field : ndarray
            observation field
        """
        f:xr.Dataset = xr.open_dataset(self.filename, decode_times=False)
        self.file_timecount += self.file_timestep
        obs_field:np.ndarray = np.concatenate([f[varname].isel(time=self.file_timecount).to_numpy().ravel()
                                              for varname in self.varnames])
        f.close()
        return obs_field

    def set_PDAFomi(self) -> None:
        """set PDAFomi obs_f object
        """
        PDAF.omi_set_doassim(self.i_obs, self.doassim)
        PDAF.omi_set_disttype(self.i_obs, self.disttype)
        PDAF.omi_set_ncoord(self.i_obs, self.ncoord)
        PDAF.omi_set_id_obs_p(self.i_obs, self.id_obs_p)
        if self.domainsize is not None:
            PDAF.omi_set_domainsize(self.i_obs, self.domainsize)
        if self.obs_err_type is not None:
            PDAF.omi_set_obs_err_type(self.i_obs, self.obs_err_type)
        if self.use_global_obs is not None:
            PDAF.omi_set_use_global_obs(self.i_obs, self.use_global_obs)
        if self.icoeff_p is not None:
            PDAF.omi_set_icoeff_p(self.i_obs, self.icoeff_p)

        self.dim_obs = PDAF.omi_gather_obs(self.i_obs,
                                          self.obs_p,
                                          self.ivar_obs_p,
                                          self.ocoord_p,
                                          0.)

    def obs_op_gridpoint(self, step : int, state_p: np.ndarray, ostate : np.ndarray) -> np.ndarray:
        """convert state vector by observation operator

        Parameters
        ----------
        step : int
            current time step
        state_p : ndarray
            PE-local state vector
        ostate : ndarray
            state vector transformed by identity matrix

        Returns
        -------
        ostate : ndarray
            state vector transformed by observation operator
        """
        ostate = PDAF.omi_obs_op_gridpoint(self.i_obs, state_p, ostate)
        return ostate

    def init_dim_obs_l(self, coords_l: np.ndarray, dim_obs_l : int, local:localisation.localisation) -> int:
        """initialise dimensions of observations for local domain

        Parameters
        ----------
        coords_l : ndarray
            local domain coordinates
        dim_obs_l : int
            dimension of observation vector for local domain

        Returns
        -------
            dim_obs_l : int
                dimension of observation vector for local domain
        """
        dim_obs_l = PDAF.omi_init_dim_obs_l_iso(self.i_obs, coords_l,
                                    local.locweight,
                                    local.cradius,
                                    local.sradius, dim_obs_l)
        return dim_obs_l
