import numpy as np

import parallelisation
import model
import localisation
import obs_writer
import obs

class obs_factory:
    """This class implements all user-supplied functions
    used by PDAFomi. These functions are called at every time steps

    Attributes
    ----------
    pe : `parallelisation.parallelisation`
        parallelization object
    model : `model.model`
        model object
    local : `localisation.localisation`
        localisation object
    obs_list : list
        list of observation types
    nobs : int
        total number of observation types
    """
    def __init__(self, config, pe:parallelisation.parallelisation, model_t:model.model,
                 local: localisation.localisation) -> None:
        # Initialise observations
        self.pe: parallelisation.parallelisation =  pe
        self.model:model.model = model_t
        self.local:localisation.localisation = local
        # create new observations
        self.obs_list:list = []
        self.nobs:int = 0
        for key, item in config.items():
            if key != 'n_obs':
                self.nobs += 1
                self.obs_list.append(obs.obs(self.nobs, item, self.pe, self.model)
                )

    def setWriter(self) -> None:
        if self.pe.filter_pe:
            self.writers: list = []
            for i in range(self.nobs):
                self.writers.append( obs_writer.obs_writer(f'MAOOAM_{self.obs_list[i].obsname}.nc', self.model, self.obs_list[i].obs_den)
                )

    def set_doassim(self, step:int) -> None:
        for i in range(self.nobs):
            self.obs_list[i].doassim = 0
            if step % self.obs_list[i].delt_obs == 0:
                self.obs_list[i].doassim = 1

    @property
    def delt_obs(self):
        return min([obs.delt_obs for obs in self.obs_list])

    def get_obs_f(self, step:int, dim_obs_f:int, observation_f:np.ndarray) -> np.ndarray:
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

        istart:int = 0
        iend: int
        for i in range(self.nobs):
            iend = istart + self.obs_list[i].dim_obs
            self.writers[i].write(step, observation_f[istart:iend])
            istart = iend

        return observation_f

    def init_dim_obs_gen_pdafomi(self, step:int, dim_obs:int) -> int:
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
        for i in range(self.nobs):
            self.obs_list[i].init_dim_obs_gen(step, self.model)
            dim_obs += self.obs_list[i].dim_obs

        return dim_obs

    def init_dim_obs_pdafomi(self, step:int, dim_obs:int) -> int:
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
        for i in range(self.nobs):
            if self.obs_list[i].doassim == 1:
                self.obs_list[i].init_dim_obs(step, self.model)
                dim_obs += self.obs_list[i].dim_obs

        return dim_obs

    def obs_op_pdafomi(self, step:int, dim_p:int, dim_obs_p:int, state_p:np.ndarray, ostate:np.ndarray) -> np.ndarray:
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

        Returns
        -------
        ostate : ndarray
            state vector in obs space
        """
        for i in range(self.nobs):
            if self.obs_list[i].doassim == 1:
                ostate = self.obs_list[i].obs_op_gridpoint(step, state_p, ostate)

        return ostate

    def init_dim_obs_l_pdafomi(self, domain_p:int, step:int, dim_obs:int, dim_obs_l:int) -> int:
        """initialise number of observation for local domain

        Parameters
        ----------
        domain_p : int
            local domain index
        step : int
            current time step
        dim_obs : int
            dimension of observation vector
        dim_obs_l : int
            dimension of observation vector for local domain

        Retruns
        -------
        dim_obs_l : int
            dimension of observation vector for local domain
        """
        n_grid:int = self.model.nx*self.model.ny
        state_index:int = domain_p - ((domain_p - 1)//n_grid)*n_grid
        dx:float = self.model.xc[0, 1] - self.model.xc[0, 0]
        dy:float = self.model.yc[1, 0] - self.model.yc[0, 0]
        coords_l:np.ndarray = np.zeros(2)
        coords_l[0] = np.ceil(state_index/self.model.ny)
        coords_l[1] = state_index - (coords_l[0] - 1)*self.model.ny
        coords_l = (coords_l - 1)* np.array([dx, dy])
        for i in range(self.nobs):
            if self.obs_list[i].doassim == 1:
                dim_obs_l = self.obs_list[i].init_dim_obs_l(coords_l, dim_obs_l, self.local)
        return dim_obs_l
