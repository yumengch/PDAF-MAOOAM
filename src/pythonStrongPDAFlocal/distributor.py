import log
import model
import parallelisation
import obs_factory
import state_vector

import numpy as np

class distributor:
    """This class implements the function where
    PDAF distributes ensemble to the model field

    Attributes
    ----------
    model: model.model
        model instance
    pe: parallelisation.parallelisation
        parallelisation instance
    sv: state_vector.state_vector
        state vector instance
    obs: obs_factory.obs_factory
        observation factory instance

    Methods
    -------
    distribute_initial_state(dim_p:int, state_p:np.ndarray) -> np.ndarray
        Distribute initial ensemble to model field, nothing is done
    distribute_atmosphere_state(dim_p:int, state_p:np.ndarray) -> np.ndarray
        Distribute atmosphere component of analysis state to model field
    distribute_ocean_and_atmosphere_state(dim_p:int, state_p:np.ndarray) -> np.ndarray
        Distribute full analysis state vector to model field
    next_observation_pdaf(stepnow:int, nsteps:int, doexit:int, time:float) -> tuple[int, int, float]
        The time for the next observation
    """
    def __init__(self, model_t:model.model,
                 pe:parallelisation.parallelisation,
                 obs:obs_factory.obs_factory,
                 sv:state_vector.state_vector) -> None:
        # get the model insta
        self.model:model.model = model_t
        self.pe:parallelisation.parallelisation = pe
        self.sv: state_vector.state_vector = sv
        self.obs: obs_factory.obs_factory = obs

    def distribute_initial_state(self, dim_p:int, state_p:np.ndarray) -> np.ndarray:
        """Distribute initial ensemble to model field, nothing is done
        """
        if (self.pe.mype_ens == 0):
            log.logger.info ('distribute_state_pdaf: starting from restart files')
        return state_p

    def distribute_atmosphere_state(self, dim_p:int, state_p:np.ndarray) -> np.ndarray:
        """Distribute atmosphere component of analysis state to model field

        Parameters
        ----------
        dim_p : int
            size of the state vector
        state_p : ndarray
            1D state vector on local PE
        """
        n_grid = self.model.nx*self.model.ny
        varnames = self.sv.varnames[:].copy()
        for i, varname in enumerate(varnames):
            self.model.fields[varname][:] = \
                state_p[i*n_grid:(i+1)*n_grid].reshape(self.model.nx,
                                                   self.model.ny, order='F')
        if self.pe.mype_ens == 0: log.logger.info ('distbute atmos')
        self.model.toFourier_A()
        return state_p

    def distribute_ocean_and_atmosphere_state(self, dim_p:int, state_p:np.ndarray) -> np.ndarray:
        """Distribute full analysis state vector to model field

        Parameters
        ----------
        dim_p : int
            size of the state vector
        state_p : ndarray
            1D state vector on local PE
        """
        n_grid = self.model.nx*self.model.ny
        varnames = self.sv.varnames[:].copy()
        for i, varname in enumerate(varnames):
            self.model.fields[varname][:] = \
                state_p[i*n_grid:(i+1)*n_grid].reshape(self.model.nx,
                                                   self.model.ny, order='F')
        if self.pe.mype_ens == 0: log.logger.info ('distbute atmos')
        self.model.toFourier_A()
        if self.pe.mype_ens == 0: log.logger.info ('distbute ocean')
        self.model.toFourier_O()
        return state_p

    def next_observation_pdaf(self, stepnow:int, nsteps:int, doexit:int, time:float) -> tuple[int, int, float]:
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
        nsteps = self.obs.delt_obs
        if stepnow + nsteps <= self.model.total_steps:
            doexit = 0
            if self.pe.mype_ens == 0:
                log.logger.info(f'{stepnow} {stepnow} Next observation at time step {stepnow + nsteps}')
        else:
            doexit = 1
            if self.pe.mype_ens == 0:
                log.logger.info(f'{stepnow} No more observations - end assimilation')

        return nsteps, doexit, time