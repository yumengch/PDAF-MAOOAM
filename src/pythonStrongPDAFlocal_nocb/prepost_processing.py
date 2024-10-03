import numpy as np
import parallelisation
import log

import numba # type: ignore

@numba.njit
def rms(ens_p):
    n_samples, n_ensembles = ens_p.shape

    # Calculate the variance for each sample
    variances = np.empty(n_samples)
    for i in range(n_samples):
        variances[i] = np.var(ens_p[i, :])*n_ensembles/(n_ensembles-1)

    # Calculate the mean of the variances and then the RMS error estimate
    mean_variance = np.mean(variances)
    rmserror_est = np.sqrt(mean_variance)

    return rmserror_est

class prepost_processor:
    """User-supplied functions for pre and post processing of the ensemble.

    Attributes
    ----------
    model : `model.model`
        model object
    pe : `parallelisation.parallelisation`
        parallelisation object

    Methods
    -------
    initial_process(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, state_p, uinv, ens_p, flag)
        pre- and post-processing of ensemble
    prepostprocess(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, state_p, uinv, ens_p, flag)
        pre- and post-processing of ensemble
    """
    def __init__(self, pe:parallelisation.parallelisation) -> None:
        self.pe:parallelisation.parallelisation = pe

    def initial_process(self, step, dim_p, dim_ens, dim_ens_p,
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
        log.logger.info('Analyze initial state ensemble')

        rmserror_est = np.sqrt(np.mean(np.var(ens_p, axis=1, ddof=1)))

        if self.pe.mype_filter == 0:
            log.logger.info(f'RMS error: {rmserror_est}')

        return state_p, uinv, ens_p

    def prepostprocess(self, step, dim_p, dim_ens, dim_ens_p,
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
        if (step < 0):
            log.logger.info('Analyze and write forecasted state ensemble')
        else:
            log.logger.info('Analyze and write assimilated state ensemble')

        if self.pe.mype_filter == 0:
            log.logger.info(f'RMS error: {rms(ens_p)}')

        return state_p, uinv, ens_p