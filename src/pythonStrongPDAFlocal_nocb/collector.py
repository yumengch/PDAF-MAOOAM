import numpy as np

import model
import state_vector


class collector:
    """This class implements functions where PDAF collects the state vector from model ensemble

    Attributes
    ----------
    model: model.model
        model instance
    sv: state_vector.state_vector
        state vector instance

    Methods
    -------
    collect_state_pdaf(dim_p:int, state_p:np.ndarray) -> np.ndarray
        generate state vector from the model fields
    init_ens_pdaf(filtertype:int, dim_p:int, dim_ens:int, state_p:np.ndarray, uinv:np.ndarray, ens_p:np.ndarray, status_pdaf:int) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]
        initialise the PDAF ensemble and nothing is done
    """
    def __init__(self, model_t:model.model, sv: state_vector.state_vector) -> None:
        # initialise the model instance
        self.model: model.model = model_t
        self.sv: state_vector.state_vector = sv

    def collect_state_pdaf(self, dim_p:int, state_p:np.ndarray) -> np.ndarray:
        """generate state vector from the model fields

        Parameters
        ----------
        dim_p : int
            size of state vector
        state_p : ndarray
            1D state vector on local PE

        Returns
        -------
        state_p : ndarray
            1D state vector on local PE
        """
        self.model.toPhysical_A()
        self.model.toPhysical_O()
        state_p[:] = np.concatenate([self.model.fields[varname].ravel(order='F')
                                     for varname in self.sv.varnames])
        return state_p

    def init_ens_pdaf(self, filtertype:int, dim_p:int, dim_ens:int, state_p:np.ndarray,
                      uinv:np.ndarray, ens_p:np.ndarray, status_pdaf:int) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
        """initialise the PDAF ensemble and nothing is done here

        Parameters
        ----------
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
            ensemble state vector on local PE (dim_p, dim_ens)
        status_pdaf : int
            status of PDAF

        Returns
        -------
        state_p : ndarray
            1D state vector on local PE
        uinv : ndarray
            2D left eigenvector with shape (dim_ens - 1,dim_ens - 1)
        ens_p : ndarray
            ensemble state vector on local PE (dim_p, dim_ens)
        status_pdaf : int
            status of PDAF
        """
        return state_p, uinv, ens_p, status_pdaf