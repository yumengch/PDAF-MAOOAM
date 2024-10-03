import configparser

import numpy as np
import pyPDAF_local.PDAF as PDAF # type: ignore

import log
import state_vector


class localisation:

    """class for localization information and user-supplied functions

    Attributes
    ----------
    loc_weight : int
        - (0) constant weight of 1
        - (1) exponentially decreasing with sradius
        - (2) use 5th-order polynomial
        - (3) regulated localization of R with mean error variance
        - (4) regulated localization of R with single-point error variance
    cradius : float
        range for local observation domain
    sradius : float
        support range for 5th order polynomial
        or radius for 1/e for exponential weighting
    local_filter : bool
        a boolean variable determining whether a local filter is used
        by default, it is false and will be later determined by PDAF function
        in init_pdaf
    sv : state_vector.state_vector
        a reference to the state vector object

    Methods
    -------
    init_n_domains(self, n_domains:int) -> int
        initialise the number of local domains
    set_lim_coords(self, xmin:float, xmax:float, ymin:float, ymax:float)
        set the coordinates of the local domain
    init_dim_l_pdaf(self, step:int, domain_p:int, dim_l:int) -> int
        initialise the local dimension of PDAF
    g2l_state_pdaf(self, step:int, domain_p:int, dim_p:int, state_p:np.ndarray, dim_l:int, state_l:np.ndarray) -> np.ndarray
        convert state vector (on each processor) to domain local state vector
    l2g_state_pdaf(self, step:int, domain_p:int, dim_l:int, state_l:np.ndarray, dim_p:int, state_p:np.ndarray) -> np.ndarray
        convert local state vector to PE-local global state vector
    """

    def __init__(self, sv:state_vector.state_vector, config:configparser.SectionProxy) -> None:
        """Initialise the localisation class.
        """
        self.locweight : int = config.getint('locweight', 0)
        self.cradius : float = config.getfloat('cradius', 0.0)
        self.sradius : float = config.getfloat('sradius', 0.0)
        # a boolean variable determining whether a local filter is used
        # by default, it is false and will be later determined by PDAF function
        # in init_pdaf
        self.local_filter:bool = False
        # a reference to the state vector object
        self.sv:state_vector.state_vector = sv

    def init_n_domains_pdaf(self, step:int, ndomains:int) -> int:
        """initialize the number of local domains

        Parameters
        ----------`
        step : int
            current time step
        ndomains: int
            number of local domains on current processor

        Returns
        -------
        n_domains_p : int
            PE-local number of analysis domains
        """
        log.logger.info(f'ndomains: ndomains {self.sv.dim_state_p}')
        return self.sv.dim_state_p

    def set_lim_coords(self, xmin:float, xmax:float, ymin:float, ymax:float) -> None:
        """Set the coordinates of the local domain.
        """
        lim_coords:np.ndarray = np.zeros((2, 2))
        lim_coords[0, 0] = xmin
        lim_coords[0, 1] = xmax
        lim_coords[1, 0] = ymax
        lim_coords[1, 1] = ymin

        PDAF.omi_set_domain_limits(lim_coords)

    def init_dim_l_pdaf(self, step:int, domain_p:int, dim_l:int) -> int:
        """initialise the local dimension of PDAF.

        The function returns
        the dimension of local state vector

        Parameters
        ----------
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
        id_lstate_in_pstate:np.ndarray = domain_p*np.ones(dim_l, dtype=np.intc)
        PDAF.local_set_indices(id_lstate_in_pstate)
        return dim_l

    def g2l_state_pdaf(self, step:int, domain_p:int, dim_p:int, state_p:np.ndarray, dim_l:int, state_l:np.ndarray) -> np.ndarray:
        """convert state vector (on each processor) to domain local state vector

        Parameters
        ----------
        step : int
            current time step
        domain_p : int
            local domain index
        dim_p : int
            dimension of state vector
        state_p : ndarray
            state vector
        dim_l : int
            dimension of domain local state vector
        state_l : ndarray
            domain local state vector for local analysis
        """
        # generic initialization
        # domain_p is the *domain_p*-th local domain on the local processor
        # The dimension of the state vector on the local domain is dim_l defined in init_dim_l_pdaf
        # Here, because we set dim_l = 1, each local domain is one element of the state vector
        # In more complex cases, it is possible to define a relationship matrix between local domain
        # and the element of the state vector.
        state_l[0] = state_p[domain_p - 1]
        return state_l

    def l2g_state_pdaf(self, step:int, domain_p:int, dim_l:int, state_l:np.ndarray, dim_p:int, state_p:np.ndarray) -> np.ndarray:
        """convert local state vector to PE-local global state vector

        Parameters
        ----------
        step : int
            current time step
        domain_p : int
            local domain index
        dim_l : int
            dimension of local state vector
        state_l : ndarray
            local state vector for local analysis
        dim_p : int
            dimension of state vector
        state_p : ndarray
            state vector
        """
        state_p[domain_p -1] = state_l[0]
        return state_p
