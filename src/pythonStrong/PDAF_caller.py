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
import pyPDAF.PDAF as PDAF
import pyPDAF.PDAF.PDAFomi as PDAFomi

class init_pdaf:

    """initialise PAF

    Attributes
    ----------
    filter_param_i : ndarray
        a list of integer options for EnKF
    filter_param_r : ndarray
        a list of float options for EnKF
    """

    def __init__(self, das):
        """constructor

        Parameters
        ----------
        das: `DAS.DAS`
            object for the DA system
        screen : int
            verbosity of PDAF screen output
        """
        # All other filters
        filter_param_i, filter_param_r = \
            self.setETKFOptions(7, 2, das.sv, das.infl, das.options)

        _, _, status = PDAF.init(das.options.filtertype,
                                 das.options.subtype,
                                 0,
                                 filter_param_i,
                                 filter_param_r,
                                 das.pe.COMM_model.py2f(),
                                 das.pe.COMM_filter.py2f(),
                                 das.pe.COMM_couple.py2f(), das.pe.task_id,
                                 das.pe.n_modeltasks, das.pe.filterpe,
                                 das.UserFuncs.init_ens_pdaf, das.screen)
        try:
            assert status == 0, \
                f'ERROR {status} \
                in initialization of PDAF - stopping! \
                (PE f{pe.mype_world})'
        except AssertionError:
            pe.abort_parallel()

        steps, time, doexit, status = PDAF.get_state(10, 10,
                                         das.UserFuncs.next_observation_pdaf,
                                         das.UserFuncs.distribute_state_pdaf,
                                         das.UserFuncs.prepoststep_ens_pdaf,
                                         status)

    def setETKFOptions(self, dim_pint, dim_preal,
                       sv, infl, options):
        """Summary

        Parameters
        ----------
        dim_pint : int
            size of integer filter options
        dim_preal : int
            size of float filter options
        sv : `StateVector.StateVector`
            an object of StateVector
        infl : `Inflation.Inflation`
            inflation object
        options : `FilterOptions.FilterOptions`
            filtering options
        """
        filter_param_i = np.zeros(dim_pint, dtype=int)
        filter_param_r = np.zeros(dim_preal)

        filter_param_i[0] = sv.dim_state_p
        filter_param_i[1] = sv.dim_ens
        filter_param_i[2] = 0
        filter_param_i[3] = options.incremental
        filter_param_i[4] = infl.type_forget
        filter_param_i[5] = options.type_trans
        filter_param_i[6] = options.type_sqrt

        filter_param_r[0] = infl.forget

        return filter_param_i, filter_param_r


class assimilate_pdaf:

    """assimilation calls
    """

    def __init__(self, das, it):
        """constructor

        Parameters
        ----------
        model : `Model.Model`
            model object
        obs : `ObsFactory.ObsFactory`
            observation object
        pe : `parallelization.parallelization`
            parallelization object
        sv : `StateVector.StateVector`
            an object of StateVector
        filtertype : int
            type of filter
        """


        if das.options.filtertype == 100:
            status = \
                PDAFomi.generate_obs(das.UserFuncs.collect_state_pdaf,
                                     das.UserFuncs.distribute_state_pdaf,
                                     das.UserFuncsO.init_dim_obs_gen_pdafomi,
                                     das.UserFuncsO.obs_op_pdafomi,
                                     das.UserFuncsO.get_obs_f,
                                     das.UserFuncs.prepoststep_ens_pdaf,
                                     das.UserFuncs.next_observation_pdaf)
        else:
            das.obs.set_doassim(it, das.sv)
            status = \
                PDAFomi.assimilate_global(das.UserFuncs.collect_state_pdaf,
                                          das.UserFuncs.distribute_state_pdaf,
                                          das.UserFuncsO.init_dim_obs_pdafomi,
                                          das.UserFuncsO.obs_op_pdafomi,
                                          das.UserFuncs.prepoststep_ens_pdaf,
                                          das.UserFuncs.next_observation_pdaf)

        if status != 0:
            print(('ERROR ', status,
                   ' in PDAF_put_state - stopping! (PE ',
                   pe.mype_world, ')'))
            pe.abort_parallel()
