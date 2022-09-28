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
import functools
import U_PDAF
import U_PDAFomi


class init_pdaf:

    """initialise PAF

    Attributes
    ----------
    filter_param_i : ndarray
        a list of integer options for EnKF
    filter_param_r : ndarray
        a list of float options for EnKF
    """

    def __init__(self, sv, infl, options,
                 local, model, pe, obs, writer, screen):
        """constructor

        Parameters
        ----------
        sv : `StateVector.StateVector`
            an object of StateVector
        infl : `Inflation.Inflation`
            inflation object
        options : `FilterOptions.FilterOptions`
            filtering options
        local : `local.local`
            local object
        model : `Model.Model`
            model object
        pe : `parallelization.parallelization`
            parallelization object
        obs : `ObsFactory.ObsFactory`
            ObsFactory object
        screen : int
            verbosity of PDAF screen output
        """
        if (options.filtertype == 2):
            # EnKF with Monte Carlo init
            filter_param_i, filter_param_r = \
                self.setEnKFOptions(6, 2, sv, infl, options)
        else:
            # All other filters
            filter_param_i, filter_param_r = \
                self.setETKFOptions(7, 2, sv, infl, options)

        U_init_ens_pdaf = functools.partial(U_PDAF.init_ens_pdaf,
                                            model)

        _, _, status = PDAF.init(options.filtertype,
                                 options.subtype,
                                 0,
                                 filter_param_i,
                                 filter_param_r,
                                 pe.COMM_model.py2f(),
                                 pe.COMM_filter.py2f(),
                                 pe.COMM_couple.py2f(), pe.task_id,
                                 pe.n_modeltasks, pe.filterpe,
                                 U_init_ens_pdaf, screen)
        try:
            assert status == 0, \
                f'ERROR {status} \
                in initialization of PDAF - stopping! \
                (PE f{pe.mype_world})'
        except AssertionError:
            pe.abort_parallel()

        U_next_observation_pdaf = \
            functools.partial(U_PDAF.next_observation_pdaf,
                              model, pe, obs.delt_obs)
        U_distribute_state_pdaf = \
            functools.partial(U_PDAF.distribute_state_pdaf,
                              model)
        U_prepoststep_ens_pdaf = \
            functools.partial(U_PDAF.prepoststep_ens_pdaf,
                              sv, model, pe, obs, writer)

        steps, time, doexit, status = PDAF.get_state(10, 10,
                                         U_next_observation_pdaf,
                                         U_distribute_state_pdaf,
                                         U_prepoststep_ens_pdaf,
                                         status)

    def setEnKFOptions(self, dim_pint, dim_preal,
                       sv, infl, options):
        """set ensemble kalman filter options

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
        filter_param_i[2] = options.rank_analysis_enkf
        filter_param_i[3] = options.incremental
        filter_param_i[4] = 0

        filter_param_r[0] = infl.forget

        return filter_param_i, filter_param_r

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

    def __init__(self, model, obs, pe, sv, local, writer, filtertype):
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
        local : `Localization.Localization`
            a localization object
        filtertype : int
            type of filter
        """
        localfilter = PDAF.get_localfilter()

        U_collect_state_pdaf = \
            functools.partial(U_PDAF.collect_state_pdaf,
                              model)
        U_next_observation_pdaf = \
            functools.partial(U_PDAF.next_observation_pdaf,
                              model, pe, obs.delt_obs)
        U_distribute_state_pdaf = \
            functools.partial(U_PDAF.distribute_state_pdaf,
                              model)
        U_prepoststep_ens_pdaf = \
            functools.partial(U_PDAF.prepoststep_ens_pdaf,
                              sv, model, pe, obs, writer)
        U_init_dim_obs_PDAFomi = \
            functools.partial(U_PDAFomi.init_dim_obs_pdafomi,
                              obs,
                              local.local_range,
                              pe.mype_filter,
                              sv.dim_state,
                              sv.dim_state_p)

        U_obs_op_PDAFomi = \
            functools.partial(U_PDAFomi.obs_op_pdafomi,
                              obs)

        if (localfilter == 1):
            U_init_n_domains_pdaf = \
                functools.partial(local.init_n_domains_pdaf,
                                  sv)
            U_init_dim_l_pdaf = \
                functools.partial(local.init_dim_l_pdaf,
                                  model, pe.mype_filter)
            U_init_dim_obs_l_pdafomi = \
                functools.partial(U_PDAFomi.init_dim_obs_l_pdafomi,
                                  obs, local)
            status = \
                PDAFomi.assimilate_local(U_collect_state_pdaf,
                                              U_distribute_state_pdaf,
                                              U_init_dim_obs_PDAFomi,
                                              U_obs_op_PDAFomi,
                                              U_prepoststep_ens_pdaf,
                                              U_init_n_domains_pdaf,
                                              U_init_dim_l_pdaf,
                                              U_init_dim_obs_l_pdafomi,
                                              local.g2l_state_pdaf,
                                              local.l2g_state_pdaf,
                                              U_next_observation_pdaf)
        else:
            if filtertype == 8:
                U_localize_covar_pdafomi = \
                    functools.partial(U_PDAFomi.localize_covar_pdafomi,
                                      obs, local, model,
                                      pe.mype_filter)
                status = \
                    PDAFomi.assimilate_lenkf(U_collect_state_pdaf,
                                             U_distribute_state_pdaf,
                                             U_init_dim_obs_PDAFomi,
                                             U_obs_op_PDAFomi,
                                             U_prepoststep_ens_pdaf,
                                             U_localize_covar_pdafomi,
                                             U_next_observation_pdaf)
            elif filtertype == 100: 
                U_get_obs_f = functools.partial(U_PDAFomi.get_obs_f, obs)
                U_init_dim_obs_gen_PDAFomi = \
                    functools.partial(U_PDAFomi.init_dim_obs_gen_pdafomi,
                                      obs,
                                      local.local_range,
                                      pe.mype_filter,
                                      sv.dim_state,
                                      sv.dim_state_p)
                status = \
                    PDAFomi.generate_obs(U_collect_state_pdaf,
                                         U_distribute_state_pdaf,
                                         U_init_dim_obs_gen_PDAFomi,
                                         U_obs_op_PDAFomi,
                                         U_get_obs_f,
                                         U_prepoststep_ens_pdaf,
                                         U_next_observation_pdaf)
            else:
                status = \
                    PDAFomi.assimilate_global(U_collect_state_pdaf,
                                              U_distribute_state_pdaf,
                                              U_init_dim_obs_PDAFomi,
                                              U_obs_op_PDAFomi,
                                              U_prepoststep_ens_pdaf,
                                              U_next_observation_pdaf)


        if status != 0:
            print(('ERROR ', status,
                   ' in PDAF_put_state - stopping! (PE ',
                   pe.mype_world, ')'))
            pe.abort_parallel()
