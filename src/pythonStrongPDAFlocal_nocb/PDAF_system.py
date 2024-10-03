
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
import pyPDAF_local.PDAF as PDAF

import config
import log
from collector import collector
from distributor import distributor
from prepost_processing import prepost_processor
from parallelisation import parallelisation
from obs_factory import obs_factory
from state_vector import state_vector
from filter_options import filter_options
from model import model
from localisation import localisation


class DAS:
    """Data assimilation system

    Attributes
    ----------
    pe : `parallelization.parallelization`
        parallelization object
    sv : `StateVector.StateVector`
        an object of StateVector
    options : `FilterOptions.FilterOptions`
        filtering options
    local : `Localization.Localization`
        localization object
    model : `Model.Model`
        model object
    obs : `ObsFactory.ObsFactory`
        a factory of observations

    screen : int
        verbosity of PDAF screen output
    isStrong : bool
        Whether strong coupling is used
    """

    def __init__(self):
        """constructor
        """

    def init(self, config_t:config.PDAFConfig):
        """initialise DA system

        Parameters
        ----------
        config_t : `config.PDAFConfig`
            configuration object
        """
        # init DAS options
        self.screen = config_t['Global'].getint('screen', 3)
        # init parallelisation
        self.pe:parallelisation = parallelisation(config_t['Ensemble'].getint('dim_ens', 1),
                                                  config_t['Ensemble'].getint('n_modeltasks', 1))
        # init model
        self.model:model = model(config_t['Model'], self.pe.task_id)
        # init options
        self.sv:state_vector = state_vector(self.model, self.pe.dim_ens)
        self.options:filter_options = filter_options(config_t['FilterOptions'])
        self.local:localisation = localisation(self.sv, config_t['Local'])
        # init obs
        self.obs:obs_factory = obs_factory(config_t['Obs'], self.pe, self.model, self.local)

        # init obs. output
        if self.options.filtertype == 100:
            self.obs.setWriter()

        # Initial Screen output
        if (self.pe.mype_ens == 0):
            log.logger.info('+++++ pyPDAF online mode +++++')
            log.logger.info('+++++ MAOOAM +++++')


    def init_pdaf(self) -> None:
        """initialise PDAF
        """
        # All other filters
        filter_param_i, filter_param_r = \
            self.setETKFOptions(7, 2, self.sv, self.options)

        colltor:collector = collector(self.model, self.sv)
        _, _, status = PDAF.init(self.options.filtertype,
                                 self.options.subtype,
                                 0,
                                 filter_param_i,
                                 filter_param_r,
                                 self.pe.comm_model.py2f(),
                                 self.pe.comm_filter.py2f(),
                                 self.pe.comm_couple.py2f(), self.pe.task_id,
                                 self.pe.n_modeltasks, self.pe.filter_pe,
                                 colltor.init_ens_pdaf, self.screen)
        assert status == 0, f'ERROR {status} in initialization of PDAF - stopping! (PE f{self.pe.mype_ens})'

        lfilter:int = PDAF.get_localfilter()
        self.local.local_filter = lfilter == 1

        PDAF.omi_init(self.obs.nobs)

        distribtor:distributor = distributor(self.model, self.pe, self.obs, self.sv)
        init_processor: prepost_processor = prepost_processor(self.pe)
        _, _, _, status = PDAF.get_state(10, 10,
                                         distribtor.next_observation_pdaf,
                                         distribtor.distribute_initial_state,
                                         init_processor.initial_process,
                                         status)

        self.local.set_lim_coords(self.model.xc[0, 0], self.model.xc[0, -1],
                                  self.model.yc[0, 0], self.model.yc[-1, 0])


    def setETKFOptions(self, dim_pint:int, dim_preal:int,
                       sv:state_vector, options:filter_options) -> tuple[np.ndarray, np.ndarray]:
        """Summary

        Parameters
        ----------
        dim_pint : int
            size of integer filter options
        dim_preal : int
            size of float filter options
        sv : `StateVector.StateVector`
            an object of StateVector
        options : `FilterOptions.FilterOptions`
            filtering options
        """
        filter_param_i = np.zeros(dim_pint, dtype=np.intc)
        filter_param_r = np.zeros(dim_preal)

        filter_param_i[0] = sv.dim_state_p
        filter_param_i[1] = sv.dim_ens
        filter_param_i[2] = 0
        filter_param_i[3] = options.incremental
        filter_param_i[4] = options.type_forget
        filter_param_i[5] = options.type_trans
        filter_param_i[6] = options.type_sqrt

        filter_param_r[0] = options.forget

        return filter_param_i, filter_param_r


    def assimilate_pdaf(self, it:int) -> None:
        """assimilation call

        Parameters
        ----------
        it : int
            current time step
        """
        colltor:collector = collector(self.model, self.sv)
        distribtor:distributor = distributor(self.model, self.pe, self.obs, self.sv)
        ensemble_processor: prepost_processor = prepost_processor(self.pe)

        self.obs.set_doassim(it)
        if all([obs.doassim == 1 for obs in self.obs.obs_list]):
            distribute_state = distribtor.distribute_ocean_and_atmosphere_state
        else:
            distribute_state = distribtor.distribute_atmosphere_state

        status:int = 0
        if self.options.filtertype == 100:
            status = \
                PDAF.omi_generate_obs(colltor.collect_state_pdaf,
                                     distribtor.distribute_ocean_and_atmosphere_state,
                                     self.obs.init_dim_obs_gen_pdafomi,
                                     self.obs.obs_op_pdafomi,
                                     self.obs.get_obs_f,
                                     ensemble_processor.prepostprocess,
                                     distribtor.next_observation_pdaf)
        elif self.local.local_filter:
            status = \
                PDAF.localomi_assimilate(colltor.collect_state_pdaf,
                                        distribute_state,
                                        self.obs.init_dim_obs_pdafomi,
                                        self.obs.obs_op_pdafomi,
                                        ensemble_processor.prepostprocess,
                                        self.local.init_n_domains_pdaf,
                                        self.local.init_dim_l_pdaf,
                                        self.obs.init_dim_obs_l_pdafomi,
                                        distribtor.next_observation_pdaf,
                                        status)
        else:
            status = \
                PDAF.omi_assimilate_global(colltor.collect_state_pdaf,
                                          distribute_state,
                                          self.obs.init_dim_obs_pdafomi,
                                          self.obs.obs_op_pdafomi,
                                          ensemble_processor.prepostprocess,
                                          distribtor.next_observation_pdaf)

        assert status == 0, f'ERROR {status} in PDAF_put_state - stopping!'

    def forward(self):
        """time forward DA system

        Parameters
        ----------
        isStrong : bool
            Whether strong coupling is used
        """
        t = self.model.t0
        for step in range(self.model.total_steps):
            # output for analysis
            if (step % self.model.tw) < self.model.dt:
                if self.pe.mype_ens == 0: log.logger.info(f'a {step}')
                self.model.writer.write(t, 'a', self.model.field_p)
            #model  forward
            t = self.model.step()
            # output for forecast
            if ((step + 1) % self.model.tw) < self.model.dt:
                if self.pe.mype_ens == 0: log.logger.info(f'f {step}')
                self.model.writer.write(t, 'f', self.model.field_p)

            self.assimilate_pdaf(step+1)


    def finalise(self):
        PDAF.print_info(11)
        if self.pe.mype_ens==0: PDAF.print_info(3)

        PDAF.deallocate()
        self.model.finalise()
