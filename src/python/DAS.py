
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
import pyPDAF.PDAF.PDAFomi as PDAFomi
import U_PDAF
import U_PDAFomi
import PDAF_caller

from parallelization import parallelization
from ObsFactory import ObsFactory
from StateVector import StateVector
from FilterOptions import FilterOptions
from Inflation import Inflation
from Localization import Localization
from Model import Model
from ModelWriter import ModelWriter


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

    def init(self, config):
        """initialise DA system

        Parameters
        ----------
        config : `Config.PDAFConfig`
            configuration object
        """
        # init DAS options
        self.screen = config['Global'].getint('screen', 3)
        self.isStrong = config['Global'].getboolean('isStrong', True)
        self.is_freerun = config['Global'].getboolean('is_freerun', False)
        # init parallelisation
        self.pe = parallelization(config['Ensemble'], self.screen)
        # init model
        self.model = Model(config['Model'], self.pe.task_id)
        # init obs
        self.obs = ObsFactory(config['Obs'], self.pe.mype_world, self.model)
        # init options
        self.sv = StateVector(config['StateVector'], self.model, self.pe.n_modeltasks, self.isStrong)
        self.options = FilterOptions(config['FilterOptions'])
        self.infl = Inflation(config['Inflation'])
        # init obs. output
        if self.options.filtertype == 100:
            self.obs.setWriter(self.pe, self.model)

        # Initial Screen output
        if (self.pe.mype_world == 0):
            print('+++++ pyPDAF online mode +++++')
            print('+++++ MAOOAM +++++')

        if self.pe.filterpe:
            # init observations
            PDAFomi.init(self.obs.nobs)

        self.UserFuncs = U_PDAF.PDAFUserFuncs(self)
        self.UserFuncsO = U_PDAFomi.PDAFomiUserFuncs(self)
        PDAF_caller.init_pdaf(self)
        self.UserFuncs.firsttime = False

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
                if self.pe.mype_world == 0: print(('a', step))
                self.model.writer.write(t, 'a', self.model.field_p)
            #model  forward
            t = self.model.step(step)
            # output for forecast
            if ((step + 1) % self.model.tw) < self.model.dt:
                if self.pe.mype_world == 0: print(('f', step))
                self.model.writer.write(t, 'f', self.model.field_p)

            if ((not self.isStrong) and (self.sv.component == 'ao')):
                self.sv.setFields('a')
                self.obs['ObsA'].doassim = 1
                self.obs['ObsO'].doassim = 0
            PDAF_caller.assimilate_pdaf(self)

            if ((not self.isStrong) and (self.sv.component == 'ao')):
                self.sv.setFields('o')
                self.obs['ObsA'].doassim = 0
                self.obs['ObsO'].doassim = 1
                PDAF_caller.assimilate_pdaf(self)


    def __del__(self):
        import pyPDAF.PDAF as PDAF
        if (self.pe.mype_world==0): PDAF.print_info(2)
        if (self.pe.mype_world==0): PDAF.print_info(11)
        if (self.pe.mype_world==0): PDAF.print_info(3)

        PDAF.deallocate()
