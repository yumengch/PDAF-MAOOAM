
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
import U_PDAF
import PDAF_caller

from parallelization import parallelization
from StateVector import StateVector
from FilterOptions import FilterOptions
from Inflation import Inflation
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
    model : `Model.Model`
        model object
    obs : `ObsFactory.ObsFactory`
        a factory of observations

    screen : int
        verbosity of PDAF screen output
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
        # init parallelisation
        self.pe = parallelization(config['Ensemble'], self.screen)
        # init model
        self.model = Model(config['Model'], self.pe.task_id)
        # init options
        self.sv = StateVector(self.model, self.pe.n_modeltasks)
        self.options = FilterOptions(config['FilterOptions'])
        self.infl = Inflation(config['Inflation'])

        # Initial Screen output
        if (self.pe.mype_world == 0):
            print('+++++ pyPDAF online mode +++++')
            print('+++++ MAOOAM +++++')

        self.UserFuncs = U_PDAF.PDAFUserFuncs(self)
        PDAF_caller.init_pdaf(self)

    def forward(self):
        """time forward DA system
        """
        t = self.model.t0
        for step in range(self.model.total_steps):
            # output for analysis
            if (step % self.model.tw) < self.model.dt:
                if self.pe.mype_world == 0: print(('a', step))
                self.model.writer.write(t, 'a', self.model.field_p)
            # model forward
            t = self.model.step()
            # output for forecast
            if ((step + 1) % self.model.tw) < self.model.dt:
                if self.pe.mype_world == 0: print(('f', step))
                self.model.writer.write(t, 'f', self.model.field_p)


    def finalise(self):
        import pyPDAF.PDAF as PDAF
        PDAF.print_info(2)
        PDAF.print_info(11)
        PDAF.print_info(3)

        PDAF.deallocate()

        self.pe.finalise()
