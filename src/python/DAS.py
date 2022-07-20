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
import PDAF_caller
from ModelWriter import ModelWriter


class DAS:

    """Data assimilation system

    Attributes
    ----------
    sv : `StateVector.StateVector`
        an object of StateVector
    options : `FilterOptions.FilterOptions`
        filtering options
    infl : `Inflation.Inflation`
        inflation object
    local : `Localization.Localization`
        localization object
    model : `Model.Model`
        model object
    obs : `ObsFactory.ObsFactory`
        a factory of observations
    pe : `parallelization.parallelization`
        parallelization object
    screen : int
        verbosity of PDAF screen output
    """

    def __init__(self, pe, model, obs, screen):
        """constructor

        Parameters
        ----------
        pe : `parallelization.parallelization`
            parallelization object
        model : `Model.Model`
            model object
        obs : `Obs.Obs`
            observation object
        screen : int
            verbosity of PDAF screen output
        """
        self.pe = pe
        self.model = model
        self.obs = obs
        self.screen = screen

    def init(self, sv, options, infl, local):
        """initialise DA system

        Parameters
        ----------
        sv : `StateVector.StateVector`
            an object of StateVector
        options : `FilterOptions.FilterOptions`
            filtering options
        infl : `Inflation.Inflation`
            inflation object
        local : `Localization.Localization`
            a localization info obejct
        """
        # init model
        self.model.init_field()


        # Initialize PDAF
        self.sv = sv
        self.options = options
        self.infl = infl
        self.local = local

        if self.pe.filterpe:
            # init observations
            PDAFomi.init(self.obs.nobs)

        self.writer = None
        if self.pe.filterpe:
            self.writer = ModelWriter('MAOOAM.nc', self.sv.dim_ens, self.model)
        if self.options.filtertype == 100:
            self.obs.setWriter(self.pe, self.model)

        PDAF_caller.init_pdaf(self.sv, self.infl,
                              self.options,
                              self.local,
                              self.model, self.pe,
                              self.obs, self.writer, self.screen)

    def forward(self, step, usePDAF):
        """time forward DA system

        Parameters
        ----------
        step : int
            current time step
        usePDAF : bool
            whether PDAF is used
        """
        self.model.step(self.pe, step, usePDAF)
        if usePDAF:
            PDAF_caller.assimilate_pdaf(self.model, self.obs, self.pe,
                                        self.sv, self.local, self.writer,
                                        self.options.filtertype)
