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
from parallelization import parallelization
from ObsFactory import ObsFactory
from StateVector import StateVector
from FilterOptions import FilterOptions
from Inflation import Inflation
from Localization import Localization
from DAS import DAS
import Config
import Model

import sys
import time
import mpi4py.MPI

def main():
    config = Config.PDAFConfig()
    usePDAF = config['Global'].getboolean('usePDAF', True)
    screen = config['Global'].getint('screen', 3)

    if usePDAF:
        pe = parallelization(config['Ensemble'], screen=screen)

    # Initial Screen output
    if (pe.mype_world == 0):
        print('+++++ PDAF online mode +++++')
        print('+++++ MAOOAM +++++')

    # Initialize model
    model = Model.Model(config['Model'], pe=pe)

    obs = ObsFactory(config['Obs'], pe.mype_filter, model)

    das = DAS(pe, model, obs, screen=screen)
    sv = StateVector(das.model, dim_ens=das.pe.n_modeltasks)
    options = FilterOptions(config['FilterOptions'])
    infl = Inflation(config['Inflation'])
    local = Localization(config['Localization'])
    das.init(sv, options, infl, local)

    for it in range(model.total_steps):
        das.forward(it, usePDAF)

    pe.finalize_parallel()

if __name__ == '__main__':
    main()