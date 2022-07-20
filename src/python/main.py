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
import matplotlib.pyplot as plt
import numpy as np

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

def global_except_hook(exctype, value, traceback):
    sys.stderr.write("except_hook. Calling MPI_Abort().\n")
    # NOTE: mpi4py must be imported inside exception handler, not globally.
    # In chainermn, mpi4py import is carefully delayed, because
    # mpi4py automatically call MPI_Init() and cause a crash on Infiniband environment.
    mpi4py.MPI.COMM_WORLD.Abort(1)
    sys.__excepthook__(exctype, value, traceback)
# sys.excepthook = global_except_hook

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

    # for step in range(model.t0):
    #     model.step(pe, step, usePDAF)

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