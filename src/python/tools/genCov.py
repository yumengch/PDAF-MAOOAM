import sys
sys.path.append('/home/yumengch/NCEO/MAOOAM/qgs-0.2.5/')

from qgs.params.params import QgParams

import datetime
import numpy as np
import xarray as xr
import scipy

import pyPDAF.PDAF as PDAF
import functools
import time

def getVarAttrs():
    fieldnames = ['psi_a', 'T_a', 'psi_o', 'T_o']
    field_standard_name = ['atmosphere_streamfunction',
                           'atmosphere_temperature',
                           'ocean_streamfunction',
                           'ocean_temperature',]
    field_long_name = ['streamfunction in the atmosphere',
                       'temperature in the atmosphere',
                       'streamfunction in the ocean',
                       'temperature in the ocean',
                      ]

    output_varnames = ['sigma', ] \
                    + [name + '_mean' for name in fieldnames] \
                    + [name + '_svd' for name in fieldnames]
    standard_name = ['singular_value', ] \
                  + ['mean_' + name for name in field_standard_name] \
                  + ['singular_vector_' + name for name in field_standard_name]
    long_name = ['singular value of the state vector', ] \
              + ['temporal mean of the ' + name for name in field_long_name] \
              + ['singular vector of the ' + name for name in field_long_name]

    dims = [('rank',), 
            ('ny', 'nx' ), ('ny', 'nx' ), ('ny', 'nx' ), ('ny', 'nx' ),
            ('rank', 'ny', 'nx', ), ('rank', 'ny', 'nx', ), 
            ('rank', 'ny', 'nx', ), ('rank', 'ny', 'nx', )]
    return output_varnames, standard_name, long_name, dims


def getFileAttrs():
    attrs = dict()
    attrs['Conventions'] = 'CF-1.8'
    attrs['title'] = 'NetCDF output of initial ensemble covariance from MAOOAM-pyPDAF'
    attrs['institution'] = 'NCEO-AWI-UoR'
    attrs['source'] = 'pyPDAF and MAOOAM'
    attrs['history'] = f'{datetime.datetime.now().isoformat(timespec="seconds")}: Data created'
    attrs['reference'] = 'https://github.com/yumengch/pyPDAF'
    return attrs


def init_model_params(n):
    # Setting some model parameters
    # Model parameters instantiation with default specs
    model_parameters = QgParams({'n' : n})
    # Mode truncation at the wavenumber 2 in both x and y spatial coordinate
    model_parameters.set_atmospheric_channel_fourier_modes(6, 6)
    # Mode truncation at the wavenumber 2 in the x and at the
    # wavenumber 4 in the y spatial coordinates for the ocean
    model_parameters.set_oceanic_basin_fourier_modes(2, 4)

    # Setting MAOOAM parameters according to the publication linked above
    model_parameters.set_params({'kd': 0.0290, 'kdp': 0.0290, 'n': n, 'r': 1.e-7,
                                 'h': 136.5, 'd': 1.1e-7})
    model_parameters.atemperature_params.set_params({'eps': 0.7, 'T0': 289.3, 'hlambda': 15.06, })
    model_parameters.gotemperature_params.set_params({'gamma': 5.6e8, 'T0': 301.46})

    model_parameters.atemperature_params.set_insolation(103.3333, 0)
    model_parameters.gotemperature_params.set_insolation(310., 0)

    return model_parameters


def toPhysical(model_parameters, coeffs, xc, yc):
    # define spatial domain
    natm, noc = model_parameters.nmod
    nt = len(coeffs)

    psi_a = np.zeros((nt, ) + xc.shape)
    T_a = np.zeros_like(psi_a)
    psi_o = np.zeros_like(psi_a)
    T_o = np.zeros_like(psi_a)

    # get atmospheric components
    basis = model_parameters.atmospheric_basis.num_functions()
    for i, b in enumerate(basis):
        psi_a += np.outer(coeffs[:, i], b(xc, yc)).reshape((nt, ) + xc.shape)
        T_a += np.outer(coeffs[:, i+natm], b(xc, yc)).reshape((nt, ) + xc.shape)

    # get ocean components
    basis = model_parameters.oceanic_basis.num_functions()
    for i, b in enumerate(basis):
        psi_o += np.outer(coeffs[:, i+2*natm], b(xc, yc)).reshape((nt, ) + xc.shape)
        T_o += np.outer(coeffs[:, i+2*natm+noc], b(xc, yc)).reshape((nt, ) + xc.shape)

    return psi_a, T_a, psi_o, T_o


def distributeData(svals, svdU, nx, ny, nrank, meanstate):
    size = nx*ny
    data = [svals, ] \
         + [meanstate[i*size:(i+1)*size].reshape(ny, nx) for i in range(4)] \
         + [svdU[:, i*size:(i+1)*size].reshape(nrank, ny, nx) for i in range(4)]
    return data


def getDataArrays(svals, svdU, nx, ny, nrank, meanstate):
    da = dict()
    attrs = dict()

    data = distributeData(svals, svdU, nx, ny, nrank, meanstate)
    VarAttrs = zip(*getVarAttrs())
    for i, attr in enumerate(VarAttrs):
        varname, standard_name, long_name, dim = attr
        attrs['standard_name'] = standard_name
        attrs['long_name'] = long_name
        da[varname] = xr.DataArray(data[i],
                                     dims=dim, name=varname, attrs=attrs)
    return da


if __name__ == '__main__':
    # initialise model configurations
    n = 1.5
    mp = init_model_params(n)
    # read model trajectory
    ds = xr.open_dataset('trajectory.nc')
    nt = ds.dims['time']
    fields = np.hstack([ds.psi_a_f, ds.T_a_f, ds.psi_o_f, ds.T_o_f])[..., 0]
    # convert model from Fourier to physical space
    x0 = 0; x1 = 2*np.pi/n
    y0 = 0; y1 = np.pi
    nx = 9
    ny = 17
    X = np.linspace(x0, x1, nx)
    Y = np.linspace(y0, y1, ny)
    xc, yc = np.meshgrid(X, Y)
    psi_a, T_a, psi_o, T_o = toPhysical(mp, fields, xc, yc)

    # do EOF decomposition using PDAF
    maxtimes = nt
    nrank = nt
    size = nx*ny
    dim_state = np.array([size, size, size, size])
    offsets = np.array([0, size, 2*size, 3*size])
    state = np.hstack([psi_a.reshape(nt, size), T_a.reshape(nt, size), psi_o.reshape(nt, size), T_o.reshape(nt, size)])
    state = state.T
    print (state.shape)
    _, stddev, svals, svdU, meanstate, status = PDAF.eofcovar(dim_state, offsets, 1, 1,
                                                              state, state.mean(axis=1), 3)
    assert status == 0, 'Error in PDAF.EOFcovar'
    # save the covariance matrix
    svdU = svdU.T
    ds = xr.Dataset(getDataArrays(svals, svdU, nx, ny, nrank, meanstate),
                    attrs=getFileAttrs())
    ds.to_netcdf('covariance.nc')
