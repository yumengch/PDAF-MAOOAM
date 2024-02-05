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

import mod_model

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
    model_parameters = QgParams({'n' : n, 'kd': 0.02, 'kdp': 0.02, 'n': n, 'r': 1.e-7,
                                 'h': 165, 'd': 9e-8})
    # Mode truncation at the wavenumber 2 in both x and y spatial coordinate
    model_parameters.set_atmospheric_channel_fourier_modes(2, 2)
    # Mode truncation at the wavenumber 2 in the x and at the
    # wavenumber 4 in the y spatial coordinates for the ocean
    model_parameters.set_oceanic_basin_fourier_modes(2, 4)

    # Setting MAOOAM parameters according to the publication linked above
    model_parameters.set_params({'kd': 0.02, 'kdp': 0.02, 'n': n, 'r': 1.e-7,
                                      'h': 165, 'd': 9e-8})
    model_parameters.atemperature_params.set_params({'eps': 0.7, 'T0': 290.2, 'hlambda': 15.06, })
    model_parameters.gotemperature_params.set_params({'gamma': 6.6e8, 'T0': 299.35})
    model_parameters.atemperature_params.set_insolation(103.3333, 0)
    model_parameters.gotemperature_params.set_insolation(310, 0)

    return model_parameters


def toPhysicalPy(model_parameters, coeffs, xc, yc):
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

    return psi_a[0], T_a[0], psi_o[0], T_o[0]

def toPhysical(coeffs, psi_a, T_a, psi_o, T_o):
    mod_model.mod_model.tophysical_a_multistep(coeffs, psi_a, T_a)
    mod_model.mod_model.tophysical_o_multistep(coeffs, psi_o, T_o)
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
    import time
    # initialise model configurations
    n = 1.5
    mp = init_model_params(n)
    # read model trajectory
    ds = xr.open_dataset('trajectory.nc', decode_times=False)
    nt = ds.dims['time']
    fields = np.hstack([ds.psi_a_f, ds.T_a_f, ds.psi_o_f, ds.T_o_f])
    fields = np.asfortranarray(fields)

    start = time.perf_counter()
    # convert model from Fourier to physical space
    nx = 129 
    ny = 129
    # avoiding excessive memory use
    psi_a = np.zeros((nt, nx, ny), order='F')
    T_a = np.zeros_like(psi_a, order='F')
    psi_o = np.zeros_like(psi_a, order='F')
    T_o = np.zeros_like(psi_a, order='F')   

    mod_model.mod_model.initialize_model()
    psi_a, T_a, psi_o, T_o = toPhysical(fields, psi_a, T_a, psi_o, T_o)
    psi_a = np.ascontiguousarray(psi_a)
    T_a = np.ascontiguousarray(T_a)
    psi_o = np.ascontiguousarray(psi_o)
    T_o = np.ascontiguousarray(T_o)

    # do EOF decomposition using PDAF
    start = time.perf_counter()
    maxtimes = nt
    nrank = nt
    size = nx*ny
    dim_state = np.array([size, size, size, size])
    offsets = np.array([0, size, 2*size, 3*size])
    state = np.hstack([psi_a.reshape(nt, size), T_a.reshape(nt, size), psi_o.reshape(nt, size), T_o.reshape(nt, size)])
    state = state.T
    _, stddev, svals, svdU, meanstate, status = PDAF.eofcovar(dim_state, offsets, 1, 1,
                                                              state, state.mean(axis=1), 3)
    assert status == 0, 'Error in PDAF.EOFcovar'
    print ('genCovar time', time.perf_counter() - start)

    start = time.perf_counter()
    # save the covariance matrix
    ds = xr.Dataset(getDataArrays(svals, svdU.T, nx, ny, nrank, meanstate),
                     attrs=getFileAttrs())
    ds.to_netcdf('covariance.nc')
    print ('write time', time.perf_counter() - start)
