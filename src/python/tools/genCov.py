from MAOOAM.model.params_maooam import natm, noc

import datetime
import numpy as np
import xarray as xr

import pyPDAF.PDAF as PDAF


def getVarAttrs():
    fieldnames = ['psi_a', 'T_a', 'psi_o', 'T_o']
    field_standard_name = ['atmosphere_streamfunction_coefficient',
                           'atmosphere_temperature_coefficient',
                           'ocean_streamfunction_coefficient',
                           'ocean_temperature_coefficient',]
    field_long_name = ['coefficient of streamfunction in the atmosphere',
                       'coefficient of temperature in the atmosphere',
                       'coefficient of streamfunction in the ocean',
                       'coefficient of temperature in the ocean',
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
            ('natm', ), ('natm', ), ('noc', ), ('noc', ),
            ('natm', 'rank'), ('natm', 'rank'), ('noc', 'rank'), ('noc', 'rank')]
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


def distributeData(svals, svdU, meanstate):
    offsets = np.array([0, natm, 2*natm, 2*natm + noc, 2*(natm + noc)])
    data = [svals, ] \
         + [meanstate[offsets[i]:offsets[i+1]] for i in range(4)] \
         + [svdU[offsets[i]:offsets[i+1]] for i in range(4)]
    return data


def getDataArrays(svals, svdU, meanstate):
    da = dict()
    attrs = dict()

    data = distributeData(svals, svdU, meanstate)
    VarAttrs = zip(*getVarAttrs())
    for i, attr in enumerate(VarAttrs):
        varname, standard_name, long_name, dim = attr
        attrs['standard_name'] = standard_name
        attrs['long_name'] = long_name
        da[varname] = xr.DataArray(data[i],
                                     dims=dim, name=varname, attrs=attrs)
    return da


if __name__ == '__main__':
    ds = xr.open_dataset('trajectory.nc')
    field = np.hstack([ds.psi_a, ds.T_a, ds.psi_o, ds.T_o])
    maxtimes, dim_state = field.shape
    dim_state = np.array([natm, natm, noc, noc])
    offsets = np.array([0, natm, 2*natm, 2*natm + noc])

    _, stddev, svals, svdU, meanstate, status = PDAF.eofcovar(dim_state, offsets, 1, 1,
                                                              field.T, field.mean(axis=0), 3)
    assert status == 0, 'Error in PDAF.EOFcovar'
    ds = xr.Dataset(getDataArrays(svals, svdU, meanstate),
                    attrs=getFileAttrs())
    ds.to_netcdf('covariance.nc')