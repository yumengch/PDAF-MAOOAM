import os
import xarray as xr


def mergeData(dirname):
    E = [xr.open_dataset('{:}/maooam_{:03}.nc'.format(dirname, i), decode_times=False) for i in range(1, 17)]
    E = [f.assign_coords({'ens': (['ens',], [i],
                                  {'long_name': 'ensemble index',
                                   'standard_name': 'ensemble number',
                                   'units': '1'
                                  }
                                 )
                         }) for i, f in enumerate(E)]
    return xr.combine_by_coords(E, combine_attrs='drop_conflicts')


if __name__ == '__main__':
    # convert to single file
    mergeData('maooam_output').to_netcdf('maooam_A1O7_For.nc')
