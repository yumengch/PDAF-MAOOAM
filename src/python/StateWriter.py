import datetime
import numpy as np
import xarray as xr
import netCDF4

# from MAOOAM.model.params_maooam import natm, noc


class StateWriter:
    def __init__(self, filename, dim_ens, model):
        self.nc = netCDF4.Dataset(filename,'w')

        attrs = self.getFileAttrs()
        for attr, value in attrs.items():
            setattr(self.nc, attr, value)

        self.nc.createDimension('time', None)
        self.nx, self.ny = model.nx, model.ny
        time=self.nc.createVariable('time', np.dtype('float64').char, 'time')
        time.long_name = 'time'
        time.units = "days since 1900-1-1 0:0:0"
        for dim, dimlen in zip(['ens', 'ny', 'nx'], [dim_ens, self.ny, self.nx]):
            self.nc.createDimension(dim, dimlen)

        for vartype, typename in zip(['f', 'a'], ['forecast', 'analysis']):
            for attr in zip(*self.getVarAttrs()):
                varname, standard_name, long_name, dim = attr
                var = self.nc.createVariable(f'{varname}_{vartype}', np.dtype('float64').char, dim)
                var.standard_name = standard_name + '_' + typename
                var.long_name = typename + ' of ' + long_name

        self.time_count = 0


    def write(self, step, typename, inputData):
        data = self.distributeData(inputData)
        self.nc['time'][self.time_count] = step
        for i, varname in enumerate([f'psi_a_{typename}', f'T_a_{typename}', 
                                     f'psi_o_{typename}', f'T_o_{typename}']):
            self.nc[varname][self.time_count] = data[i]
        if int(self.nc['time'][-1]) != step:
            self.time_count += 1


    def getVarAttrs(self):
        fieldnames = ['psi_a', 'T_a', 'psi_o', 'T_o']
        field_standard_name = ['atmosphere_streamfunction',
                               'atmosphere_temperature',
                               'ocean_streamfunction',
                               'ocean_temperature']
        field_long_name = ['streamfunction in the atmosphere',
                           'temperature in the atmosphere',
                           'streamfunction in the ocean',
                           'temperature in the ocean'
                          ]

        dims = [('time', 'ny', 'nx', 'ens'), ('time', 'ny', 'nx', 'ens'), 
                ('time', 'ny', 'nx', 'ens'), ('time', 'ny', 'nx', 'ens')]
        return fieldnames, field_standard_name, field_long_name, dims


    def getFileAttrs(self):
        attrs = dict()
        attrs['Conventions'] = 'CF-1.8'
        attrs['title'] = 'NetCDF output from MAOOAM-pyPDAF'
        attrs['institution'] = 'NCEO-AWI-UoR'
        attrs['source'] = 'MAOOAM-pyPDAF'
        attrs['history'] = f'{datetime.datetime.now().isoformat(timespec="seconds")}: Data created'
        attrs['reference'] = 'https://github.com/yumengch/pyPDAF'
        return attrs


    def distributeData(self, inputData):
        size = self.nx*self.ny
        output = [inputData[i*size:(i+1)*size].reshape(self.ny, self.nx) for i in range(4)]
        return output


    def getDataArrays(self, inputData):
        da = dict()
        attrs = dict()

        data = self.distributeData(inputData)
        VarAttrs = zip(*self.getVarAttrs())
        for i, attr in enumerate(VarAttrs):
            varname, standard_name, long_name, dim = attr
            attrs['standard_name'] = standard_name
            attrs['long_name'] = long_name
            da[varname] = xr.DataArray(data[i],
                                         dims=dim, name=varname, attrs=attrs)
        return da