import datetime
import numpy as np
import xarray as xr
import netCDF4

# from MAOOAM.model.params_maooam import natm, noc


class ObsWriter:
    def __init__(self, filename, model):
        self.nc = netCDF4.Dataset(filename,'w')

        attrs = self.getFileAttrs()
        for attr, value in attrs.items():
            setattr(self.nc, attr, value)

        self.nc.createDimension('time', None)
        time=self.nc.createVariable('time', np.dtype('float64').char, 'time')
        time.long_name = 'time'
        time.units = "days since 1900-1-1 0:0:0"

        self.natm, self.noc = model.model_parameters.nmod
        for dim, dimlen in zip(['natm', 'noc'], [self.natm, self.noc]):
            self.nc.createDimension(dim, dimlen)

        for attr in zip(*self.getVarAttrs()):
            varname, standard_name, long_name, dim = attr
            var = self.nc.createVariable(varname, np.dtype('float64').char, dim)
            var.standard_name = standard_name
            var.long_name = long_name

        self.time_count = 0


    def write(self, step, inputData):
        data = self.distributeData(inputData)
        self.nc['time'][self.time_count] = step
        for i, varname in enumerate(['psi_a', 'T_a', 
                                     'psi_o', 'T_o']):
            self.nc[varname][self.time_count] = data[i]
        self.time_count += 1


    def getVarAttrs(self):
        fieldnames = ['psi_a', 'T_a', 'psi_o', 'T_o']
        field_standard_name = ['atmosphere_streamfunction_coefficient',
                               'atmosphere_temperature_coefficient',
                               'ocean_streamfunction_coefficient',
                               'ocean_temperature_coefficient']
        field_long_name = ['coefficient of streamfunction in the atmosphere',
                           'coefficient of temperature in the atmosphere',
                           'coefficient of streamfunction in the ocean',
                           'coefficient of temperature in the ocean'
                          ]

        dims = [('time', 'natm'), ('time', 'natm'), ('time', 'noc'), ('time', 'noc')]
        return fieldnames, field_standard_name, field_long_name, dims


    def getFileAttrs(self):
        attrs = dict()
        attrs['Conventions'] = 'CF-1.8'
        attrs['title'] = 'NetCDF output of synthetic observations from MAOOAM-pyPDAF'
        attrs['institution'] = 'NCEO-AWI-UoR'
        attrs['source'] = 'pyPDAF and MAOOAM'
        attrs['history'] = f'{datetime.datetime.now().isoformat(timespec="seconds")}: Data created'
        attrs['reference'] = 'https://github.com/yumengch/pyPDAF'
        return attrs


    def distributeData(self, inputData):
        offsets = np.array([0, self.natm, 2*self.natm, 2*self.natm + self.noc, 2*(self.natm + self.noc)])
        output = [inputData[offsets[i]:offsets[i+1]] for i in range(4)]
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

    def __del__(self):
        self.nc.close()