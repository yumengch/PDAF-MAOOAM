import configparser


class PDAFConfig:
    """This class controls the PDAF-MAOOAM configuration/namelist reading
    """
    _initialised = False
    def __init__(self, configfile='configuration.ini'):
        assert not PDAFConfig._initialised, 'PDAF configuration is already set up!!!'
        PDAFConfig._initialised = True

        self.config = configparser.ConfigParser()
        self.config.read(configfile)

    def __getitem__(self, key):
        return self.config[key]