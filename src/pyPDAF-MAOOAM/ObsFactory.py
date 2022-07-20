from Obs import Obs
from ObsWriter import ObsWriter


class ObsFactory:
    def __init__(self, config, mype_filter, model):
        # create new observations
        self.obs = dict()
        for i, (key, item) in enumerate(config.items(), 1):
            if key != 'n_obs':
                self.obs[item] = Obs(i, item, mype_filter, model.nx)

    def setWriter(self, pe, model):
        if pe.filterpe:
            self.writer = dict()
            for obsname in self.obs:
                self.writer[obsname] = ObsWriter(f'MAOOAM_{obsname}.nc', model)

    def __getitem__(self, key):
        return self.obs[key]

    def items(self):
        return self.obs.items()

    def values(self):
        return self.obs.values()

    @property
    def nobs(self):
        return len(self.obs)

    @property
    def delt_obs(self):
        return min([obs.delt_obs for obs in self.obs.values()])