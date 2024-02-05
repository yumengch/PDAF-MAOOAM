from Obs import Obs


class ObsFactory:
    def __init__(self, config, mype, model):
        # create new observations
        self.obs = dict()
        for i, (key, item) in enumerate(config.items(), 1):
            if key != 'n_obs':
                self.obs[key] = Obs(i, item, mype, model.nx)

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
