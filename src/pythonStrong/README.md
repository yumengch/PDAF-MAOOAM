# MAOOAM-pyPDAF

This is an online Data Assimilation System (DAS) for [the Modular Arbitrary-Order Ocean-Atmosphere Model (MAOOAM)](https://github.com/Climdyn/MAOOAM.git) using pyPDAF. The MAOOAM is solved in a Fourier space and the state vector is the coefficients of Fourier mode. The twin experiments observe and assimilate in the Fourier modes instead of the physical space.

## Generating trajectory
In MAOOAM directory, run:
```bash
python tools/genTrajectory.py
```

## Generating initial covariance matrix for the ensemble
After generating model trajectory, in MAOOAM directory, run:
```bash
python tools/genCov.py
```

## Setting PDAF and observation options
The options for PDAF can be set at [configuration.ini](configuration.ini). In particular,
```
[Obs]
obs1 = ObsA
obs2 = ObsB
...
```
The `ObsA` is the name of the corresponding observation configuration files, which corresponds to [ObsA.ini](ObsA.ini). The observation filename is given by default as `MAOOAM_{obsname}.nc`, such as `MAOOAM_ObsA.nc`. It can also be given by `ObsA.ini`.

### Generating synthetic observations
Setting the `filtertype = 100`.

## Running the DA
`python main.py`
