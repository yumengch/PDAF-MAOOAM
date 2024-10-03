# A suite of experiments using PDAF and pyPDAF with MAOOAM

## Fortran system
In the Fortran PDAF DA system, running the experiment requires compile the code with Makefile subjected to the path to required libraries. The experiment also requires a link or copy of the MAOOAM model.
The experiments are:
- [fortran-alpha](src/fortran-alpha): a strongly coupled DA system where the amount of assimilated observation is adjusted by an alpha coefficient
- [fortran-deflate](src/fortran-deflate): a strongly coupled DA where only the ensemble spread of observed component is inflated.
- [fortranFree](src/fortranFree): Free run ensemble
- [fortranGenCov](src/fortranGenCov): generation of error covariance matrix from a long model trajectory.
- [fortranStrong](src/fortranStrong): A standard strongly coupled DA system with PDAF
- [fortranStrongSingle](src/fortranStrongSingle): A standard strongly coupled DA system where only one observation is assimilated
- [fortranWeak](src/fortranWeak): A weakly coupled DA system with PDAF

## Python system
In the pyPDAF DA systems, running the experiment requires the Python library in [model](src/model) compiled by Makefile subjected to the path to required libraries.
The Python model library read the namelist of the original MAOOAM model as well.
The experiments are:
- [pythonFree](src/pythonFree): pyPDAF ensemble free run system
- [pythonGenCov](src/pythonGenCov): pyPDAF calls to generating error covariance matrix from model trajectories 
- [pythonStrong](src/pythonStrong): Strongly coupled DA system in Python
- [pythonStrongPDAFlocal](src/pythonStrongPDAFlocal): Strongly coupled DA system using `PDAFlocal` module
- [pythonStrongSingle](src/pythonStrongSingle): Strongly coupeld DA with only one observaiton being assimilated
- [pythonWeak](src/pythonWeak): Weakly coupled DA system
