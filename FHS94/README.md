# FHS94 diagnostic script
A diagnostic script to analyze the HELD-SUAREZ test case which replaces the CAM physics package with a simple relaxation of the temperature field toward a zonally symmetric equilibrium profile and simple linear drag at the lower boundary. This test case follows the specifications outlined by Held and Suarez (1994), Bull. Amer. Met. Soc., 75, 1825-1830. 

The analysis script looks for the U,V and T density (rho): it computes pressure from the density (rho) field using the ideal gas law. The test case and compset is described in more detail at: https://www2.cesm.ucar.edu/models/simpler-models/held-suarez.html

This script is an attempt to replace FHS94_test.ncl script provided on that site with one based on Python and Xarrays.

The script presupposes that the NCAR geocat environment is available with Conda. See https://geocat.ucar.edu for information about using the geocat environment.

To install the geocat environment, follow these directions:

>conda create -n geocat -c conda-forge -c ncar geocat-comp geocat-datafiles

>conda activate geocat

>conda install -c conda-forge matplotlib cartopy jupyter

>conda install netcdf4

Running on your laptop would follow the pattern:

>conda activate geocat

>python fhs94_diag.py PATH_TO_FILES (optional arguments)

Typing

>fhs94_diag.py -h 

will give a list of the optional arguments:

optional arguments:
  -h, --help  show this help message and exit
  --saveplt   save the plot to a file (default: display plot)
  --verbose   enable verbose output

The plot will be saved to your home directory.

Leaving the conda environment can be accomplished by:

> conda deactivate

The basic run pattern in an HPC environment (for example on casper) is:

> module load conda

> conda activate geocat

> python fhs94_diag.py PATH_TO_FILES (optional arguments)

> conda deactivate

KNOWN ISSUES:

1. The script was developed using test files without a proper time axis. Therefore the analyis time window is determined crudely, using a hardwired data file name pattern, namely:

"latlon_val.FHS94.mpasa120.che.gnu.cam.h1.0002-01*"

For a full 1200 day run, this results in the window starting after 365 days (rather than day 200), with only 835 days (rather than the prescribed 1000 days) being analyzed. The results are still useful for validation of the model.

2. There is a know minor plotting bug: "interpolation point out of data bounds" that occurs at the bottom of the atmosphere:

/Users/loft/.conda/envs/geocat/lib/python3.9/site-packages/metpy/interpolate/one_dimension.py:137: UserWarning: Interpolation point out of data bounds encountered
  warnings.warn('Interpolation point out of data bounds encountered')
  
3. There is an option in the script called interp_early that should be left as False, as it will result in painfully long script run times. As it is, with interp_early=False set, a full analysis run of an FHS94 dataset will still take 30 minutes on a single node NCAR's Casper at 120 km/ 32 levels resolution. We plan to address performance of the script using Dask parallelism in a future release. 

