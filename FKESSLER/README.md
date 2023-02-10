# FKESSLER diagnostic script
A diagnostic script to analyze the FKESSLER moist baroclinic wave test case with KESSLER microphysics and is specifically designed to validate experimental data produced by the CAM-MPAS model. The analysis script looks for the surface pressure (PS) and precipitation rate variables (PRECL) in day 10 of the history file. This is prescribed by: https://www2.cesm.ucar.edu/models/simpler-models/fkessler/index.html

See that URL for background on this test.
This script is an attempt to replace fkessler.ncl script provided on that site with one based on Python and Xarrays.

The script presupposes that the NCAR geocat environment is available with Conda. See https://geocat.ucar.edu for information about using the geocat environment.

To install the geocat environment, follow these directions:

>conda create -n geocat -c conda-forge -c ncar geocat-comp geocat-datafiles

>conda activate geocat

>conda install -c conda-forge matplotlib cartopy jupyter

>conda install netcdf4

Running on your laptop would follow the pattern:

>conda activate geocat

>python fkessler_diag.py PATH_TO_FILES (optional arguments)

Typing

>fkessler_diag.py -h 

will give a list of the optional arguments:

optional arguments:
  -h, --help  show this help message and exit
  --saveplt   save the plot to your home directory (default: display plot)
  --gblplt    plot global field (default: plot prescribed region - (20:70N, 150W:30E)
  --verbose   enable verbose output

Leaving the conda environment can be accomplished by:

> conda deactivate

The basic run pattern in an HPC environment (for example on casper) is:

> module load conda

> conda activate geocat

> python fkessler_diag.py PATH_TO_FILES (optional arguments)

> conda deactivate

KNOWN ISSUES:

1) A sharp eyed user will notice that the Longitude location of the baroclinic disturbances are located 180 degrees out of phase with those shown on the CESM website. To date, we have been unable to identify the origin of this difference. If anyone figures it out let us know.


