# utilities

## checkout init.sh
.. to clean workspace before starting testcase

## simulation must be run in two separate terminal windows

### openCFS
- ´´´cd openfoam_opencfs_coupling/openCFS´´´
- ´´´cfs testcase´´´

### openFOAM
- ´´´cd openfoam_opencfs_coupling/openFoam´´´
- maybe you have to make openFOAM accessable beforehand ´´´source /usr/lib/openfoam/openfoam2412/etc/bashrc´´´
- ´´´icoFoam´´´

## at this point openFOAM/preCICE must be closed via CTRL-C