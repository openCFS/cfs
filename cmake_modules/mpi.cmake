# here we concentrate the check for the proper MPI installation
# call it in the USE_MPI case!

# We need for petsc. 

# We don't build MPI by cfsdeps but require it to set CC, CXX and FC to mpic, mpicxx and mpif90 version!
# This we test here commonly.
#
# There can be several MPI installations (e.g. openmpi, mpich, intel)
# The current must be selected by setting the CC, CXX and FC environment variables!
#
# to make use of the mpi libs (petsc), cfs needs to be called by mpirun. 
# If you have no scaling, try mpirun.hydra from mpic

# Detect MPI using CMake's built‐in FindMPI module
include(FindMPI)
find_package(MPI REQUIRED COMPONENTS C CXX Fortran)
if(NOT MPI_C_FOUND OR NOT MPI_CXX_FOUND OR NOT MPI_Fortran_FOUND)
  message(FATAL_ERROR "MPI compiler wrappers (mpicc, mpicxx, mpif90) not found. Please install or load an MPI distribution.")
endif()

# search for the mpi version in the path
# CC=/usr/lib64/mpi/gcc/openmpi/bin/mpicc
get_filename_component(MPI_BIN ${MPI_C_COMPILER} DIRECTORY) # is exported in def_use_mpi.hh.in
get_filename_component(MPI_BASE_DIR ${MPI_BIN} DIRECTORY)
get_filename_component(MPI_BASE ${MPI_BASE_DIR} NAME) # is exported in def_use_mpi.hh.in
# MPI_BASE shall be openmpi, mpich or intel64 or anything local installed with a version on it
# we don't check for common base for CC, CXX and FC but it is necessary for MPI!

