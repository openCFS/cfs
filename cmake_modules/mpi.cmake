# here we concentrate the check for the proper MPI installation
# call it in the USE_MPI case!

# We need MPI for phist and petsc. 

# We don't build MPI by cfsdeps but require it to set CC, CXX and FC to mpic, mpicxx and mpif90 version!
# This we test here commonly.
#
# There can be several MPI installations (e.g. openmpi, mpich, intel)
# The current must be selected by setting the CC, CXX and FC environment variables!
#
# to make use of the mpi libs (petsc and phist), cfs needs to be called by mpirun. 
# If you have no scaling, try mpirun.hydra from mpic
if(NOT CMAKE_C_COMPILER MATCHES "mpi")
  message(FATAL_ERROR "For MPI usage CC needs to be set to a mpicc compiler. Recreate the build directory!")
endif() 
if(NOT CMAKE_CXX_COMPILER MATCHES "mpi")
  message(FATAL_ERROR "For MPI usage CXX needs to be set to a mpicxx compiler. Recreate the build directory!")
endif()
if(NOT CMAKE_Fortran_COMPILER MATCHES "mpi") # case sensitive! workfs only like this with cmake 3.9.4
 message(FATAL_ERROR "For MPI usage FC needs to be set to a mpif90 compiler. Recreate the build directory!")
endif()

# search for the mpi version in the path
# CC=/usr/lib64/mpi/gcc/openmpi/bin/mpicc
get_filename_component(MPI_BIN ${CMAKE_C_COMPILER} DIRECTORY) # is exported in def_use_mpi.hh.in
get_filename_component(MPI_BASE_DIR ${MPI_BIN} DIRECTORY)
get_filename_component(MPI_BASE ${MPI_BASE_DIR} NAME) # is exported in def_use_mpi.hh.in
# MPI_BASE shall be openmpi, mpich or intel64 or anything local installed with a version on it
# we don't check for common base for CC, CXX and FC but it is necessary for MPI!

