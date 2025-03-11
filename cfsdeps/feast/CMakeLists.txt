cmake_minimum_required(VERSION 3.0)

project(FEAST Fortran)
set(FEAST_VERSION 4.0)

# we remove all p* source and header files as we don't support mpi
set(FEAST_SRCS
  ${FEAST_VERSION}/src/kernel/libnum.f90
  ${FEAST_VERSION}/src/kernel/feast_aux.f90
  ${FEAST_VERSION}/src/kernel/feast_tools.f90
  ${FEAST_VERSION}/src/kernel/dzfeast.f90
  ${FEAST_VERSION}/src/sparse/dzfeast_pev_sparse.f90
  ${FEAST_VERSION}/src/sparse/dzifeast_pev_sparse.f90
  ${FEAST_VERSION}/src/sparse/dzlsprim.f90 
  ${FEAST_VERSION}/src/sparse/dzfeast_sparse.f90
  ${FEAST_VERSION}/src/sparse/dzifeast_sparse.f90
  ${FEAST_VERSION}/src/sparse/sclsprim.f90
  ${FEAST_VERSION}/src/dense/dzfeast_dense.f90
  ${FEAST_VERSION}/src/dense/dzfeast_pev_dense.f90
  ${FEAST_VERSION}/src/banded/dzfeast_banded.f90 )

# we have our own header with guard and name mangling maintained in cfsdeps/feast
set(FEAST_HEADERS ${FEAST_VERSION}/include/feast_cpp.h )

include_directories(${FEAST_VERSION}/src/kernel)

add_library(feast STATIC ${FEAST_SRCS})

# set special CMAKE_Fortran_FLAGS flags like -cpp and for gfortran -ffree-line-length-none -ffixed-line-length-none -fallow-argument-mismatch

install(TARGETS feast ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${FEAST_HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
