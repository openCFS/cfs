CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

PROJECT(FEAST Fortran)
SET(FEAST_VERSION 3.0)

SET(FEAST_SRCS
  ${FEAST_VERSION}/src/kernel/f90_functions_wrapper.f90
  ${FEAST_VERSION}/src/kernel/feast_tools.f90
  ${FEAST_VERSION}/src/kernel/dzfeast.f90
  ${FEAST_VERSION}/src/kernel/scfeast.f90
  ${FEAST_VERSION}/src/sparse/dzfeast_sparse.f90
  ${FEAST_VERSION}/src/sparse/scfeast_sparse.f90
  ${FEAST_VERSION}/src/kernel/feast_aux.f90
)
SET(FEAST_HEADERS
  feast_cpp.h
)

SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /fpp")

INCLUDE_DIRECTORIES(
  ${FEAST_VERSION}/src/kernel
)

ADD_LIBRARY(feast STATIC ${FEAST_SRCS})

INSTALL(TARGETS feast ARCHIVE DESTINATION "${LIB_SUFFIX}")

INSTALL(FILES ${FEAST_HEADERS} DESTINATION include)
