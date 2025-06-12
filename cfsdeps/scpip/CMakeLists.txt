cmake_minimum_required(VERSION 3.16) # cmake 4.0 does not allow < 3.5

project(scpip Fortran)

set(SCPIP_SRCS Scpip30.f Sparsecz.f scpblas.f scpdummy.f)

set(SCPIP_LIBNAME scpip)

add_library(${SCPIP_LIBNAME} STATIC ${SCPIP_SRCS})

install(TARGETS ${SCPIP_LIBNAME} ARCHIVE DESTINATION "${LIB_SUFFIX}")

