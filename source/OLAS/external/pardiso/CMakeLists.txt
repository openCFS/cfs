
IF(EXPLICIT_TEMPLATE_INSTANTIATION)
    SET(PARDISO_SRCS pardisosolver.cc)
ELSE(EXPLICIT_TEMPLATE_INSTANTIATION)
    SET(PARDISO_SRCS PardisoDummy.cc)
ENDIF(EXPLICIT_TEMPLATE_INSTANTIATION)

ADD_LIBRARY(pardiso-olas STATIC ${PARDISO_SRCS})

# CFS_FORTRAN_LIBS has been set in cmake_modules/distro.cmake
# LAPACK_LIBRARY and BLAS_LIBRARY are defined in
# cmake_modules/FindFortranLibs.cmake
TARGET_LINK_LIBRARIES(pardiso-olas
  cfsgeneral
  datainout
  graph-olas
  logging
  matvec
  paramh
  solver-olas
  pardiso-lib1
)
# needed as solver-olas is cyclic dependent on us 
SET_PROPERTY(TARGET pardiso-olas PROPERTY LINK_INTERFACE_MULTIPLICITY 3)

SET(DEPENDS
  "${LAPACK_LIBRARY}"
  "${BLAS_LIBRARY}"
  "${CFS_FORTRAN_LIBS}"
)

ADD_DEPENDENT_LIBRARY(pardiso-lib "${PARDISO_LIBRARY}" "${DEPENDS}")
