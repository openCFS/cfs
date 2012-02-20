
IF(EXPLICIT_TEMPLATE_INSTANTIATION)
    SET(PARDISO_SRCS pardisosolver.cc)
ELSE(EXPLICIT_TEMPLATE_INSTANTIATION)
    SET(PARDISO_SRCS PardisoDummy.cc)
ENDIF(EXPLICIT_TEMPLATE_INSTANTIATION)

ADD_LIBRARY(pardiso-olas STATIC ${PARDISO_SRCS})

SET(TARGET_LL
  cfsgeneral
  datainout
  graph-olas
  logging
  matvec
  paramh
  solver-olas
)
# needed as solver-olas is cyclic dependent on us 
SET_PROPERTY(TARGET pardiso-olas PROPERTY LINK_INTERFACE_MULTIPLICITY 3)

SET(DEPENDS
  "${LAPACK_LIBRARY}"
  "${BLAS_LIBRARY}"
  "${CFS_FORTRAN_LIBS}"
)

# the pardiso library can be empty, as mkl 12 does not need a solver library
IF(PARDISO_LIBRARY)
  ADD_DEPENDENT_LIBRARY(pardiso-lib "${PARDISO_LIBRARY}" "${DEPENDS}")
  SET(TARGET_LL
    ${TARGET_LL}
    pardiso-lib1)
ENDIF(PARDISO_LIBRARY)

TARGET_LINK_LIBRARIES(pardiso-olas ${TARGET_LL})