INCLUDE_DIRECTORIES(${CHOLMOD_INCLUDE_DIR})

IF(EXPLICIT_TEMPLATE_INSTANTIATION)
    SET(CHOLMOD_SRCS CholMod.cc)
ELSE(EXPLICIT_TEMPLATE_INSTANTIATION)
    SET(CHOLMOD_SRCS CholModDummy.cc)
ENDIF(EXPLICIT_TEMPLATE_INSTANTIATION)

#SET_SOURCE_FILES_PROPERTIES(${ILUPACK_SRCS}
#  PROPERTIES COMPILE_FLAGS "-fno-strict-aliasing")

ADD_LIBRARY(cholmod-olas STATIC ${CHOLMOD_SRCS})

# CFS_FORTRAN_LIBS has been set in cmake_modules/distro.cmake
# LAPACK_LIBRARY and BLAS_LIBRARY are defined in
# cmake_modules/FindFortranLibs.cmake
TARGET_LINK_LIBRARIES(cholmod-olas
  cfsgeneral
  datainout
  logging
  matvec
  paramh
  solver-olas
  cholmod-lib1
)

SET(DEPENDS
  "${LAPACK_LIBRARY}"
  "${BLAS_LIBRARY}"
  "${METIS_LIBRARY}"
)

ADD_DEPENDENT_LIBRARY(cholmod-lib "${CHOLMOD_LIBRARY}" "${DEPENDS}")
