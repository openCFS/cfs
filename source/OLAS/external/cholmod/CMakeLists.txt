INCLUDE_DIRECTORIES(${SUITESPARSE_INCLUDE_DIR})

SET(CHOLMOD_SRCS CholMod.cc)

ADD_LIBRARY(cholmod-olas STATIC ${CHOLMOD_SRCS})

# CFS_FORTRAN_LIBS has been set in cmake_modules/distro.cmake
# LAPACK_LIBRARY and BLAS_LIBRARY are defined in
# cmake_modules/FindFortranLibs.cmake
set(TARGET_LL ${SUITESPARSE_LIBRARY} ) # has also umpfpack what we don't need here

SET(TARGET_LL
  ${TARGET_LL}
  ${LAPACK_LIBRARY}
  ${BLAS_LIBRARY}
  ${METIS_LIBRARY}
  )

TARGET_LINK_LIBRARIES(cholmod-olas ${TARGET_LL} )

add_dependencies(cholmod-olas suitesparse)
