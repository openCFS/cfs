INCLUDE_DIRECTORIES(${CFSDEPS_INCLUDE_DIR})

SET(ILUPACK_SRCS Ilupack.cc)

SET_SOURCE_FILES_PROPERTIES(${ILUPACK_SRCS}
  PROPERTIES COMPILE_FLAGS "-fno-strict-aliasing")

ADD_LIBRARY(ilupack-olas STATIC ${ILUPACK_SRCS})

# CFS_FORTRAN_LIBS has been set in cmake_modules/distro.cmake
# LAPACK_LIBRARY and BLAS_LIBRARY are defined in
# cmake_modules/FindFortranLibs.cmake
SET(TARGET_LL
  ${ILUPACK_LIBRARY}
  )

TARGET_LINK_LIBRARIES(ilupack-olas ${TARGET_LL} )

ADD_DEPENDENCIES(ilupack-olas cfsdeps)
ADD_DEPENDENCIES(ilupack-olas ilupack-double)
ADD_DEPENDENCIES(ilupack-olas ilupack-complex)
