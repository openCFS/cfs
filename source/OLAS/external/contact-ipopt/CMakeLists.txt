INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})

# check whether this is all needed
#IF(EXPLICIT_TEMPLATE_INSTANTIATION)
SET(IPOPT_SRCS ContactIpopt.cc)
#ELSE(EXPLICIT_TEMPLATE_INSTANTIATION)
#    SET(IPOPT_SRCS IpoptDummy.cc)
#ENDIF(EXPLICIT_TEMPLATE_INSTANTIATION)

ADD_LIBRARY(contact-ipopt STATIC ${IPOPT_SRCS})

# CFS_FORTRAN_LIBS has been set in cmake_modules/distro.cmake
# LAPACK_LIBRARY and BLAS_LIBRARY are defined in
# cmake_modules/FindFortranLibs.cmake
TARGET_LINK_LIBRARIES(contact-ipopt
  cfsgeneral
  datainout
  logging
  matvec
  paramh
  solver-olas
  ipopt-lib1
)

# the same thing is already done in optimization, however, cmake does know nothing about that here
SET(COMBINED ${IPOPT_LIBRARY} ${IPOPT_HSL_LIBRARY} ${METIS_LIBRARY})
ADD_DEPENDENT_LIBRARY(ipopt-lib "${COMBINED}" "${DEPENDS}" )
