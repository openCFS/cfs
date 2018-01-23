#---------------------------------------------------
# The ADOLC Library
#
#Project Homepage
# http://projects.coin-or.org/ADOL-C
#---------------------------------------------------
# 
SET(ADOLC_FOUND 0)

#-------------------------------------------------------------------------------
# Look for ADOLC header.
#-------------------------------------------------------------------------------

#BUILD_EXTLIB("SGpp"
#  "${CFS_BINARY_DIR}/include/sgpp_base.hpp"
#  "${CFS_DEPS_ROOT}/sgpp/build_sgpp.pl"
#  "build_sgpp.log")

#SET(HOME "/home/bjurgel")  # this needs to be set manually
SET(HOME "~")
#SET(ADOLC_INCLUDE_DIR "${HOME}/adolc_base/include")
SET(ADOLC_INCLUDE_DIR "${CFS_SOURCE_DIR}/include")
MESSAGE( STATUS "bjurgel comment: ADOLC_INCLUDE_DIR is:   " ${ADOLC_INCLUDE_DIR}  )
#MESSAGE( STATUS "bjurgel comment: CFSHOME is:   " ${CFS_HOME}  )

#SET(ADOLC_INCLUDE_DIR ${ADOLC_INCLUDE_DIR} "/home/users/bjurgel/adic_adolc_interface/adic_runtime")

#-------------------------------------------------------------------------------
# Determine paths of ADOLC libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")

SET(ADOLC_LIBRARY
  #"${LD}/libadolc.so" # die Statische version *.a funktioniert nicht und wird in der Zukunft von ADOL-C nicht mehr erzeugt
  "${HOME}/adolc_base/lib64/libadolc.so"
  CACHE FILEPATH "ADOLC library.")
MARK_AS_ADVANCED(ADOLC_LIBRARY)
SET(ADIC_RUNTIME_LIBRARY
  "${LD}/libruntimedense.so"
  CACHE FILEPATH "ADIC_RUNTIME_LIBRARY")
MARK_AS_ADVANCED(ADIC_RUNTIME_LIBRARY)



MESSAGE( STATUS "bjurgel comment: ADOLC_LIBRARY is:   " ${ADOLC_LIBRARY}  )

# Add a switch for the compiler.
# As the ADOLC toolbox is not available for everyone, we have to exclude the corresponding source code from compiling.
IF(USE_ADOLC)
  ADD_DEFINITIONS(-DUSE_ADOLC)
ENDIF()

SET(ADOLC_FOUND 1)
