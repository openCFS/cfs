#-------------------------------------------------------------------------------
# The GNU MP Bignum Library
# Needed by MPFR and CGAL
#
# Project Homepage
# http://gmplib.org
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to gmp sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(gmp_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/gmp")
set(gmp_source  "${gmp_prefix}/src/gmp")
set(gmp_install  "${CMAKE_CURRENT_BINARY_DIR}")


if(CFS_ARCH STREQUAL "I386")
  SET(ABI "32")
else(CFS_ARCH STREQUAL "I386")
  SET(ABI "64")
endif(CFS_ARCH STREQUAL "I386")

SET(CONF_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/gmp/gmp-configure.sh.in")
SET(CONF "${gmp_prefix}/gmp-configure.sh")
CONFIGURE_FILE("${CONF_TEMPL}" "${CONF}" @ONLY) 

#-------------------------------------------------------------------------------
# The gmp external project
#-------------------------------------------------------------------------------
ExternalProject_Add(gmp
  PREFIX "${gmp_prefix}"
  SOURCE_DIR "${gmp_source}"
  URL ${GMP_URL}/${GMP_BZ2}
  URL_MD5 ${GMP_MD5}
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND sh ${CONF}
  BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} -f Makefile
  INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} -f Makefile install
)

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/gmp/gmp-patch.cmake.in")
SET(PFN "${gmp_prefix}/gmp-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# We do not use the PATCH_COMMAND  of ExternalProject_Add since we do not only
# want to apply the patch script  during configuration time but also if it has
# changed.  Therefore,   we  need  a   dependency  on  the   configured  patch
# script. This can be achieved by  adding an additional build step between the
# download and configure steps.
#
# NOTE: The  patch script should  be designed  in such a  way, that it  can be
# applied to  an already patched  source tree. This  is due to the  fact, that
# ExternalProject_Add only extracts the source if the MD5 sum has has changed.
#-------------------------------------------------------------------------------
ExternalProject_Add_Step(gmp custom_patch
   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
   DEPENDEES download
   DEPENDERS configure
   DEPENDS "${PFN}"
   WORKING_DIRECTORY ${gmp_source}
)


SET(GMP_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-----------------------------------------------------------------------------
# Determine paths of GMP libraries.
#-----------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(GMP_LIBRARY
  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}gmp${CMAKE_STATIC_LIBRARY_SUFFIX}"
  CACHE FILEPATH "GMP library.")

MARK_AS_ADVANCED(GMP_LIBRARY)
MARK_AS_ADVANCED(GMP_INCLUDE_DIR)

