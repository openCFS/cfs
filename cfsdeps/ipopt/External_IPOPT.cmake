#-------------------------------------------------------------------------------
# IPOPT (interior point optimizer) is a general purpose open source optimizer
#
# IPOPT needs a solver. This can be HSL, which can be requested for academic purpose
# for free or MKL-Pardiso.
# * HSL is provied encrypted via ipopt_hsl.zip and requires CFS_KEY_IPOPT_HSL
# * MKL pardiso is only recently possible for IPOPT and not used here. You might 
#   change this file to use MKL pardiso!
#-------------------------------------------------------------------------------


set(IPOPT_PREFIX  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/ipopt")
set(IPOPT_SOURCE  "${IPOPT_PREFIX}/src/ipopt")
# we use this as temporary install directory to alllow packing for precompiled cfsdeps
set(IPOPT_INSTALL  "${IPOPT_PREFIX}/install")

# snopt is build by configure, therefore no cmake args needed

SET(MIRRORS
  "http://www.coin-or.org/download/source/Ipopt/${IPOPT_TGZ}"
  "${CFS_DS_SOURCES_DIR}/ipopt/${IPOPT_TGZ}"
)

SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/ipopt/${IPOPT_TGZ}")
SET(MD5_SUM ${IPOPT_MD5})

# precompiled cfsdeps not yet implemented!

# the encrypted ipopt_hsl.zip is so small and for everyone available on the web if you register, we provide it
# in the cfs/cfsdeps/ipopt itself. You just need the key.
if(NOT CFS_KEY_IPOPT_HSL)
  message(FATAL_ERROR "Key for encrypted ipopt_hsl.zip required in CFS_KEY_IPOPT_HSL. E.g. set in your ~/.cfs_platform_defaults.cmake")
endif() 
# The patch triggers downloading netlib blas and lapack such that ipopt compiles with blas support.
# However netlib blas will be removed after installation and the cfs stuff will be used
# Furthermore ipopt_hsl is extracted 
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/ipopt/ipopt-patch.cmake.in")
SET(PFN "${IPOPT_PREFIX}/ipopt-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

# the whole build projectll
ExternalProject_Add(ipopt
  PREFIX "${IPOPT_PREFIX}"
  SOURCE_DIR "${IPOPT_SOURCE}"
  # the stepp cfsdeps_download ensures the file will be found in cfsdepscache
  URL ${LOCAL_FILE}
  URL_MD5 ${IPOPT_MD5}
  # note that when unzipping the source, the Ipopt-3.11.9 directory is omitted
  BUILD_IN_SOURCE 1
  # handle mirror
  # Fill ThirdParty content (blas/ lapack/ HSLold)
  PATCH_COMMAND  ${CMAKE_COMMAND} -P "${PFN}"
  # let it install to the temporay directory where we can remove libcoinblas and libcoinlapack and prepare to copy to precompiled cfsdeps
  CONFIGURE_COMMAND ${IPOPT_SOURCE}/configure --prefix=${IPOPT_INSTALL} --libdir=${IPOPT_INSTALL}/lib64/${CFS_ARCH_STR} --disable-shared --disable-linear-solver-loader --with-metis-lib=${METIS_LIBRARY} --with-metis-incdir=${CMAKE_CURRENT_BINARY_DIR}/include --disable-pkg-config   
)

# this handles the mirrors and cfsdepscache source stuff.  
SET(DLFN "${IPOPT_PREFIX}/ipopt-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

# doing cfsdescache stuff such that the default external project download takes the file from cfsdescache
ExternalProject_Add_Step(ipopt cfsdeps_download
   COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
   DEPENDERS download
   DEPENDS "${DLFN}"
   WORKING_DIRECTORY ${ipopt_prefix}
)


# after the installation we remove not necessary libs, copy to cfs and might copy to precompiled cfsdeps
SET(PI_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/ipopt/ipopt-post_install.cmake.in")
SET(PI "${IPOPT_PREFIX}/ipopt-post_install.cmake")
CONFIGURE_FILE("${PI_TEMPL}" "${PI}" @ONLY) 

# execute the stuff from snopt-post_install.cmake after installation
ExternalProject_Add_Step(ipopt post_install COMMAND ${CMAKE_COMMAND} -P "${PI}" DEPENDEES install)

# Add project to global list of CFSDEPS, this allows "make snopt"
SET(CFSDEPS ${CFSDEPS} ipopt)

# Determine paths of IPOPT libraries. We have only libs include/coin
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(IPOPT_LIBRARY "${LD}/libipopt.a;${LD}/libcoinhsl.a" CACHE FILEPATH "IPOPT library.")
SET(IPOPT_INCLUDE_DIR "${CFS_BINARY_DIR}/include" CACHE FILEPATH "IPOPT library.")

MARK_AS_ADVANCED(IPOPT_LIBRARY)
MARK_AS_ADVANCED(IPOPT_INCLUDE_DIR)
