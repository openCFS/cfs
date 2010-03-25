# - Find SnOpt installation 
# This module finds if SnOpt is installed and determines where 
# the libraries are.
#
#  SNOPT_FOUND   = system has SnOpt 
#  SNOPT_LIBRARY = path to the SnOpt libraries
#
# Fabian Schury (09/09)
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SET(SNOPT_FOUND 0)

SET(F2C "/usr/bin" CACHE FILEPATH "Path to f2c binary.")
MARK_AS_ADVANCED(F2C)

SET(F2CINCLUDE "/usr/include" CACHE FILEPATH "Path to f2c include directory.")
MARK_AS_ADVANCED(F2CINCLUDE)

SET(F2CLIBRARY "/usr/lib" CACHE FILEPATH "Path to f2c libraries.")
MARK_AS_ADVANCED(F2CLIBRARY)

#-------------------------------------------------------------------------------
# Look for f2c variables
#-------------------------------------------------------------------------------
CONFIGURE_FILE("${CFS_DEPS_ROOT}/snopt/build_snopt_vars.pl.in"
  "${CFS_TEMP_DIR}/build_snopt_vars.pl"
  @ONLY )


BUILD_EXTLIB("SnOpt"
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libsnopt_c.a"
  "${CFS_DEPS_ROOT}/snopt/build_snopt.pl"
  "build_snopt.log")

#-------------------------------------------------------------------------------
# Determine paths of SNOPT libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(SNOPT_LIBRARY "${LD}/libsnopt_c.a;${LD}/libsnopt-blas_c.a;${LD}/libsnprint_c.a;libf2c.a;-lm" 
    CACHE FILEPATH "SnOpt library.")

MARK_AS_ADVANCED(SNOPT_LIBRARY)

SET(SNOPT_FOUND 1)
