# In this file we set the default configuration for open CFS
# It is a good idea always to set the *_DEFAULT, because setting the 
# actual variable (e.g as in the main cmake file), a warning
# "Policy CMP0077 is not set: option() honors normal variables" can appear
# when a cache variable is set there but the variable already exists. 

# This settings can be overwritten in the following ways:
# - calling cmake -D<variable>_DEFAULT=*
# - configuring a $HOME/.cfs_platform_defaults.cmake (read at the end)
# - on Windows from platform_defaults_win.cmake (read at the end)
# - on specific hosts from cfs_platform_defaults_${CFS_BUILD_HOST}.cmake
# - configuring a $HOME/.cfs_platform_defaults_${CFS_BUILD_HOST}.cmake
# - simply manually via ccmake

set(BUILD_TESTING_DEFAULT OFF)
set(BUILD_CFSTOOL_DEFAULT ON)
set(BUILD_CFSDAT_DEFAULT ON)
set(BUILD_UNIT_TESTS_DEFAULT ON)

set(DEBUG_DEFAULT OFF)

set(TESTSUITE_DIR_DEFAULT "${CFS_SOURCE_DIR}/Testsuite")
set(CFS_DEPS_CACHE_DIR_DEFAULT "${CFS_BINARY_DIR}/cfsdeps/cache")
set(CFS_DEPS_PRECOMPILED_DEFAULT ON)
set(CFS_NATIVE_DEFAULT OFF)
set(CFS_REORDERING_DEFAULT "default")

set(USE_GIDPOST_DEFAULT ON)
set(USE_GMSH_DEFAULT ON)
set(USE_ENSIGHT_DEFAULT ON)

set(USE_BLAS_LAPACK_DEFAULT "MKL")
set(CFS_PARDISO_DEFAULT "MKL")

# the established one, libxml2 is lighter
set(USE_XML_READER_DEFAULT "xerces")

set(USE_PARDISO_DEFAULT ON)
set(USE_ARPACK_DEFAULT ON)
set(USE_EMBEDDED_PYTHON_DEFAULT OFF) # not easily portable binaries
set(USE_PHIST_CG_DEFAULT OFF) # phist and ghost are difficult to build
set(USE_PHIST_EV_DEFAULT OFF)
set(BUILD_GHOST_DEFAULT OFF)
set(USE_FEAST_DEFAULT ON)
set(USE_ILUPACK_DEFAULT OFF) # is closed soure
set(USE_SUITESPARSE_DEFAULT OFF) # is configured as GPL
set(USE_LIS_DEFAULT ON)
set(USE_SUPERLU_DEFAULT ON)
set(USE_CGNS_DEFAULT ON)
set(USE_METIS_DEFAULT ON)

# commercial or condifential code required for the developer building the modules
set(USE_SCPIP_DEFAULT OFF)  
set(USE_SNOPT_DEFAULT OFF)
set(USE_SGPP_DEFAULT OFF) # rarely used
set(USE_IPOPT_DEFAULT ON)
set(USE_CGAL_DEFAULT OFF) # GPL
set(USE_LIBFBI_DEFAULT OFF)
set(USE_FLANN_DEFAULT ON)

set(USE_OPENMP_DEFAULT ON) # for Debug usually OFF

# ----------------------------------------------------------------------
# by the following specific platform_defaults values can be overwritten
# ----------------------------------------------------------------------

if(WIN32)
  include("cmake_modules/platform_defaults_win.cmake")
endif()

set(CMAKE_HOST_DEFAULTS_INC "${CFS_SOURCE_DIR}/cmake_modules/platform_defaults_${CFS_BUILD_HOST}.cmake") 
if(EXISTS "${CMAKE_HOST_DEFAULTS_INC}")
  include("${CMAKE_HOST_DEFAULTS_INC}")
endif()

IF(WIN32)
	set(CMAKE_HOST_DEFAULTS_INC_HOME "$ENV{HOMEDRIVE}/Users/$ENV{USERNAME}/.cfs_platform_defaults.cmake") 
ELSE()
  set(CMAKE_HOST_DEFAULTS_INC_HOME "$ENV{HOME}/.cfs_platform_defaults.cmake") 
ENDIF(WIN32)

# message(DEBUG "local file is ${CMAKE_HOST_DEFAULTS_INC_HOME}") # VERBOSE-TRACE are only supported from cmake 3.15
if(EXISTS "${CMAKE_HOST_DEFAULTS_INC_HOME}")
  message(STATUS "including ${CMAKE_HOST_DEFAULTS_INC_HOME}")
  include("${CMAKE_HOST_DEFAULTS_INC_HOME}")
else()
  message(STATUS "no personal configuration ${CMAKE_HOST_DEFAULTS_INC_HOME} found")
endif()

set(CMAKE_HOST_DEFAULTS_INC_HOME_HOST "$ENV{HOME}/.cfs_platform_defaults_${CFS_BUILD_HOST}.cmake") 
if(EXISTS "${CMAKE_HOST_DEFAULTS_INC_HOME_HOST}")
  include("${CMAKE_HOST_DEFAULTS_INC_HOME_HOST}")
endif()
