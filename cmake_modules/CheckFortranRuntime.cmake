#=============================================================================
# The first part of this file is dedicated to finding the Fortran runtime
# libraries of our various supported compilers. These are needed, since we
# link from C++ to Fortran libraries out of CFSDEPS which have been built
# using these Fortran compiler.
#
# The second part is dedicated to obtaining informations on the name mangling
# that is used by the Fortran compiler. Usually this task can be handled by
# the macros in FortranCInterface.cmake which is provided by CMake >= 2.8.7.
# In some cases (e.g. while cross-compiling) we have to specify the name
# mangling definitions header by hand.
#=============================================================================

#-----------------------------------------------------------------------------
# Initialize CFS_FORTRAN_LIBS variable.
#-----------------------------------------------------------------------------
SET(CFS_FORTRAN_LIBS "")

#-----------------------------------------------------------------------------
# Lets find the runtime libraries of gfortran. This branch should be taken
# on most Linuxes and MacOS.
#-----------------------------------------------------------------------------
if(CFS_FORTRAN_COMPILER_NAME STREQUAL "GNU")
  message(STATUS "GNU Fortran")
  #---------------------------------------------------------------------------
  # On Unix we may encounter the situation, that we have Intel as Fortran
  # compiler but some system libs depend on GFortran. Therefore, we assume,
  # that the system GFortran compiler is on the PATH and is named gfortran.
  #---------------------------------------------------------------------------
  IF(NOT CFS_FORTRAN_COMPILER_NAME STREQUAL "GNU")
    SET(GFORTRAN_EXECUTABLE "gfortran")
  ELSE(NOT CFS_FORTRAN_COMPILER_NAME STREQUAL "GNU")
    SET(GFORTRAN_EXECUTABLE "${CMAKE_Fortran_COMPILER}")
  ENDIF(NOT CFS_FORTRAN_COMPILER_NAME STREQUAL "GNU")

  #---------------------------------------------------------------------------
  # Let gfortran print its internal search directories, so that we can use
  # them to find its runtime libs.
  #---------------------------------------------------------------------------
  EXECUTE_PROCESS(
    COMMAND "${GFORTRAN_EXECUTABLE}" -print-search-dirs
    WORKING_DIRECTORY "${CFS_BINARY_DIR}"
    OUTPUT_VARIABLE GFORTRAN_SEARCH_DIRS
    RESULT_VARIABLE RETVAL
    )

  #---------------------------------------------------------------------------
  # Prepare search paths for input to the FIND_LIBRARY function
  #---------------------------------------------------------------------------
  SET(DIRSEP ":")
  STRING(REPLACE ";" "#____#" GFORTRAN_SEARCH_DIRS "${GFORTRAN_SEARCH_DIRS}")
  STRING(REPLACE "\n" ";" SEARCH_DIRS "${GFORTRAN_SEARCH_DIRS}")
  IF(WIN32)
    SET(DIRSEP "#____#")
  ENDIF()

  foreach(line IN LISTS SEARCH_DIRS)
    IF(line MATCHES "libraries")
      STRING(REPLACE "=" ";" line "${line}")
      LIST(GET line 1 line)
      STRING(REPLACE ${DIRSEP} ";" GFORTRAN_SEARCH_DIRS "${line}")
     # MESSAGE(FATAL_ERROR "GFORTRAN_SEARCH_DIRS ${GFORTRAN_SEARCH_DIRS}")
    ENDIF()
  endforeach()

  IF(NOT RETVAL EQUAL 0)
    SET(GFORTRAN_SEARCH_DIRS "")
  ENDIF(NOT RETVAL EQUAL 0)

  #---------------------------------------------------------------------------
  # Let's find the shared version of the gfortran runtime lib.
  #---------------------------------------------------------------------------
  FIND_LIBRARY(GFORTRAN_LIBRARY
    NAMES gfortran 
    PATHS ${GFORTRAN_SEARCH_DIRS}
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH 
  )

  MARK_AS_ADVANCED(GFORTRAN_LIBRARY)

  #---------------------------------------------------------------------------
  # Now, let's see if we can also find the static runtime libs.
  #---------------------------------------------------------------------------
  FIND_LIBRARY(QUADMATH_LIBRARY_STATIC
    NAMES libquadmath${CMAKE_STATIC_LIBRARY_SUFFIX}
    PATHS ${GFORTRAN_SEARCH_DIRS}
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH
  )
  #message("QUADMATH_LIBRARY_STATIC=${QUADMATH_LIBRARY_STATIC}")
  MARK_AS_ADVANCED(QUADMATH_LIBRARY_STATIC)

  FIND_LIBRARY(GFORTRAN_LIBRARY_STATIC
    NAMES libgfortran${CMAKE_STATIC_LIBRARY_SUFFIX}
    PATHS ${GFORTRAN_SEARCH_DIRS}
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH
  )
  #message("GFORTRAN_LIBRARY_STATIC=${GFORTRAN_LIBRARY_STATIC}")
  MARK_AS_ADVANCED(GFORTRAN_LIBRARY_STATIC)

  #---------------------------------------------------------------------------
  # Prefer static runtime libs over shared ones.
  #---------------------------------------------------------------------------
  if(NOT GFORTRAN_LIBRARY_STATIC MATCHES "NOTFOUND" AND NOT QUADMATH_LIBRARY_STATIC MATCHES "NOTFOUND")
    
    LIST(APPEND CFS_FORTRAN_LIBS "${GFORTRAN_LIBRARY_STATIC}" "${QUADMATH_LIBRARY_STATIC}")
    
    # clang on macOS complains about -static-libgfortran"
    if(NOT(APPLE AND CFS_CXX_COMPILER_NAME STREQUAL "CLANG"))
      LIST(APPEND CFS_FORTRAN_LIBS "-static-libgfortran") 
    endif()  
  else()
    message(WARNING "no static Fortran library found - this build will not run on systems missing compatible Fortran runtime libs!")
    LIST(APPEND CFS_FORTRAN_LIBS "${GFORTRAN_LIBRARY}" "-lquadmath")
  endif()
endif(CFS_FORTRAN_COMPILER_NAME STREQUAL "GNU")

#-----------------------------------------------------------------------------
# Lets find the runtime libraries of the Intel Fortran compiler.
#-----------------------------------------------------------------------------
IF(CFS_FORTRAN_COMPILER_NAME STREQUAL "IFORT")
  message(STATUS "Intel Fortran compiler")
  #---------------------------------------------------------------------------
  # Search explicitely for implicitely defined Fortran libs
  # CMake-magic: https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_IMPLICIT_LINK_LIBRARIES.html
  #---------------------------------------------------------------------------
  foreach(lib IN LISTS CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)
    FIND_LIBRARY(IFORT_${lib}_LIBRARY
      NAMES ${lib}
      PATHS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}
      NO_DEFAULT_PATH
      NO_CMAKE_ENVIRONMENT_PATH
      NO_CMAKE_PATH
      NO_SYSTEM_ENVIRONMENT_PATH
      NO_CMAKE_SYSTEM_PATH 
      )
    MARK_AS_ADVANCED(IFORT_${lib}_LIBRARY)
    #message(STATUS "${lib}")
    # copy over intel shared libs
    if("${lib}" MATCHES "(ifport)|(ifcoremt)|(imf)|(svml)|(irc)$")
      # the dynamic version is usually in the same dir as the static one
      get_filename_component(inteldir "${IFORT_${lib}_LIBRARY}" DIRECTORY)
      FIND_LIBRARY(IFORT_${lib}_LIBRARY_STATIC
        NAMES lib${lib}${CMAKE_STATIC_LIBRARY_SUFFIX}
        PATHS ${inteldir}
        NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH
      	)
      mark_as_advanced(IFORT_${lib}_LIBRARY_STATIC)
      if("${IFORT_${lib}_LIBRARY_STATIC}" MATCHES "NOTFOUND")
	message(WARNING "no static version of ${IFORT_${lib}_LIBRARY} found in ${inteldir} using dynamic version.\n
	It will be installed by cpack, but not copied over into the build directory.")
        #one can't simply copy over because cmake will complain about "Cannot generate a safe linker search path ..."
        install(FILES "${IFORT_${lib}_LIBRARY}" DESTINATION "lib64/LINUX_${CFS_ARCH}")
        list(APPEND CFS_FORTRAN_LIBS "${IFORT_${lib}_LIBRARY}")
      else()
	message(STATUS "  found static lib: ${IFORT_${lib}_LIBRARY_STATIC} and using it.")
	list(APPEND CFS_FORTRAN_LIBS "${IFORT_${lib}_LIBRARY_STATIC}")
      endif()
    endif()
  endforeach()
ENDIF(CFS_FORTRAN_COMPILER_NAME STREQUAL "IFORT")

# support for the Open64 Fortran compiler removed. Check svn version 15997
#==============================================================================
# Create a header file include/def_cfs_fortran_interface.hh with the correct
# mangling for the Fortran routines called in CFS++.
# 
# How LAPACK library enables Microsoft Visual Studio support with CMake
# and LAPACKE:
# http://www.netlib.org/lapack/lawnspdf/lawn270.pdf and
#
# Fortran for C/C++ developers made easier with CMake:
# http://www.kitware.com/blog/home/post/231
#==============================================================================

message(STATUS "setting CFS_FORTRAN_LIBS=${CFS_FORTRAN_LIBS}")

include(FortranCInterface)

# This stuff is important such that the linker finds the functions if the names
# are different (uppercase, tailing underline, ...).
# Each fortran function which are used in CFS have to be listed here else one
# will get an undefined reference to 'symbol' error.
# Also make sure the function is defined in MatVec/BLASLAPACKInterface.hh
FortranCInterface_HEADER("${CFS_BINARY_DIR}/include/def_cfs_fortran_interface.hh"
  MACRO_NAMESPACE "CFS_FORTRAN_INTERFACE_"
  SYMBOLS
  # BLAS and LAPACK
  dgecon
  dgemm
  zgemm
  dgemv
  zgemv
  dgetrf
  dgetri
  dpotrf
  dpotri
  dposv
  zsysv
  zhesv
  zgesv
  zheev
  dgbtrf
  zgbtrf
  dgbequ
  zgbequ
  dlaqgb
  zlaqgb
  dgbtrs
  zgbtrs
  dgbrfs
  zgbrfs
  dlartg
  zlartg
  dpbtrf
  zpbtrf
  dpbtrs
  zpbtrs
  dsyev
  zgetrf
  zgetri
  ilaver
  cblas_zdscal
  cblas_dznrm2
  # ARPACK
  dsaupd
  dseupd
  debug
  znaupd
  zneupd
  # Pardiso
  pardisoinit
  pardiso
  # additionally for PALM
  dscal
  dgesvd
  # scpip
  scpip30
  # snopt
  snopta
  sninit
  sngeti
  sngetr
  snset
  sngetc
  snseti
  snsetr
  snspec
  snmema
  snjac
  snopenappend
  snclose
)

