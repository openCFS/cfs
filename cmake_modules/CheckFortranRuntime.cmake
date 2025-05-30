#============================================================================r
# The first part of this file is dedicated to finding and copying the Fortran 
# runtime libraries of our various supported compilers. These are needed, since 
# we link from C++ to Fortran libraries out of CFSDEPS which have been built
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
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  #---------------------------------------------------------------------------
  # Let gfortran print its internal search directories, so that we can use
  # them to find its runtime libs.
  #---------------------------------------------------------------------------
  EXECUTE_PROCESS(
    COMMAND "${CMAKE_Fortran_COMPILER}" -print-search-dirs
    WORKING_DIRECTORY "${CFS_BINARY_DIR}"
    OUTPUT_VARIABLE GFORTRAN_SEARCH_DIRS
    RESULT_VARIABLE RETVAL )

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

  # message(STATUS "using GFORTRAN_SEARCH_DIRS=${GFORTRAN_SEARCH_DIRS}")
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
    NO_CMAKE_SYSTEM_PATH)
  MARK_AS_ADVANCED(GFORTRAN_LIBRARY)

  find_library(QUADMATH_LIBRARY
    NAMES quadmath
    PATHS ${GFORTRAN_SEARCH_DIRS}
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH)
  mark_as_advanced(QUADMATH_LIBRARY)
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
    NO_CMAKE_SYSTEM_PATH)
  MARK_AS_ADVANCED(QUADMATH_LIBRARY_STATIC)

  FIND_LIBRARY(GFORTRAN_LIBRARY_STATIC
    NAMES libgfortran${CMAKE_STATIC_LIBRARY_SUFFIX}
    PATHS ${GFORTRAN_SEARCH_DIRS}
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH)
  MARK_AS_ADVANCED(GFORTRAN_LIBRARY_STATIC)

  # at least for macOS Ventura (13) we need to link libgcc (e.g. /usr/local/gfortran/lib/gcc/aarch64-apple-darwin22/12.2.0/libgcc.a)
  # otherwise stuff like ___aarch64_cas8_sync, ___divtf3, ___unordtf2, ... is missing for libgfortran.a
  find_library(GFORTAN_LIBGCC_STATIC
    NAMES libgcc${CMAKE_STATIC_LIBRARY_SUFFIX}
    PATHS ${GFORTRAN_SEARCH_DIRS}
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH)
  mark_as_advanced(GFORTAN_LIBGCC_STATIC)

  # add gfortran
  if(NOT GFORTRAN_LIBRARY_STATIC MATCHES "NOTFOUND")
    list(APPEND CFS_FORTRAN_LIBS "${GFORTRAN_LIBRARY_STATIC}")
    # for apple add libgcc.a, if you have the above missing link issues, check this option also for other systems
    if(APPLE AND NOT GFORTAN_LIBGCC_STATIC MATCHES "NOTFOUND")
      list(APPEND CFS_FORTRAN_LIBS "${GFORTAN_LIBGCC_STATIC}")
    endif()
    # clang complains about -static-libgfortran"
    if(NOT(CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
      list(APPEND CFS_FORTRAN_LIBS "-static-libgfortran") 
    endif()
  elseif(NOT GFORTRAN_LIBRARY MATCHES "NOTFOUND")
    message(STATUS "Warnig: no static 'gfortran' library found - this build might not run on systems missing compatible Fortran runtime libs!")
    list(APPEND CFS_FORTRAN_LIBS "${GFORTRAN_LIBRARY}")
  else()
    message(FATAL_ERROR "neiter static nor dynmaic 'gfrotran' library found in GFORTRAN_SEARCH_DIRS=${GFORTRAN_SEARCH_DIRS}")
  endif()

  # now add quadmath
  if(NOT QUADMATH_LIBRARY_STATIC MATCHES "NOTFOUND")
    list(APPEND CFS_FORTRAN_LIBS "${QUADMATH_LIBRARY_STATIC}")
  elseif(NOT QUADMATH_LIBRARY MATCHES "NOTFOUND")
    message(STATUS "Warnig: no static 'quadmath' library found - this build might not run on systems missing compatible Fortran runtime libs!")
    if(CFS_ARCH MATCHES "X86_64")
    # gfortan for arm64 seems to not know 128 bits https://github.com/Homebrew/homebrew-core/issues/73949
      list(APPEND CFS_FORTRAN_LIBS "${QUADMATH_LIBRARY}")
    endif()
  else()
    message(FATAL_ERROR "neiter static nor dynmaic 'quadmath' library found in GFORTRAN_SEARCH_DIRS=${GFORTRAN_SEARCH_DIRS}")
  endif()

  set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "") # disable dynamic linking of gfortran and quadmath caused by this variable
endif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

# Lets find the runtime libraries of the Intel Fortran compiler.
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel") # ifort and ifx
  # we determine the location for Intel redistributables from the Intel Fortran compiler (ifort and ifx)  

  # An alternative would be to use CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES to find relevant libs, 
  # use find_library on the libs of interest, use get_filename_component() to get the path, and then possibly 
  # create the static and dynamic variant

  # compiler_dir: 
  # ifx: C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\bin
  #      /opt/intel/oneapi/compiler/2023.2.1/linux/bin/ifx
  # ifort: C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\bin\intel64
  #       /opt/intel/oneapi/compiler/2023.2.1/linux/bin/intel64/ifort
  get_filename_component(INTEL_COMPILER_DIR ${CMAKE_Fortran_COMPILER} PATH)    

  # Windows redistributable: C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\redist\intel64_win\compiler
  # Linux libs: /opt/intel/oneapi/compiler/2023.2.1/linux/compiler/lib/intel64  
  if(INTEL_COMPILER_DIR MATCHES "intel64")
    set(INTEL_BASE "${INTEL_COMPILER_DIR}/../..") # ifort 
  else()
    set(INTEL_BASE "${INTEL_COMPILER_DIR}/..") # ifx 
  endif()
  
  # we statically link only on Linux/legacy macOS with Intel 
  if(UNIX)
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 2024)
      set(INTEL_LINK_DIR "${INTEL_BASE}/compiler/lib/intel64")
    else()
      #cmake_print_variables(INTEL_COMPILER_DIR)
      set(INTEL_LINK_DIR "${INTEL_COMPILER_DIR}/../lib")
    endif()

    # this are the libs we link statically on Windows and UNIX
    set(INTEL_LINK_LIBS ifport ifcore imf svml irc) # one might want to check if all are still necesarray? 

    foreach(LIB IN LISTS INTEL_LINK_LIBS)
      set(FILE ${INTEL_LINK_DIR}/lib${LIB}.a)
      if(NOT EXISTS ${FILE})
        message(FATAL_ERROR "cannot find ${FILE}, check code and possibly link shared lib?!") 
      endif()
      list(APPEND CFS_FORTRAN_LIBS ${FILE})
    endforeach()
  endif() # UNIX  

  if(WIN32)
    # Copy compiler runtime .dlls to bin directory
    set(INTEL_DLLS libifcoremd.dll libifportmd.dll libiomp5md.dll libmmd.dll svml_dispmd.dll)
    if(DEBUG)
      list(APPEND INTEL_DLLS libifcoremdd.dll libmmdd.dll)
    endif()

    set(IFORT_DLL_PATHS
      "$ENV{ONEAPI_ROOT}/compiler/latest/bin"
      "${INTEL_COMPILER_DIR}"
      "${INTEL_BASE}/redist/intel64_win/compiler"
    )
  
    message(STATUS "searching in IFORT_DLL_PATHS=${IFORT_DLL_PATHS}")
   
    
    # Windows requires the libs in bin
    # message(STATUS "Copying Intel redistributable files from ${INTEL_REDIST_DIR} to ${CFS_BINARY_DIR}/bin")
    foreach(lib IN LISTS INTEL_DLLS)
      find_file(${lib}_path NAMES "${lib}" PATHS ${IFORT_DLL_PATHS} NO_CACHE)
      if(NOT "${${lib}_path}" MATCHES "NOTFOUND")
        install(PROGRAMS ${${lib}_path} DESTINATION "bin")
        message(STATUS "  ${lib} will be installed from ${${lib}_path} via cpack")
      else()
        message(WARNING "  ${lib} not found!")
      endif()
      #file(INSTALL ${INTEL_REDIST_DIR}/${LIB} DESTINATION ${CFS_BINARY_DIR}/bin)
    endforeach()
  endif() # WIN32
endif() # Intel Fortran compilers ifort and ifx

#==============================================================================
# Create a header file include/def_cfs_fortran_interface.hh with the correct
# mangling for the Fortran routines called in openCFS.
# 
# How LAPACK library enables Microsoft Visual Studio support with CMake
# and LAPACKE:
# http://www.netlib.org/lapack/lawnspdf/lawn270.pdf and
#
# Fortran for C/C++ developers made easier with CMake:
# http://www.kitware.com/blog/home/post/231
#==============================================================================

include(FortranCInterface)

# This stuff is important such that the linker finds the functions if the names
# are different (uppercase, tailing underline, ...).
# Each fortran function which are used in CFS have to be listed here else one
# will get an undefined reference to 'symbol' error.
# Also make sure the function is defined in MatVec/BLASLAPACKInterface.hh
#
# the cmake module FortranCInterface_HEADER() is kind of sensitive, therefore
# we provide def_cfs_fortran_interface_fallback.hh. Extend this file if you
# modifiy the list below

include("${CFS_SOURCE_DIR}/cmake_modules/FortranManglingModules.cmake")
 
# in case cmake name mangling fails our reference might work. Only because cmake is confused does not
# mean that it would not for the compiler! A hint for this issue, is to use the proper cmake (e.g. architecture)
if(NOT FortranCInterface_GLOBAL_FOUND OR NOT FortranCInterface_GLOBAL_FOUND)
  configure_file("${CFS_SOURCE_DIR}/include/def_cfs_fortran_interface_fallback.hh" "${CFS_BINARY_DIR}/include/def_cfs_fortran_interface.hh" COPYONLY)
  message(STATUS "because cmake's Fortran name mangling failed, we use def_cfs_fortran_interface_fallback.hh")
endif()
