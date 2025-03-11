# Use the MKLConfig.cmake file provided by recent oneAPI distributions

# use the "traditional" openCFS policy as default
if(UNIX) # static on Linux
  set(MKL_LINK "static" CACHE STRING "MKL linking model")
else() # dynamic on Windows
  set(MKL_LINK "dynamic" CACHE STRING "MKL linking model")
endif()
# both should work for all platforms
message(STATUS "selected ${MKL_LINK} linking for MKL")

if(NOT USE_OPENMP)
  set(MKL_THREADING "sequential" CACHE STRING "MKL threading layer: intel_thread (Intel OpenMP), gnu_thread (GNU OpenMP)")
else()
  set(MKL_THREADING "intel_thread" CACHE STRING "MKL threading layer: intel_thread (Intel OpenMP), gnu_thread (GNU OpenMP)")
  # this section searches for libiomp5 to save it in OMP_LIBRARY 
  # it should not be necessary at all - maybe it can be removed ...
  # TODO: check if this section is necessay on recent oneMKL installations, with/without setvars environment
  if("${MKL_THREADING}" MATCHES "intel_thread") # at the moment the only choice, 'if' put there to be extensible
    if($ENV{SETVARS_COMPLETED}) # we are in an oneAPI environemnt OMP_LIBRARY
      message(STATUS "Intel oneAPI environemnt detected, iomp5 should be found automatically ...")
    else() # we need to help a bit and set OMP_LIBRARY
      set(OMP_POSSIBLE_PATHS
        "/opt/intel/oneapi/compiler/latest/lib/intel64_lin"
        "/opt/intel/oneapi/compiler/latest/lib"
        )
      #message("MKL_ROOT=${MKL_ROOT} ${LIB_SUFFIX}iomp5${CMAKE_SHARED_LIBRARY_SUFFIX}" )
      find_library(OMP_LIBRARY NAMES "${LIB_SUFFIX}iomp5${CMAKE_SHARED_LIBRARY_SUFFIX}" PATHS ${OMP_POSSIBLE_PATHS})
      if("${OMP_LIBRARY}" MATCHES "NOTFOUND")
        message(STATUS "still ${OMP_LIBRARY} even without our hints. Please search for it on your system (or install it) and improve the hints in OMP_POSSIBLE_PATHS.")
      else()
        message(STATUS "found OMP_LIBRARY=${OMP_LIBRARY}, possibly used by MKLConfig.cmake")
      endif()
      mark_as_advanced(OMP_LIBRARY)
    endif() # threading type
  endif()
endif()

# these need to be chache variables since they are exported as cache variables by MKLConfig.cmake
# we need BLAS and LAPACK from MKL
set(ENABLE_BLAS95 ON CACHE BOOL "Enables BLAS Fortran95 API") # probably could be forced to ON
set(ENABLE_LAPACK95 ON CACHE BOOL "Enables LAPACK Fortran95 API") # probably could be forced to ON

# select integer interface for MKL
set(MKL_INTERFACE "lp64")

# possible paths to search for MKLConfig.cmake
set(MKL_POSSIBLE_PATHS
  # common
  "$ENV{MKLROOT}"
  "$ENV{ONEAPI_ROOT}/mkl/latest"
  # windows
  "C:/Program Files (x86)/Intel/oneAPI/mkl/latest"
  "$ENV{INTEL_INSTALL_DIR}/mkl/latest" # shared runners in ci pipeline
  # linux
  "/opt/intel/oneapi/mkl/latest"
)
find_package(MKL PATHS ${MKL_POSSIBLE_PATHS} NO_DEFAULT_PATH) # finds and runs MKLConfig.cmake

if(MKL_FOUND) # TODO: can be removed once we set find_package(MKL REQUIRED ... ) and remove legacy behaviour
if(WIN32 AND ${MKL_VERSION} VERSION_LESS "2024")
  message(WARNING "MKL version ${MKL_VERSION} < 2024 found on Windows - this might fail with \"don't know how to make 'C:.../libiomp5md.lib'\"")
endif()

message(STATUS "Found MKL version ${MKL_VERSION} via MKL_CONFIG=${MKL_CONFIG}")

# hide advanced cache options from MKLConfig.cmake 
mark_as_advanced(ENABLE_BLACS ENABLE_BLAS95 ENABLE_CDFT ENABLE_CPARDISO ENABLE_LAPACK95 ENABLE_OMP_OFFLOAD ENABLE_SCALAPACK)
mark_as_advanced(MKL_DIR MKL_ARCH MKL_INCLUDE MKL_INTERFACE_FULL MKL_THREADING MKL_VERSION_H)
mark_as_advanced(mkl_core_file mkl_gf_ilp64_file mkl_intel_thread_file mkl_gf_lp64_file mkl_blas95_lp64_file mkl_lapack95_lp64_file mkl_sequential_file mkl_gnu_thread_file)
set_property(CACHE MKL_LINK PROPERTY STRINGS ${MKL_LINK_LIST})
mark_as_advanced(MKL_LINK)
set_property(CACHE MKL_THREADING PROPERTY STRINGS ${MKL_THREADING_LIST})
mark_as_advanced(MKL_THREADING)

# set CFS variables to the values from MKLConfig.cmake
set(MKL_INCLUDE_DIR ${MKL_INCLUDE}) # TODO: homogenise to use original name from MKLConfig once legacy-mkl-finding code is removed

# MKL_LIB_DIR is used by some cfsdeps, find it based on the imported target
get_target_property(mkl_core_lib MKL::mkl_core IMPORTED_LOCATION)
get_filename_component(MKL_LIB_DIR "${mkl_core_lib}" DIRECTORY)

# some (extensive) debug output - TODO: remove once this is considered stable
dump_variables("MKL_")
#cmake_print_properties(TARGETS ${MKL_IMPORTED_TARGETS} PROPERTIES BINARY_DIR IMPORTED IMPORTED_LOCATION SOURCE_DIR INTERFACE_LINK_LIBRARIES)
cmake_print_variables(BLAS_LIBRARY MKL_LIBS MKL_LIB_DIR)
cmake_print_variables(MKL_LINK_LINE)
cmake_print_variables(MKL_THREAD_LIB)
cmake_print_variables(MKL_SUPP_LINK)

# if the threading lib is dynamic (e.g. Intel OMP under Linux) we need to install it so we have a portable installation
if(UNIX)
  set(INSTALL_DESTINATION "lib")
else()
  set(INSTALL_DESTINATION "bin")
endif()
# threading lib
if("${MKL_THREAD_LIB}" MATCHES "${CMAKE_SHARED_LIBRARY_SUFFIX}")
  message(STATUS "dynamic linking for MKL_THREAD_LIB = ${MKL_THREAD_LIB}")
  # this is needed for distributable builds
  install(PROGRAMS ${MKL_THREAD_LIB} DESTINATION "${INSTALL_DESTINATION}")
  message(STATUS "  will install ${MKL_THREAD_LIB} to '${INSTALL_DESTINATION}' upon installation (cpack)")
endif()
# now for each target
message(STATUS "searching for dynamically linked libraries to install ...")
foreach(lib ${MKL_LIBRARIES})
  set(tgt MKL::${lib})
  get_target_property(il ${tgt} IMPORTED_LOCATION)
  if("${il}" MATCHES "${CMAKE_SHARED_LIBRARY_SUFFIX}")
    message(STATUS "  will install ${il} to '${INSTALL_DESTINATION}' upon installation (cpack)")
    install(IMPORTED_RUNTIME_ARTIFACTS ${tgt} DESTINATION "${INSTALL_DESTINATION}" NAMELINK_COMPONENT)
    get_filename_component(il_real "${il}" REALPATH) # modify to show real path
    if(NOT "${il}" STREQUAL "${il_real}") # must be a symlink
      message(STATUS "  will install ${il_real} to '${INSTALL_DESTINATION}' upon installation (cpack)")
      install(PROGRAMS ${il_real} DESTINATION "${INSTALL_DESTINATION}")  
    endif()
  else()
    message(STATUS "  it seems ${tgt}.IMPORTED_LOCATION=${il} does not need to be installed")
  endif()
endforeach()
# now copy over missing dlls (probably only needed on Windows)
if("${MKL_LINK}" MATCHES "dynamic")
  message(STATUS "searching for necessary dynamic libs in ${MKL_LIB_DIR}")
  foreach(libname "mkl_avx" "mkl_def" "mkl_vml")
    set(globname "${CMAKE_SHARED_LIBRARY_PREFIX}${libname}*${CMAKE_SHARED_LIBRARY_SUFFIX}")
    if(UNIX)
      set(globname "${globname}.*")# numbered versions
    endif()
    file(GLOB libs "${MKL_LIB_DIR}/${globname}")
    message(STATUS "  '${globname}': copy ${libs} to '${INSTALL_DESTINATION}' upon installation (cpack)") 
    install(PROGRAMS ${libs} DESTINATION "${INSTALL_DESTINATION}")
  endforeach()
endif()

# some final checks
if(USE_OPENMP AND NOT EXISTS "${MKL_THREAD_LIB}")
  message(WARNING "MKL_THREAD_LIB=${MKL_THREAD_LIB} from MKLConfig (above) does not exist!")
endif()

set(MKL_BLAS_LIB "$<LINK_ONLY:MKL::MKL>") # used in TARGET_LL of all openCFS targets

else(MKL_FOUND) # legacy bahaviour below TODO: remove once we do not require old systems in the pipeline any more
#if(NOT MKL_FOUND)
message(WARNING "MKL was not found via MKLConfig.cmake! 
Either your MKL installation is too old (< 2021.3) or the possible paths (MKL_POSSIBLE_PATHS) above need to be refined.
please install a recent MKL >= 2021.3 and/or set the correct hints in MKL_POSSIBLE_PATHS to find MKLConfig.cmake
Trying the old way instead - not guaraneed to work and not maintained any more ...")
#-------------------------------------------------------------------------------
# This script  is responsible for the  determination of the correct  include and
# linker  parameters for  the Intel  Math Kernel  Library (MKL).  We distinguish
# between Windows, MacOS and Linux.
#
# On Linux  we determine the  required linker flags by compiling one of the MKL
# examples and reading the flags from the Makefile output (cf. below).
# another option: https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
#-------------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Function to read MKL version from header file mkl.h or mkl_version.h
#-----------------------------------------------------------------------------
function(MKL_VERSION_FROM_HEADER)
  if(EXISTS "${MKL_ROOT_DIR}/include/mkl_version.h")
    file(STRINGS "${MKL_ROOT_DIR}/include/mkl_version.h" MKL_HEADER)
  elseif(EXISTS "${MKL_ROOT_DIR}/include/mkl.h")
    file(STRINGS "${MKL_ROOT_DIR}/include/mkl.h" MKL_HEADER)
  else()
    message(FATAL_ERROR "MKL header could not be found for MKL_ROOT_DIR=${MKL_ROOT_DIR} at include/mkl_version.h or include/mkl.h")
  endif()

  foreach(line IN LISTS MKL_HEADER)
    IF(line MATCHES "#define __INTEL_MKL")
      STRING(REPLACE "_" ";" MKL_LIST "${line}")
      LIST(GET MKL_LIST 4 ITEM_FOUR)

      IF(ITEM_FOUR STREQUAL "")
        LIST(GET MKL_LIST 5 MKL_MAJOR_VERSION)
        STRING(STRIP "${MKL_MAJOR_VERSION}" MKL_MAJOR_VERSION)
      ELSEIF(ITEM_FOUR STREQUAL "MINOR")
        LIST(GET MKL_LIST 6 MKL_MINOR_VERSION)
        STRING(STRIP "${MKL_MINOR_VERSION}" MKL_MINOR_VERSION)
      ELSEIF(ITEM_FOUR STREQUAL "UPDATE")
        LIST(GET MKL_LIST 6 MKL_UPDATE)
        STRING(STRIP "${MKL_UPDATE}" MKL_UPDATE)
      ENDIF()
    ENDIF()
  endforeach()

  SET(MKL_MAJOR_VERSION ${MKL_MAJOR_VERSION} PARENT_SCOPE)
  SET(MKL_MINOR_VERSION ${MKL_MINOR_VERSION} PARENT_SCOPE)
  SET(MKL_UPDATE ${MKL_UPDATE} PARENT_SCOPE)
  IF(${MKL_MAJOR_VERSION} LESS 2019)
    MESSAGE(FATAL_ERROR "MKL version has to be at least 2019.0.0 MKL_VERSION: ${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}.${MKL_UPDATE}")
  ENDIF()
endfunction(MKL_VERSION_FROM_HEADER)

if(MSVC) # this is also icx in Windows!
  #-----------------------------------------------------------------------------
  # If not specified by the user, try to determine proper MKL root directory.
  #-----------------------------------------------------------------------------
  if(NOT MKL_ROOT_DIR)
    set(MKL_POSSIBLE_ROOT_DIRS
      "e:/dev/intel/MKL/composer_xe_2013"
      "e:/dev/intel/MKL/10.0.5.025"
      "$ENV{MKLROOT}"
      "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2019.4.228/windows/mkl"
      "C:/Program Files (x86)/Intel/oneAPI/mkl/latest" )

    find_file(MKL_H
      "include/mkl.h"
      PATHS ${MKL_POSSIBLE_ROOT_DIRS}
      DOC "MKL root dir"
      NO_DEFAULT_PATH
      NO_CMAKE_ENVIRONMENT_PATH
      NO_CMAKE_PATH
      NO_SYSTEM_ENVIRONMENT_PATH
      NO_CMAKE_SYSTEM_PATH
      NO_CMAKE_FIND_ROOT_PATH )

    STRING(REPLACE "/include/mkl.h" "" MKL_ROOT_DIR "${MKL_H}")
    #  MESSAGE(FATAL_ERROR "MKL_ROOT_DIR ${MKL_ROOT_DIR}")
  endif()
 
  # e.g. MKL_ROOT_DIR=C:/Program Files (x86)/Intel/oneAPI/mkl/latest
  assert_dir_exists(${MKL_ROOT_DIR})

  #---------------------------------------------------------------------------
  # Read mkl version from header files
  #---------------------------------------------------------------------------
  MKL_VERSION_FROM_HEADER()
  set(MKL_INCLUDE_DIR "${MKL_ROOT_DIR}/include" CACHE PATH "mkl include dir")
  
  if(MKL_MAJOR_VERSION VERSION_GREATER_EQUAL 2024)
    set(MKL_LIB_DIR "${MKL_ROOT_DIR}/lib" CACHE PATH "here we assume the mkl libs")
  else()
    # this should depend on the MKL installation, not on CXX compiler
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL 19.1)
      set(MKL_LIB_DIR "${MKL_ROOT_DIR}/lib/intel64_win" CACHE PATH "here we assume the mkl libs")
    else()
     set(MKL_LIB_DIR "${MKL_ROOT_DIR}/lib/intel64" CACHE PATH "here we assume the mkl libs")
    endif()
  endif()

  assert_dir_exists(${MKL_LIB_DIR})
  set(MKL_ROOT_DIR ${MKL_ROOT_DIR} CACHE PATH "Directory of MKL.")

  # mkl_solver* libs are deprecated as of v10
  # see: https://web.archive.org/web/20091027212420/https://software.intel.com/en-us/articles/mkl_solver_libraries_are_deprecated_libraries_since_version_10_2_Update_2/

  set(MKL_BLAS_LIB
    ${MKL_LIB_DIR}/mkl_intel_lp64_dll.lib
    ${MKL_LIB_DIR}/mkl_intel_thread_dll.lib
    ${MKL_LIB_DIR}/mkl_core_dll.lib )
	
  if(USE_OPENMP)	
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL 19.1)
      set(LIBOMP ${MKL_ROOT_DIR}/../compiler/lib/intel64_win/libiomp5md.lib) 
      assert_file_exists(${LIBOMP})
    else() 
      # there seems to be bo static openmp-lib ?!
      # set(LIBOMP "${MKL_ROOT_DIR}/../../compiler/latest/windows/redist/intel64_win/compiler/libiomp5md.dll")
    endif()
    list(APPEND MKL_BLAS_LIB ${LIBOMP})
  endif() 	

  set(MKL_LAPACK_LIB ${MKL_BLAS_LIB})

  set(DEPS_SEQUENTIAL
    ${MKL_LIB_DIR}/mkl_intel_lp64.lib
    ${MKL_LIB_DIR}/mkl_intel_thread.lib
    ${MKL_LIB_DIR}/mkl_core.lib
    ${MKL_LIB_DIR}/mkl_sequential.lib
    ${LIBOMP} ) # LIBOMP is empty without openMP and with too new MSVC

  # Copy MKL redistributable dlls to bin/ directory
  
  # redistributable directory
  if(MKL_MAJOR_VERSION VERSION_LESS 2024)
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL 19.1)
      set(MKL_REDIST_DIR "${MKL_ROOT_DIR}/../redist/intel64/mkl/")
    else()
      # C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\redist\intel64_win\compiler>
      set(MKL_REDIST_DIR "${MKL_ROOT_DIR}/redist/intel64/")
    endif()
  else()
    set(MKL_REDIST_DIR "${MKL_ROOT_DIR}/lib/")
  endif()
  assert_dir_exists(${MKL_REDIST_DIR})

  set(LIB_DEST_DIR "${CFS_BINARY_DIR}/bin/")
  # copy all dlls to binary dir
  message(STATUS "Copying MKL redistributable files from ${MKL_REDIST_DIR} to ${LIB_DEST_DIR}")
  file(COPY ${MKL_REDIST_DIR} DESTINATION ${LIB_DEST_DIR} FILES_MATCHING PATTERN "*.dll")

# END MSVC
elseif(APPLE) # note tha APPLE is ALSO UNIX! 
  # for APPLE we do it as with MSVC, we forget about the complex LINUX compile_mkl_test.sh.in stuff
  # this is possibly also an option for Linux to get rid of the ugly stuff. 
  # base is the interactive configurator: https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html
  
  if(NOT MKL_ROOT_DIR)
    set(MKL_POSSIBLE_ROOT_DIRS 
       "/opt/intel/mkl"
       $ENV{MKLROOT}) # set by compilervars.sh intel64 

    find_file(MKL_H
      "include/mkl.h"
      PATHS ${MKL_POSSIBLE_ROOT_DIRS}
      DOC "MKL root dir"
      NO_DEFAULT_PATH
      NO_CMAKE_ENVIRONMENT_PATH
      NO_CMAKE_PATH
      NO_SYSTEM_ENVIRONMENT_PATH
      NO_CMAKE_SYSTEM_PATH
      NO_CMAKE_FIND_ROOT_PATH)

    string(REPLACE "/include/mkl.h" "" MKL_ROOT_DIR "${MKL_H}")
    #message(STATUS "MKL_ROOT_DIR set from include/mkl.h: ${MKL_ROOT_DIR}")
  endif()


  #---------------------------------------------------------------------------
  # Read mkl version from header files
  #---------------------------------------------------------------------------
  MKL_VERSION_FROM_HEADER()
  set(MKL_INCLUDE_DIR "${MKL_ROOT_DIR}/include" CACHE PATH "mkl include directory")
  
  set(MKL_LIB_DIR "${MKL_ROOT_DIR}/lib" CACHE PATH "here we assume the mkl libs")

  # https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html
  if(USE_OPENMP)
    set(MKL_THREAD_LIB ${MKL_LIB_DIR}/libmkl_intel_thread.a)
  else()
    set(MKL_THREAD_LIB ${MKL_LIB_DIR}/libmkl_sequential.a)
  endif()
  
  # lp64 is 32bit integer, ilp64 is 64bit integer 

  # -liomp5 is only required in the USE_OPENMP case but it should not harm
  set(MKL_BLAS_LIB 
    ${MKL_LIB_DIR}/libmkl_intel_lp64.a
    ${MKL_THREAD_LIB}
    ${MKL_LIB_DIR}/libmkl_core.a
    -liomp5 
    -lpthread 
    -lm 
    -ldl)

  # for macOS and clang and MKL we need to link -lgcc_s.1.dylib from the libgfortran.a location
  # to prevent undefined references ___divtf3, ___eqtf2, ___gttf2, ..., ___subtf3, ___unordtf2
  # according to nm there is no static lib to provide this (used by libgfortran.a)
  # Edit: at least on macOS there is a static libgcc.a, see GFORTAN_LIBGCC_STATIC in CheckFortranRuntime.cmake
  if(NOT GFORTRAN_LIBRARY)
    message(FATAL_ERROR "neeed path to libgcc_s.1.dylib from GFORTRAN_LIBRARY")
  endif()
  
  get_filename_component(GFORTRAN_LIB_DIR ${GFORTRAN_LIBRARY} DIRECTORY)
  set(MKL_BLAS_LIB
    ${MKL_BLAS_LIB}
    -L${GFORTRAN_LIB_DIR}
    -lgcc_s.1)
     
  set(MKL_LAPACK_LIB ${MKL_BLAS_LIB})

  set(MKL_ROOT_DIR ${MKL_ROOT_DIR} CACHE PATH "Directory of MKL." FORCE)


  mark_as_advanced(MKL_H)
  mark_as_advanced(MKL_ROOT_DIR)

elseif(UNIX AND NOT APPLE) # neither MSVC and neither APPLE. Hence UNIX and Linux. Note that APPLE is also UNIX
  if(NOT MKL_ROOT_DIR)
    set(MKL_POSSIBLE_ROOT_DIRS
       "/opt/intel/mkl"
       "/opt/intel/oneapi/mkl/latest"
       "/share/programs/intel-mkl/latest/mkl"
       $ENV{MKLROOT}) # set by compilervars.sh intel64

    find_file(MKL_H
      "include/mkl.h"
      PATHS ${MKL_POSSIBLE_ROOT_DIRS}
      DOC "MKL root dir"
      NO_DEFAULT_PATH
      NO_CMAKE_ENVIRONMENT_PATH
      NO_CMAKE_PATH
      NO_SYSTEM_ENVIRONMENT_PATH
      NO_CMAKE_SYSTEM_PATH
      NO_CMAKE_FIND_ROOT_PATH)
    mark_as_advanced(MKL_H)

    string(REPLACE "/include/mkl.h" "" MKL_ROOT_DIR "${MKL_H}")
    message(STATUS "MKL_ROOT_DIR set from include/mkl.h: ${MKL_ROOT_DIR}")
  endif()
  set(MKL_ROOT_DIR ${MKL_ROOT_DIR} CACHE PATH "Directory of MKL.")
  mark_as_advanced(MKL_ROOT_DIR)
  
  if(NOT EXISTS ${MKL_ROOT_DIR}/lib/intel64/libmkl_core.a)
    message(FATAL_ERROR "MKL_ROOT_DIR=${MKL_ROOT_DIR} but ${MKL_ROOT_DIR}/lib/intel64/libmkl_core.a is not found") 
  endif()
  set(MKL_LIB_DIR "${MKL_ROOT_DIR}/lib/intel64" CACHE PATH "here we assume the mkl libs")
  message(STATUS "setting MKL_LIB_DIR=${MKL_LIB_DIR} (validated)")

  #---------------------------------------------------------------------------
  # Read mkl version from header files
  #---------------------------------------------------------------------------
  MKL_VERSION_FROM_HEADER()

  # set include dir
  set(MKL_INCLUDE_DIR "${MKL_ROOT_DIR}/include" CACHE PATH "here we assume the mkl includes")

  if(USE_OPENMP)
    set(MKL_OMP_LIB "${MKL_LIB_DIR}/libiomp5.so")
    # strange that we search here for alternative locations of libiomp5.so but not for the other libs later ...
    if(NOT EXISTS "${MKL_OMP_LIB}")
      set(MKL_OMP_LIB "${MKL_ROOT_DIR}/../compiler/lib/intel64_lin/libiomp5.so")
    endif()
    # path for oneAPI 2020.2
    if(NOT EXISTS "${MKL_OMP_LIB}")
      set(MKL_OMP_LIB "${MKL_ROOT_DIR}/../../compiler/latest/linux/compiler/lib/intel64_lin/libiomp5.so")
    endif()
    # path for oneAPI 2024.0: There is a new "Unified Directory Layout", https://www.intel.com/content/www/us/en/developer/articles/release-notes/onemkl-release-notes.html
    if(NOT EXISTS "${MKL_OMP_LIB}")
      set(MKL_OMP_LIB "${MKL_ROOT_DIR}/../../compiler/${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}/lib/libiomp5.so")
    endif()
    cmake_print_variables(MKL_OMP_LIB)
    # copy over
    file(COPY ${MKL_OMP_LIB} DESTINATION ${LIBRARY_OUTPUT_PATH})
    set(MKL_OMP_LIB_LINE "-L${LIBRARY_OUTPUT_PATH} -liomp5")
    set(MKL_THREADING_LIB "${MKL_LIB_DIR}/libmkl_intel_thread.a")
  else()
    set(MKL_THREADING_LIB "${MKL_LIB_DIR}/libmkl_sequential.a")
  endif()

  # set the link line, essentailly from
  # see https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html
  # regarting the correct Fortran interface library (libmkl_gf_lp64):
  # see: https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/560573
  # Select interface layer: Fortran API with 32 bit integer (lp64, not ilp64)
  set(MKL_FORTRAN_INTERFACE_LIB "${MKL_LIB_DIR}/libmkl_gf_lp64.a")
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    message(STATUS "Using Intel Fortran compiler (ifort or ifx): use libmkl_intel_lp64.a")
    set(MKL_FORTRAN_INTERFACE_LIB "${MKL_LIB_DIR}/libmkl_intel_lp64.a")
  endif()
  set(MKL_BLAS_LIB
     ${MKL_LIB_DIR}/libmkl_blas95_lp64.a
     ${MKL_LIB_DIR}/libmkl_lapack95_lp64.a 
     -Wl,--start-group
     ${MKL_FORTRAN_INTERFACE_LIB}
     ${MKL_THREADING_LIB}
     ${MKL_LIB_DIR}/libmkl_core.a
     -Wl,--end-group
     ${MKL_OMP_LIB_LINE}
     -lpthread
     -lm
     -ldl)
  set(MKL_LAPACK_LIB ${MKL_BLAS_LIB})
  set(MKL_LIBS "${MKL_FORTRAN_INTERFACE_LIB},${MKL_THREADING_LIB},${MKL_LIB_DIR}/libmkl_core.a")
else() # end UNIX
  message(FATAL_ERROR "unhandled system type")
endif()  


#-------------------------------------------------------------------------------
# Set BLAS, LAPACK and PARDISO libraries depending on the MKL version.
#-------------------------------------------------------------------------------

# see also External_OpenBLAS and External_LAPACK, where the setting is quite different:
# e.g. BLAS_LIBRARY=-Wl,--start-group;/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl/lib/intel64/libmkl_intel_lp64.a;/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl/lib/intel64/libmkl_gnu_thread.a;/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl/lib/intel64/libmkl_core.a;-Wl,--end-group;-L/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl/../compiler/lib/intel64;-liomp5;-lpthread;-lm;-ldl
# e.g. LAPACK_LIBRARY=
if(USE_BLAS_LAPACK STREQUAL "MKL")
  set(BLAS_LIBRARY "${MKL_BLAS_LIB}")
  if(MKL_MAJOR_VERSION LESS 10)
    set(LAPACK_LIBRARY "${MKL_LAPACK_LIB}")
  endif(MKL_MAJOR_VERSION LESS 10)  
endif()
set(PARDISO_LIBRARY "${MKL_PARDISO_LIB}")

#-------------------------------------------------------------------------------
# Status message of found MKL
#-------------------------------------------------------------------------------
if(DEFINED MKL_UPDATE)
  set(MKL_VERSION "${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}.${MKL_UPDATE}")
else()
  set(MKL_VERSION "${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}")
endif()

# some debug output
# cmake_print_variables(CMAKE_LINK_GROUP_USING_RESCAN)
# cmake_print_variables(CMAKE_CXX_LINK_GROUP_USING_RESCAN_SUPPORTED)
endif(MKL_FOUND)

message(STATUS "defining MKL link-line via MKL_BLAS_LIB=${MKL_BLAS_LIB}")# VERBOSE-TRACE are only supported from cmake 3.15
mark_as_advanced(MKL_LIB_DIR)
mark_as_advanced(MKL_INCLUDE_DIR)
message(STATUS "Using Intel MKL version ${MKL_VERSION}")