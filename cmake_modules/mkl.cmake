# see https://www.intel.com/content/www/us/en/docs/onemkl/developer-guide-windows/2023-0/cmake-config-for-onemkl.html

# see MKLConfig.cmake from your oneAPI installation /opt/intel/oneapi/mkl/latest/lib/cmake/mkl/MKLConfig.cmake

set(MKL_LINK "static") # when we are static, we do not need to copy/install mkl dlls.
set(MKL_INTERFACE "lp64") # lp64 = 4 bytes integer, ilp64 = 8 bytes integer - the wrong setting makes cfs randomly crash. Default seems to change
if(NOT USE_OPENMP)
  set(MKL_THREADING "sequential" CACHE STRING "The OpenMP implementation we want to use")
elseif(WIN32 OR OpenMP_iomp5_LIBRARY OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  # Intel OpenMP - for Windows or for Intel compilers or with clang on Linux (not AppleClang!) are used
  # OpenMP_iomp5_LIBRARY is e.g. /opt/intel/oneapi/compiler/2021.4.0/linux/compiler/lib/intel64_lin/libiomp5.so
  # Using Clang on Linux with gnu_thread links libomp.so and libgomp.so and segfaults in parallel mode
  set(MKL_THREADING "intel_thread" CACHE STRING "The OpenMP implementation we want to use")
else()
  # standard Linux case with gcc/gfortran - strictly necessary for the musl build, in practice one can also use intel_thread on Linux with gcc
  set(MKL_THREADING "gnu_thread" CACHE STRING "The OpenMP implementation we want to use")
endif()  
set_property(CACHE MKL_THREADING PROPERTY STRINGS gnu_thread intel_thread sequential)

# possible paths to search for MKLConfig.cmake. With oneAPI's setvars, MKLROOT is set. 
set(POSSIBLE_PATHS "$ENV{MKLROOT}" "/opt/intel/oneapi/mkl/latest" "C:/Program Files (x86)/Intel/oneAPI/mkl/latest" "$ENV{INTEL_INSTALL_DIR}/mkl/latest") # last is shared runners in ci pipeline

# searches and runs MKLConfig.cmake which sets most necessary cmake variables "magically"
find_package(MKL CONFIG PATHS ${POSSIBLE_PATHS})   

if(MKL_FOUND)
  # cmake_print_properties(TARGETS MKL::MKL PROPERTIES INTERFACE_COMPILE_OPTIONS INTERFACE_INCLUDE_DIRECTORIES)
  # set CFS variables with the values from MKLConfig.cmake
  set(MKL_INCLUDE_DIR ${MKL_INCLUDE}) # TODO: homogenise to use original name from MKLConfig once legacy_mkl.cmake is removed
  # MKL_LIB_DIR (cfs name) is used by some cfsdeps, find it based on the imported target. 
  get_target_property(mkl_core_lib MKL::mkl_core IMPORTED_LOCATION)
  get_filename_component(MKL_LIB_DIR "${mkl_core_lib}" DIRECTORY)
  # TODO: remove MKL_BLAS_LIB and use BLAS_LIBRARY once we don't need legacy_mkl.cmake any more
  set(MKL_BLAS_LIB "$<LINK_ONLY:MKL::MKL>") # used in TARGET_LL of all openCFS targets
  set(BLAS_LIBRARY "${MKL_BLAS_LIB}")
  cmake_print_variables(MKL_LINK_LINE)
  # generate a resolved MKL-link-line for cfsdeps
  set(MKL_LIBRARIES "${MKL_LINK_LINE}")
  list(FILTER MKL_LIBRARIES INCLUDE REGEX "MKL::") # only keep parts called MKL::
  # now replace all targets with the full paths
  set(MKL_LINK_LINE_CFSDEPS "${MKL_LINK_LINE} ${MKL_THREAD_LIB} ${MKL_SUPP_LINK}")
  foreach(lib ${MKL_LIBRARIES})
    get_target_property(loc ${lib} IMPORTED_LOCATION)
    # message(STATUS "  ${lib} -> ${loc}")
    string(REPLACE "${lib}" "${loc}" MKL_LINK_LINE_CFSDEPS "${MKL_LINK_LINE_CFSDEPS}")
  endforeach()
  string(REPLACE ";" " " MKL_LINK_LINE_CFSDEPS "${MKL_LINK_LINE_CFSDEPS}")
  cmake_print_variables(MKL_LINK_LINE_CFSDEPS)
else()
  # this is the legacy case - probably only for centos6-gcc runner - remove if possible!
  include("${CFS_SOURCE_DIR}/cmake_modules/legacy_mkl.cmake")
endif()


#Linux: MKL_LINK_LINE=$<IF:$<BOOL:OFF>,,>;-Wl,--start-group;MKL::mkl_intel_lp64;MKL::mkl_intel_thread;MKL::mkl_core;-Wl,--end-group
#Linux: MKL_LINK_LINE=$<IF:$<BOOL:OFF>,,>;-Wl,--start-group;MKL::mkl_gf_lp64;MKL::mkl_gnu_thread;MKL::mkl_core;-Wl,--end-group  -lgomp;-lm;-ldl;-lpthread
#Win32: MKL_LINK_LINE=$<IF:$<BOOL:OFF>,,>;MKL::mkl_intel_lp64;MKL::mkl_intel_thread;MKL::mkl_core

# no copy/install of mkl redistributables necessary as we link statically.

mark_as_advanced(ENABLE_BLACS ENABLE_BLAS95 ENABLE_CDFT ENABLE_CPARDISO ENABLE_LAPACK95 ENABLE_OMP_OFFLOAD ENABLE_SCALAPACK ENABLE_TRY_SYCL_COMPILE)
mark_as_advanced(MKL_DIR MKL_ARCH MKL_INCLUDE MKL_INTERFACE_FULL MKL_THREADING MKL_VERSION_H MKL_MPI OMP_LIBRARY OMP_DLL_DIR)
mark_as_advanced(mkl_core_file mkl_cdft_core_file mkl_gf_ilp64_file mkl_intel_thread_file mkl_gf_lp64_file mkl_blas95_lp64_file mkl_blacs_intelmpi_lp64_file mkl_lapack95_lp64_file mkl_sequential_file mkl_gnu_thread_file mkl_scalapack_lp64_file)
