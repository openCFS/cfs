# finds first OpenMP and then sets blas. 
# to be called before compiler.cmake

if(MSVC)
  # use OpenMP beyond the original outdated 2.0 from MSVC. Note that OpenMP_CXX_VERSION is still 2.0 (bug) also for --version
  # https://stackoverflow.com/questions/71910329/update-openmp-to-higher-versions-in-visual-studio-2022
  # note that we use the intel omp lib and therefore have no the redistributable issue
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.30) # the else is further below
    set(OpenMP_RUNTIME_MSVC "llvm") # input for FindOpenMP
  endif()
endif()

# note that MKL might use libiomp5 even if gnu/msvc is found here
# when Intel compilers are used, here usually libiomp5 is found in first place
find_package(OpenMP) # properties might be used from deps even with USE_OPENMP=OFF

if(USE_OPENMP)
  if(NOT OpenMP_FOUND)
    message(FATAL_ERROR "USE_OPENMP selected but not found in system")
  endif()

  if(MSVC AND CMAKE_VERSION VERSION_LESS 3.30)
    if(${OpenMP_CXX_FLAGS} STREQUAL "-openmp") # cmake uses -openmp instead of /openmp 
      set(OpenMP_CXX_FLAGS "") # reset to prevent warnings that /openmp:llvm will overwrite /openmp  
    endif()
    set(OpenMP_CXX_FLAGS "${OpenMP_CXX_FLAGS} /openmp:llvm")
  endif() 

  set(CFS_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${CFS_CXX_FLAGS}")  

  if(APPLE)
    # homebrew uses since Okt 2022 not the system path and we need to help cfs and lis
    # according to cmake docu, OpenMP_<lang>_INCLUDE_DIR is input and OpenMP_<lang>_INCLUDE_DIRS is output, wherever the OpenMP_<lang>_INCLUDE_DIR is set?!
    # sometimes the C stuff is set and is valid for C++, give it a try
    if(NOT OpenMP_CXX_INCLUDE_DIR)
      # dump_variables("OpenMP")
      set(OpenMP_CXX_INCLUDE_DIR ${OpenMP_C_INCLUDE_DIR})
    endif()
    assert_set(OpenMP_CXX_INCLUDE_DIR) # /opt/homebrew[/opt/libomp]/include
    assert_set(OpenMP_libomp_LIBRARY) # /opt/homebrew[/opt/libomp]/lib/libomp.dylib"
    get_filename_component(OpenMP_LIBDIR ${OpenMP_libomp_LIBRARY} DIRECTORY) # also use for lis
    include_directories(AFTER SYSTEM "${OpenMP_CXX_INCLUDE_DIR}")
    set(CFS_LINKER_FLAGS "${CFS_LINKER_FLAGS} -lomp -L${OpenMP_LIBDIR} ")
  endif()   
endif() # USE_OPENMP

assert_set(CFS_SOURCE_DIR)

# Apple's Accelerate Framework is a system lib - nothing to build
if(USE_BLAS_LAPACK STREQUAL "ACCELERATE")
  set(USE_ACCELERATE 1)
  include("${CFS_SOURCE_DIR}/cmake_modules/FindAppleAccelerate.cmake")
endif()

# load mkl from oneAPI or oneMKL - it will be linked statically
if(USE_BLAS_LAPACK STREQUAL "MKL")
  set(USE_MKL 1)
  include("${CFS_SOURCE_DIR}/cmake_modules/mkl.cmake")
endif()

# NETLIB or OPENBLAS are dependcies and set in FindCFSDEPS.cmake

