# this shall collect all cmake install() of openCFS and redistributale system stuff (fortran, libiomp, ..)
# the stuff is only copied when calling cpack -G ZIP (note that most other generators don't work)

# TODO: move stuff from CheckFortranRuntime() here

# if the threading lib is dynamic (Intel OpenMP) we need to install it to we have a portable installation
if(UNIX)
  set(INSTALL_DESTINATION "lib")
else()
  set(INSTALL_DESTINATION "bin")
endif()

# only the Intel threading lib is dynamic and needs to be installed. On Linux with gnu compilers mkl uses gnu
#MKL_THREAD_LIB=/opt/intel/oneapi/compiler/2024.2/lib/libiomp5.so
#MKL_THREAD_LIB=C:/Program Files (x86)/Intel/oneAPI/compiler/latest/lib/libiomp5md.lib
if(MKL_THREAD_LIB) # empty in the gnu openmp case and legacy_mkl case (which installs itself)
  install(PROGRAMS ${MKL_THREAD_LIB} DESTINATION "${INSTALL_DESTINATION}")
  message(STATUS "  will install ${MKL_THREAD_LIB} to '${INSTALL_DESTINATION}' upon installation (cpack)")
endif()


# choose what to install
# see https://cmake.org/cmake/help/v3.9/command/install.html#installing-directories

# tailing / in source leads to rename to destionation
#install(PROGRAMS "${CFS_BINARY_DIR}/bin" DESTINATION "." PATTERN "[A-Z]*" EXCLUDE)
install(DIRECTORY "${EXECUTABLE_OUTPUT_PATH}/" DESTINATION "bin" USE_SOURCE_PERMISSIONS
  PATTERN "*.lib" EXCLUDE
  PATTERN "*.cmake" EXCLUDE
  PATTERN "*.in" EXCLUDE
  PATTERN "*.exp" EXCLUDE
  PATTERN "*/vtkHashSource-*" EXCLUDE
  PATTERN "*.manifest" EXCLUDE)

install(DIRECTORY "${LIBRARY_OUTPUT_PATH}/" DESTINATION "lib" USE_SOURCE_PERMISSIONS
  PATTERN "*.a" EXCLUDE
  PATTERN "*.la" EXCLUDE
  PATTERN "*.pc" EXCLUDE
  PATTERN "*.lib" EXCLUDE
  PATTERN "*/cmake/*" EXCLUDE
  PATTERN "*.in" EXCLUDE
  PATTERN "pkgconfig" EXCLUDE)

# the output 'Copying "distro.sh" script to /Users/fwein/code/cfs/release/share/scripts folder...' 
# ist via add_custom_target("distro.sh") in share/scripts/CMakeLists.txt 

# license install
file(COPY "${CFS_SOURCE_DIR}/LICENSE" DESTINATION "${CFS_BINARY_DIR}/license/")
install(DIRECTORY "${CFS_BINARY_DIR}/license/" DESTINATION "license" USE_SOURCE_PERMISSIONS)

include(CPack)

