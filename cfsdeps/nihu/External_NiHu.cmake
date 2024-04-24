include(ExternalProject)

# Builds external third party projects.
# The parent script CMakeLists.txt defines the "GLOBAL_OUTPUT_PATH" variable,
# which will be used as output directory for all *.lib, *.dll, *.a, *.so, *.pdb files.

###############################################################################
# NiHu
###############################################################################

###link_directories(${GLOBAL_OUTPUT_PATH}) # Where else to put!?
###file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/external/NIHU/)

set(NIHU_SRC_DIR ${CMAKE_BINARY_DIR}/cfsdeps/nihu/src/NiHu/src)

include(ExternalProject)

ExternalProject_Add(
  NiHu

  PREFIX ${CMAKE_BINARY_DIR}/cfsdeps/nihu    # external project's base path
  GIT_REPOSITORY "git://last.hit.bme.hu/toolbox/nihu.git"
  GIT_TAG "nightly"

  UPDATE_COMMAND "" # specify if needed
  PATCH_COMMAND ""  # specify if needed

  ## path to NiHu's CMakeLists
  CMAKE_ARGS -H${CMAKE_BINARY_DIR}/cfsdeps/nihu/src
  CMAKE_ARGS -H${NIHU_SRC_DIR}
             -DNIHU_INSTALL_DIR=${CMAKE_BINARY_DIR}/cfsdeps/nihu/src/NiHu-install
            #  -DNIHU_EIGEN_INSTALL=1
             -DNIHU_EIGEN_PATH=${CMAKE_BINARY_DIR}/cfsdeps/eigen/install/include/Eigen
             -DNIHU_MATLAB_PATH=/usr/local/MATLAB/R2022b
)