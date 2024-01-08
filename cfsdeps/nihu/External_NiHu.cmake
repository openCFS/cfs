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
  CMAKE_ARGS -H${CMAKE_BINARY_DIR}/cfsdeps/nihu/src/NiHu/src
  # TODO: make that an alias to further use in NiHu-executable builds

  # Here we could utilize an alias e.g.
  # CMAKE_ARGS -DNIHU_INSTALL_DIR=${CMAKE_BINARY_DIR}/cfsdeps/nihu/src/NiHu-build/install_dir
  # CMAKE_ARGS -H${CMAKE_BINARY_DIR}/cfsdeps/nihu/src/NiHu/src
  #            -DNIHU_INSTALL_DIR=${CMAKE_BINARY_DIR}/cfsdeps/nihu/src/NiHu-build/install_dir
  #            -DNIHU_EIGEN_INSTALL=1
  CMAKE_ARGS -H${NIHU_SRC_DIR}
            #  -DNIHU_INSTALL_DIR=${CMAKE_BINARY_DIR}/cfsdeps/nihu/src/NiHu-build/install_dir
             -DNIHU_INSTALL_DIR=${CMAKE_BINARY_DIR}/cfsdeps/nihu/src/NiHu-install
             -DNIHU_EIGEN_INSTALL=1
)