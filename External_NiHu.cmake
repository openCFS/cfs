include(ExternalProject)

# Builds external third party projects.
# The parent script CMakeLists.txt defines the "GLOBAL_OUTPUT_PATH" variable,
# which will be used as output directory for all *.lib, *.dll, *.a, *.so, *.pdb files.

###############################################################################
# NiHu
###############################################################################

###link_directories(${GLOBAL_OUTPUT_PATH}) # Where else to put!?

###if(USE_NIHU)

  ###file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/external/NIHU/)
  ###file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/external/TEST/)
  ###file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/external/NIHU/src/NiHu/build_dir)

  ##include(CMakeLists_NiHu.txt)

  message("STARTING EXTERNAL_NIHU!!")
  include(ExternalProject)

  ExternalProject_Add(
    NiHu
    ###PREFIX ${CMAKE_BINARY_DIR}/external/NIHU    # external project's base path
    PREFIX ${CMAKE_BINARY_DIR}/cfsdeps/nihu    # external project's base path
    GIT_REPOSITORY "git://last.hit.bme.hu/toolbox/nihu.git"
    GIT_TAG "nightly"

    UPDATE_COMMAND "" # specify if needed
    PATCH_COMMAND ""  # specify if needed

    ## path to NiHu's CMakeLists
    ###CMAKE_ARGS -H${CMAKE_BINARY_DIR}/external/NIHU/src/NiHu/src   # NiHu's CMakeLists location
                                                                  # always one directory deeper
                                                                  # maybe there is a more "general"
                                                                  # less hard-coded way
    CMAKE_ARGS -H${CMAKE_BINARY_DIR}/cfsdeps/nihu/src/NiHu/src
    # Here we could utilize an alias e.g.
    # CMAKE_ARGS -H${NIHU_BINARY_DIR}/src 

    ## path to NiHu's install directory
    # CMAKE_ARGS -DNIHU_INSTALL_DIR="{CMAKE_BINARY_DIR}/external/NIHU/src/NiHu/build_dir/install_dir"
    ###CMAKE_ARGS -DNIHU_INSTALL_DIR="/install_dir"
    ###CMAKE_ARGS -DNIHU_INSTALL_DIR=/install_dir
    CMAKE_ARGS -DNIHU_INSTALL_DIR=${CMAKE_BINARY_DIR}/cfsdeps/nihu/src/NiHu-build/install_dir
    ###CMAKE_ARGS -DNIHU_INSTALL_DIR=${CMAKE_BINARY_DIR}/cfsdeps/nihu/src/NiHu-install # potentially || potentially anywhere

    ## abandone install step for this external project
    # message("STOP EXTERNAL_NIHU INSTALLATION!!")
    # INSTALL_COMMAND cmake -E echo "SKIPPING ALL INSTALL STEPS!!"
  )

###endif()