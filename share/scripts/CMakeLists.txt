CONFIGURE_FILE(hdf5_import_mesh.sh.in ${CFS_BINARY_DIR}/bin/hdf5_import_mesh.sh)
CONFIGURE_FILE(hdf5_import_dataset.sh.in ${CFS_BINARY_DIR}/bin/hdf5_import_dataset.sh)
CONFIGURE_FILE(hdf5_aux.sh.in ${CFS_BINARY_DIR}/bin/hdf5_aux.sh)
CONFIGURE_FILE(hdf5_concat.sh.in ${CFS_BINARY_DIR}/bin/hdf5_concat.sh)

#-------------------------------------------------------------------------------
# Custom target for copying common.sh to ${CFS_BINARY_DIR}/bin
#-------------------------------------------------------------------------------
IF(WIN32)
  SET(COMMON_SCRIPT "common.bat")
ELSE(WIN32)
  SET(COMMON_SCRIPT "common.sh")
ENDIF(WIN32)

ADD_CUSTOM_TARGET(common_script ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${COMMON_SCRIPT} ${CFS_BINARY_DIR}/bin
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${COMMON_SCRIPT}
  COMMENT "Copying \"${COMMON_SCRIPT}\" script to ${CFS_BINARY_DIR}/bin folder...")

IF(MINGW)
  ADD_CUSTOM_TARGET(common_script_sh ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/common.sh ${CFS_BINARY_DIR}/bin
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/common.sh
    COMMENT "Copying \"common.sh\" script to ${CFS_BINARY_DIR}/bin folder...")
ENDIF()

#-------------------------------------------------------------------------------
# Custom target for copying cfs to ${CFS_BINARY_DIR}/bin
#-------------------------------------------------------------------------------
IF(WIN32)
  SET(CFS_SCRIPT "cfs.bat")
ELSE(WIN32)
  SET(CFS_SCRIPT "cfs")
ENDIF(WIN32)

ADD_CUSTOM_TARGET(cfs_script ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${CFS_SCRIPT} ${CFS_BINARY_DIR}/bin
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${CFS_SCRIPT}
  COMMENT "Copying CFS start script \"${CFS_SCRIPT}\" to
  ${CFS_BINARY_DIR}/bin folder...")

IF(MINGW)
  ADD_CUSTOM_TARGET(cfs_script_sh ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/cfs ${CFS_BINARY_DIR}/bin
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cfs
    COMMENT "Copying \"cfs\" script to ${CFS_BINARY_DIR}/bin folder...")
ENDIF()

#-------------------------------------------------------------------------------
# Custom target for copying cfs-debug to ${CFS_BINARY_DIR}/bin (only LINUX)
#-------------------------------------------------------------------------------
IF(WIN32)
  #nothing to be done here
ELSE(WIN32)
  ADD_CUSTOM_TARGET(cfs_debug_script ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/cfs-debug ${CFS_BINARY_DIR}/bin
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cfs-debug
    COMMENT "Copying CFS-debug start script \"${CFS_SCRIPT}\" to
    ${CFS_BINARY_DIR}/bin folder...")
ENDIF(WIN32)

#-------------------------------------------------------------------------------
# Custom target for copying cfstool to ${CFS_BINARY_DIR}/bin
#-------------------------------------------------------------------------------
IF(CFSTOOL)
  IF(WIN32)
    SET(CFSTOOL_SCRIPT "cfstool.bat")
  ELSE(WIN32)
    SET(CFSTOOL_SCRIPT "cfstool")
  ENDIF(WIN32)

  ADD_CUSTOM_TARGET(cfstool_script ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${CFSTOOL_SCRIPT} ${CFS_BINARY_DIR}/bin
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${CFSTOOL_SCRIPT}
    COMMENT "Copying CFS Tool start script \"${CFSTOOL_SCRIPT}\" to
    ${CFS_BINARY_DIR}/bin folder...")

  IF(MINGW)
    ADD_CUSTOM_TARGET(cfstool_script_sh ALL
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/cfstool ${CFS_BINARY_DIR}/bin
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cfstool
      COMMENT "Copying \"cfstool\" script to ${CFS_BINARY_DIR}/bin folder...")
  ENDIF()

ENDIF(CFSTOOL)

#-------------------------------------------------------------------------------
# Custom target for copying cfsdat to ${CFS_BINARY_DIR}/bin
#-------------------------------------------------------------------------------
IF(CFSDAT)
  IF(WIN32)
    SET(CFSDAT_SCRIPT "cfsdat.bat")
  ELSE(WIN32)
    SET(CFSDAT_SCRIPT "cfsdat")
  ENDIF(WIN32)

  ADD_CUSTOM_TARGET(cfsdat_script ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${CFSDAT_SCRIPT} ${CFS_BINARY_DIR}/bin
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${CFSDAT_SCRIPT}
    COMMENT "Copying CFSDat start script \"${CFSDAT_SCRIPT}\" to
    ${CFS_BINARY_DIR}/bin folder...")

  IF(MINGW)
    ADD_CUSTOM_TARGET(cfsdat_script_sh ALL
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/cfsdat ${CFS_BINARY_DIR}/bin
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cfsdat
      COMMENT "Copying \"cfsdat\" script to ${CFS_BINARY_DIR}/bin folder...")
  ENDIF()
ENDIF(CFSDAT)

# we skip to copy UNIT_TEST cfstests here, I usually directly link to the bin directory
# whoever needs it might add it :)


#-------------------------------------------------------------------------------
# Custom target for copying h5tool to ${CFS_BINARY_DIR}/bin
#-------------------------------------------------------------------------------
IF(WIN32)
  SET(H5TOOL_SCRIPT "h5tool.bat")
ELSE(WIN32)
  SET(H5TOOL_SCRIPT "h5tool")
ENDIF(WIN32)

ADD_CUSTOM_TARGET(h5tool_script ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${H5TOOL_SCRIPT} ${CFS_BINARY_DIR}/bin
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${H5TOOL_SCRIPT}
  COMMENT "Copying HDF5 Tool start script \"${H5TOOL_SCRIPT}\" to
  ${CFS_BINARY_DIR}/bin folder...")

#-------------------------------------------------------------------------------
# Custom target for copying distro.sh to ${CFS_BINARY_DIR}/share/scripts
#-------------------------------------------------------------------------------
IF(WIN32)
  SET(DISTRO_SCRIPT "winver.bat")
ELSE(WIN32)
  SET(DISTRO_SCRIPT "distro.sh")
ENDIF(WIN32)

ADD_CUSTOM_TARGET(distro_script ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${DISTRO_SCRIPT} ${CFS_BINARY_DIR}/share/scripts
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${DISTRO_SCRIPT}
  COMMENT "Copying \"${DISTRO_SCRIPT}\" script to ${CFS_BINARY_DIR}/share/scripts folder..."
  )

IF(MINGW)
  ADD_CUSTOM_TARGET(distro_script_sh ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/distro.sh ${CFS_BINARY_DIR}/share/scripts
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/distro.sh
    COMMENT "Copying \"distro.sh\" script to ${CFS_BINARY_DIR}/share/scripts folder...")
ENDIF()
