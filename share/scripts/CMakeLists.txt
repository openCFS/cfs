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
install(PROGRAMS ${COMMON_SCRIPT} DESTINATION bin)

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
install(PROGRAMS ${CFS_SCRIPT} DESTINATION bin)


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
  install(PROGRAMS ${CFS_SCRIPT} DESTINATION bin)
ENDIF(WIN32)

#-------------------------------------------------------------------------------
# Custom target for copying cfstool to ${CFS_BINARY_DIR}/bin
#-------------------------------------------------------------------------------
IF(BUILD_CFSTOOL)
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
  install(PROGRAMS ${CFSTOOL_SCRIPT} DESTINATION bin)

ENDIF(BUILD_CFSTOOL)

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
  install(PROGRAMS ${CFSDAT_SCRIPT} DESTINATION bin)

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
  install(PROGRAMS ${H5TOOL_SCRIPT} DESTINATION bin)

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
install(PROGRAMS ${DISTRO_SCRIPT} DESTINATION bin)

# THIS IS A PROTOTYPE FOR NSIS INSTALLER --> not yet included/supported
#----------------------------------------------------------------------------------
# Custom target for configuring NSIS script and copying it to CMAKE_INSTALL_PREFIX
#----------------------------------------------------------------------------------
#IF(WIN32)
#  IF(${CFS_VERSION_SUFFIX} MATCHES "SR[1-9]")
#    STRING(REPLACE "SR" "" SR_NUM ${NACS_VERSION_SUFFIX})
#  ELSE(${CFS_VERSION_SUFFIX} MATCHES "SR[1-9]")
#    SET(SR_NUM 0)
#  ENDIF(${CFS_VERSION_SUFFIX} MATCHES "SR[1-9]")
#  SET(CFS_FULL_VERSION "${CFS_VERSION_YEAR}.${CFS_VERSION_MONTH}.0.0")
#  CONFIGURE_FILE(
#    ${CMAKE_CURRENT_SOURCE_DIR}/nsis/CFS-setup.nsi.in
#    ${CFS_BINARY_DIR}/CFS-win${SYSTEM_BITS}-setup.nsi
#    @ONLY)
#  INSTALL(FILES ${CFS_BINARY_DIR}/CFS-win${SYSTEM_BITS}-setup.nsi
#                ${CMAKE_CURRENT_SOURCE_DIR}/nsis/defines.nsh
#                ${CMAKE_CURRENT_SOURCE_DIR}/nsis/EnvVarUpdate.nsh
#		${CMAKE_CURRENT_SOURCE_DIR}/nsis/CFS.ico
#		${CMAKE_CURRENT_SOURCE_DIR}/nsis/vcredist_2017_en_x64.exe
#          DESTINATION .
#         )
#ENDIF(WIN32)
