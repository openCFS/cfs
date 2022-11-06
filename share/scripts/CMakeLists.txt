
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
install(PROGRAMS ${DISTRO_SCRIPT} DESTINATION share/scripts)
