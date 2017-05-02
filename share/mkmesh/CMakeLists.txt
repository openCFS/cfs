# configure the files to adapt to installation paths
CONFIGURE_FILE(start.ans.in ${CFS_BINARY_DIR}/share/mkmesh/start.ans)
CONFIGURE_FILE(mkmesh.mlib.in ${CFS_BINARY_DIR}/share/mkmesh/mkmesh.mlib)

#-------------------------------------------------------------------------------
# Custom target for copying mkmesh to ${CFS_BINARY_DIR}/bin
#-------------------------------------------------------------------------------
SET(MKMESH_SCRIPT "mkmesh")

ADD_CUSTOM_TARGET(mkmesh_script ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${MKMESH_SCRIPT} ${CFS_BINARY_DIR}/bin
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${MKMESH_SCRIPT}
  COMMENT "Copying \"${MKMESH_SCRIPT}\" script to ${CFS_BINARY_DIR}/bin folder...")
