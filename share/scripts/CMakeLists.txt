
# Custom target for copying distro.sh to ${CFS_BINARY_DIR}/share/scripts
if(UNIX)
  SET(DISTRO_SCRIPT "distro.sh")
  add_custom_target(distro_script ALL
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/${DISTRO_SCRIPT} ${CFS_BINARY_DIR}/share/scripts
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${DISTRO_SCRIPT}
    COMMENT "Copying \"${DISTRO_SCRIPT}\" script to ${CFS_BINARY_DIR}/share/scripts folder..." )
  install(PROGRAMS ${DISTRO_SCRIPT} DESTINATION share/scripts)
endif()
