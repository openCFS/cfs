SET(SRC_DIR ${CFS_SOURCE_DIR}/share/SMSM_directions)
SET(DST_DIR ${CFS_BINARY_DIR}/share/SMSM_directions)

ADD_CUSTOM_TARGET(share_smsm_directions ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/TABSPHEREI4S4_2562.txt ${DST_DIR}/TABSPHEREI4S4_2562.txt
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/TABSPHERE34635_test.txt ${DST_DIR}/TABSPHERE34635_test.txt
  COMMENT "Updating SMSM direction tables..."
)

install(DIRECTORY "${CFS_BINARY_DIR}/share/SMSM_directions" DESTINATION "share/" USE_SOURCE_PERMISSIONS PATTERN "CMakeFiles" EXCLUDE PATTERN "Makefile" EXCLUDE PATTERN "*.dir" EXCLUDE PATTERN "*.cmake" EXCLUDE)
