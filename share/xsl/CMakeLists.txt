SET(SRC_DIR ${CFS_SOURCE_DIR}/share/xsl)
SET(DST_DIR ${CFS_BINARY_DIR}/share/xsl)

ADD_CUSTOM_TARGET(share_xsl ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SRC_DIR}/check_for_docu_tag.xsl ${DST_DIR}/check_for_docu_tag.xsl
  COMMENT "Updating XSL style sheets..."
)

CONFIGURE_FILE("${SRC_DIR}/cfs_info_output_html.xsl.in" "${DST_DIR}/cfs_info_output_html.xsl")