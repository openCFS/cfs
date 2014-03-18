CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/doc-common.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/doc-common.cmake"
  @ONLY
)

ADD_SUBDIRECTORY(developer)
ADD_SUBDIRECTORY(user)

ADD_CUSTOM_TARGET(doc
  DEPENDS doc-devel doc-user
  )
