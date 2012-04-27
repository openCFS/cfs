IF (LATEX_COMPILER)

  CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/doc-advprepost.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/doc-advprepost.cmake"
    @ONLY )

  ADD_CUSTOM_TARGET(doc-advprepost
                    ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_BINARY_DIR}/doc-advprepost.cmake"
		    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/doc-advprepost.cmake"
		    )

  CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/doc-xmlref.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/doc-xmlref.cmake"
    @ONLY )

  ADD_CUSTOM_TARGET(doc-xmlref
                    ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_BINARY_DIR}/doc-xmlref.cmake"
		    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/doc-xmlref.cmake"
		    )

  CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/doc-xmlmat.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/doc-xmlmat.cmake"
    @ONLY )

  ADD_CUSTOM_TARGET(doc-xmlmat
                    ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_BINARY_DIR}/doc-xmlmat.cmake"
		    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/doc-xmlmat.cmake"
		    )
  
  ADD_CUSTOM_TARGET(doc-user
    DEPENDS doc-advprepost doc-xmlref doc-xmlmat
    )                    

ELSE (LATEX_COMPILER)
  MESSAGE(STATUS "LaTeX not found")
ENDIF(LATEX_COMPILER)
