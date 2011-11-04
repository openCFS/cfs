IF (LATEX_COMPILER)
  
  IF(WIN32)
    SET(GEN_DOCU_IN "${CMAKE_CURRENT_SOURCE_DIR}/makedoc.bat.in")
    SET(GEN_DOCU "${CMAKE_CURRENT_BINARY_DIR}/makedoc.bat")
  ELSE(WIN32)
    SET(GEN_DOCU_IN "${CMAKE_CURRENT_SOURCE_DIR}/makedoc.sh.in")
    SET(GEN_DOCU "${CMAKE_CURRENT_BINARY_DIR}/makedoc.sh")
  ENDIF(WIN32)

  CONFIGURE_FILE(${GEN_DOCU_IN}
    ${GEN_DOCU}
    @ONLY )
  MESSAGE(STATUS "configured ${GEN_DOCU_IN} --> ${GEN_DOCU}")

  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/../auxilliaries/pyg.sh.in
    ${CMAKE_CURRENT_BINARY_DIR}/../auxilliaries/pyg.sh
    @ONLY )

  ADD_CUSTOM_TARGET(doc-user
                    ${GEN_DOCU} "build"
		    )
                    
  ADD_CUSTOM_TARGET(doc-user-html
                    ${GEN_DOCU} "build-html"
		    )

  ADD_CUSTOM_TARGET(doc-user-rebuild
                    ${GEN_DOCU} "rebuild"
		    )

  ADD_CUSTOM_TARGET(doc-user-clean
                    ${GEN_DOCU} "clean"
		    )


ELSE (LATEX_COMPILER)
  MESSAGE(STATUS "LaTeX not found")
ENDIF(LATEX_COMPILER)
