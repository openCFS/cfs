IF(USE_GMSH)
  SET(SIMINOUTGMSH_SRCS SimInputGmsh.cc SimOutputGmsh.cc SimOutputParsed.cc GmshHelper.cc)

  IF(CMAKE_COMPILER_IS_GNUCXX) 
    SET_SOURCE_FILES_PROPERTIES(SimInputGmsh.cc SimOutputGmsh.cc GmshHelper.cc
      PROPERTIES COMPILE_FLAGS "-fno-strict-aliasing")
  ENDIF()

  ADD_LIBRARY(siminoutgmsh STATIC ${SIMINOUTGMSH_SRCS})
  
  TARGET_LINK_LIBRARIES(siminoutgmsh febasis)
ENDIF(USE_GMSH)
