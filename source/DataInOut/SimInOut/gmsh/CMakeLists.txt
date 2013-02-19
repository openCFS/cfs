IF(USE_GMSH)
  SET(SIMINOUTGMSH_SRCS SimInputGmsh.cc SimOutputGmsh.cc SimOutputParsed.cc GmshHelper.cc)

  ADD_LIBRARY(siminoutgmsh STATIC ${SIMINOUTGMSH_SRCS})
  
  TARGET_LINK_LIBRARIES(siminoutgmsh febasis)
ENDIF(USE_GMSH)

