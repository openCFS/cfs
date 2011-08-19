IF(USE_GMSH)
  SET(SIMINOUTGMSH_SRCS simInputGmsh.cc simOutputGmsh.cc GmshHelper.cc)

  ADD_LIBRARY(siminoutgmsh STATIC ${SIMINOUTGMSH_SRCS})
ENDIF(USE_GMSH)

