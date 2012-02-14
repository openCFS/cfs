IF(USE_GMSH)
SET(SIMINOUTGMSH_SRCS
  simInputGmsh.cc
  simOutputGmsh.cc
  GmshHelper.cc
)

  ADD_LIBRARY(siminoutgmsh STATIC ${SIMINOUTGMSH_SRCS})
  
  TARGET_LINK_LIBRARIES(siminoutgmsh
    cfsgeneral
    datainout
    domain
    logging
    matvec
    paramh
    pde
    utils
  )
  
ENDIF(USE_GMSH)

