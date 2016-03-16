SET(DESIGN_SRCS
  DesignSpace.cc
  DesignElement.cc
  DesignStructure.cc
  DensityFile.cc
  AuxDesign.cc
  ShapeDesign.cc
  DesignMaterial.cc
  Filter.cc
)

ADD_LIBRARY(design STATIC ${DESIGN_SRCS})

SET(TARGET_LL
  cfsgeneral
  datainout
  domain
  logging
  matvec
  optimization
  optimizer
  paramh
  pde
  utils
)

IF(USE_SGPP)
  SET(TARGET_LL ${TARGET_LL} ${SGPP_LIBRARY})
ENDIF(USE_SGPP)

TARGET_LINK_LIBRARIES(design ${TARGET_LL})
