SET(DESIGN_SRCS
  AuxDesign.cc
  DensityFile.cc
  DesignElement.cc
  DesignMaterial.cc
  DesignSpace.cc
  DesignStructure.cc
  Filter.cc
  LocalElementCache.cc
  ShapeDesign.cc
  ShapeMapDesign.cc
  SplineBoxDesign.cc
)

ADD_LIBRARY(design STATIC ${DESIGN_SRCS})

SET(TARGET_LL
  optimization
  cfsgeneral
  datainout
  domain
  logging
  matvec
  optimizer
  paramh
  pde
  utils
)

IF(USE_SGPP)
  SET(TARGET_LL ${TARGET_LL} ${SGPP_LIBRARY})
ENDIF(USE_SGPP)

TARGET_LINK_LIBRARIES(design ${TARGET_LL})
