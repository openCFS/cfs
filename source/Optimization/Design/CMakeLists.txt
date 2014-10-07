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

TARGET_LINK_LIBRARIES(design ${TARGET_LL})
