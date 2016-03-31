SET(DESIGN_SRCS
  DesignSpace.cc
  DesignElement.cc
  DesignStructure.cc
  DensityFile.cc
  AuxDesign.cc
  ShapeMapDesign.cc
  ShapeDesign.cc
  DesignMaterial.cc
  Filter.cc
)

ADD_LIBRARY(design STATIC ${DESIGN_SRCS})

SET(TARGET_LL
  design
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

  # Add a switch for the compiler.
  # As the SGpp toolbox is not available for everyone, we have to exclude the corresponding source code from compiling.
  ADD_DEFINITIONS(-DUSE_SGPP)
ENDIF(USE_SGPP)

TARGET_LINK_LIBRARIES(design ${TARGET_LL})
