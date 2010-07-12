SET(DESIGN_SRCS DesignElement.cc DesignMaterial.cc DesignSpace.cc DesignStructure.cc
                Filter.cc ShapeDesign.cc DensityFile.cc)

ADD_LIBRARY(design STATIC ${DESIGN_SRCS})

SET(TARGET_LL
  domain
  ${TARGET_LL}
)

TARGET_LINK_LIBRARIES(design ${TARGET_LL})
