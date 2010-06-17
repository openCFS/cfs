SET(DESIGN_SRCS DesignElement.cc DesignMaterial.cc DesignSpace.cc DesignStructure.cc
                Filter.cc ShapeDesign.cc)

ADD_LIBRARY(design STATIC ${DESIGN_SRCS})

TARGET_LINK_LIBRARIES(design optimization ${TARGET_LL})
