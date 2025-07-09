set(DESIGN_SRCS
  AuxDesign.cc
  DensityFile.cc
  DesignElement.cc
  DesignMaterial.cc
  DesignSpace.cc
  DesignStructure.cc
  FeaturedDesign.cc
  Filter.cc
  LocalElementCache.cc
  ShapeDesign.cc
  ShapeMapDesign.cc
  FeatureMappingDesign.cc
  SplineBoxDesign.cc
  MaterialTensor.cc)
  # for SpahettiDesign.cc see below 

set(TARGET_LL
  optimization
  cfsgeneral
  datainout
  domain
  logging
  matvec
  optimizer
  paramh
  pde
  utils)

# Currently SpaggettiDesign.cc works only with Python. 
# The Python part shall be split later as optional part
if(USE_EMBEDDED_PYTHON)
  # see FindPrograms.cmake on how the =ython_ stuff is set
  list(APPEND TARGET_LL ${Python_LIBRARIES})
  set(DESIGN_SRCS ${DESIGN_SRCS} SpaghettiDesign.cc)
endif()

add_library(design STATIC ${DESIGN_SRCS})

if(USE_SGPP)
  set(TARGET_LL ${TARGET_LL} ${SGPP_LIBRARY})
endif()

target_link_libraries(design ${TARGET_LL})
