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
  SplineBoxDesign.cc)
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
  # see FindPrograms.cmake on how the PYTHON_ stuff is set
  set(TARGET_LL ${TARGET_LL} ${PYTHON_LIBRARY})
  # already in FindCFSDEPS.cmake: include_directories(${PYTHON_INCLUDE_DIR})
  include_directories(${PYTHON_SITE_PACKAGES_DIR})
  set(DESIGN_SRCS ${DESIGN_SRCS} SpaghettiDesign.cc)
endif()

add_library(design STATIC ${DESIGN_SRCS})

if(USE_SGPP)
  set(TARGET_LL ${TARGET_LL} ${SGPP_LIBRARY})
endif()

target_link_libraries(design ${TARGET_LL})
