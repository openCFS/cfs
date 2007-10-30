# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(MATERIALS_SRCS acousticMaterial.cc baseMaterial.cc
  electroMagneticMaterial.cc electrostaticMaterial.cc flowMaterial.cc
  heatMaterial.cc mechanicMaterial.cc piezoMaterial.cc
  thermoelasticMaterial.cc pyroelectricMaterial.cc)
# MaterialData.cc is removed (?)

ADD_LIBRARY(materials STATIC ${MATERIALS_SRCS})

TARGET_LINK_LIBRARIES(materials
  utils
  domain
)
