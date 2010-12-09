# prevent strict aliasing for due to bug in boost signal library v1.35
IF(CFS_CXX_COMPILER_NAME STREQUAL "GCC")
set_source_files_properties(mechanicMaterial.cc
                            PROPERTIES COMPILE_FLAGS
                            -fno-strict-aliasing)
ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "GCC")

SET(MATERIALS_SRCS acousticMaterial.cc baseMaterial.cc
  electroMagneticMaterial.cc electrostaticMaterial.cc flowMaterial.cc
  heatMaterial.cc mechanicMaterial.cc piezoMaterial.cc
  thermoelasticMaterial.cc pyroelectricMaterial.cc magstrictMaterial.cc)
# MaterialData.cc is removed (?)

ADD_LIBRARY(materials STATIC ${MATERIALS_SRCS})

TARGET_LINK_LIBRARIES(materials
  utils
  domain
)
