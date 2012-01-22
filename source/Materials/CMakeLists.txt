# prevent strict aliasing for due to bug in boost signal library v1.35
#IF(CFS_CXX_COMPILER_NAME STREQUAL "GCC")
#set_source_files_properties(mechanicMaterial.cc
#                            PROPERTIES COMPILE_FLAGS
#                            -fno-strict-aliasing)
#ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "GCC")

SET(MATERIALS_SRCS
  BaseMaterial.cc
  ElectroMagneticMaterial.cc 
  ElectroStaticMaterial.cc 
  HeatMaterial.cc 
  AcousticMaterial.cc
  Models/Hysteresis.cc
  Models/Jiles.cc
  Models/PiezoMicroModelBK.cc
  Models/PiezoMicroModelHF.cc
  Models/Preisach.cc
  Models/SimplePreisachInv.cc 
  )


ADD_LIBRARY(materials STATIC ${MATERIALS_SRCS})

TARGET_LINK_LIBRARIES(materials
  utils
  domain
  mesh
  odesolver
)
