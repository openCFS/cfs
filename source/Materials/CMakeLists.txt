SET(MATERIALS_SRCS
  AcousticMaterial.cc
  BaseMaterial.cc
  ElectroMagneticMaterial.cc 
  ElectroStaticMaterial.cc 
  HeatMaterial.cc
  MechanicMaterial.cc
  FlowMaterial.cc 
  PiezoMaterial.cc
  magStrictMaterial.cc
  ElectricConductionMaterial.cc
  TestMaterial.cc
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
  odesolver
)
