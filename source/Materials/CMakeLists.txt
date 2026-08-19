SET(MATERIALS_SRCS
  AcousticMaterial.cc
  BaseMaterial.cc
  ElectroMagneticMaterial.cc 
  ElectroStaticMaterial.cc 
  ElecQuasistaticMaterial.cc
  HeatMaterial.cc
  MechanicMaterial.cc
  SmoothMaterial.cc
  FlowMaterial.cc 
  PiezoMaterial.cc
  MagStrictMaterial.cc
  ElectricConductionMaterial.cc
  TestMaterial.cc
  Models/Hysteresis.cc
  Models/EBHysteresis.cc
  Models/invEBHysteresis.cc
  Models/SMSM.cc
  Models/Jiles.cc
  Models/Model.cc
  Models/PiezoMicroModelBK.cc
  Models/PiezoMicroModelHF.cc
  Models/VectorPreisachSutor.cc
  Models/Preisach.cc
  Models/VectorPreisachMayergoyz.cc
  Models/SimplePreisachInv.cc 
  )

ADD_LIBRARY(materials STATIC ${MATERIALS_SRCS})

TARGET_LINK_LIBRARIES(materials
  utils
  odesolver
  coeffunction
)
