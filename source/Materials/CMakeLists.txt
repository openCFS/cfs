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
  Models/VectorPreisachv10.cc
  #Models/VectorPreisachv7.cc
  #Models/VectorPreisach.cc
  Models/Preisach.cc
  Models/SimplePreisachInv.cc 
  )

ADD_LIBRARY(materials STATIC ${MATERIALS_SRCS})

TARGET_LINK_LIBRARIES(materials
  utils
  odesolver
)
IF(USE_ADOLC)
  INCLUDE_DIRECTORIES(${ADOLC_INCLUDE_DIR})
  #LINK_DIRECTORIES(${ADOLC_LIBRARY})
ENDIF(USE_ADOLC)


