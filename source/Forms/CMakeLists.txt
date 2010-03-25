# Just redefine __builtin_isnan to ::isnan for Intel compiler for
# nonConformingInt.cc
IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
    GET_SOURCE_FILE_PROPERTY(OPTCONDFLAGS abcInt.cc
        COMPILE_FLAGS)
    #    MESSAGE("OPTCONDFLAGS ${OPTCONDFLAGS}")
    IF(OPTCONDFLAGS STREQUAL NOTFOUND)
      SET(OPTCONDFLAGS "")
    ENDIF(OPTCONDFLAGS STREQUAL NOTFOUND)
    
    SET_SOURCE_FILES_PROPERTIES(nonConformingInt.cc
      PROPERTIES COMPILE_FLAGS "${OPTCONDFLAGS} -D__builtin_isnan=::isnan -D__builtin_isinf=::isinf")
ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")


SET(FORMS_SRCS
  abcInt.cc
  acouMechInt.cc
  acouMechNcInt.cc
  acouPowerDensityOp.cc
  acouRHSLinForm.cc
  adbInt.cc
  baseForceOp.cc
  baseForm.cc
  baseoperator.cc
  bdbInt.cc
  bdInt.cc
  curlCurlEdgeInt.cc
  curlCurlNodeInt.cc
  curlfieldop.cc
  elecchargeop.cc
  elecforceop.cc
  FlatShellElecInt.cc
  FlatShellInt.cc
  FlatShellMassInt.cc
  FlatShellPiezoInt.cc
  FlatShellStiffInt.cc
  gradfieldop.cc
  laplaceInt.cc
  linearForm.cc
  linElastInt.cc
  linElecInt.cc
  linNeumannInt.cc
  linPiezoCoupling.cc
  linPressureInt.cc
  linStrainOp.cc
  linSurfForm.cc
  linSurfStressInt.cc
  linThermoElectricInt.cc 
  linThermoMechInt.cc
  linThermoElecDampInt.cc 
  linThermoMechDampInt.cc 
  linHeatCondInt.cc
  linViscoElastInt.cc
  magforceop.cc
  massEdgeInt.cc
  massInt.cc
  mechStressStrain.cc
  nLincurlCurlNodeInt.cc
  nLinElastInt.cc
  nLinElecHystInt.cc
  nLinElecMaterial.cc
  nlinLaplaceInt.cc
  nLinMassInt.cc
  nLinPiezoCoupling.cc
  nLinPiezoHyst.cc
  nLinPiezoHystRHS.cc
  nonConformingInt.cc
  pierceInt.cc
  piezoPolarizationMatrixRHSInt.cc
  pmlBasics.cc
  PMLInt.cc
  singleEntryInt.cc
  smoothInt.cc
  smoothNLInt.cc
  nLinMagHystInt2D.cc
  mixedInt.cc
  fluidMechInt.cc
  fluidMechShearStress.cc
  fluidMechMassInt.cc
  fluidMechStiffInt.cc
  nLinPiezoMicro.cc
  nLinPiezoMicroRHS.cc
  )

ADD_LIBRARY(forms STATIC ${FORMS_SRCS})


TARGET_LINK_LIBRARIES(forms
  domain
  materials
  mathparser
  matvec
  utils
)
