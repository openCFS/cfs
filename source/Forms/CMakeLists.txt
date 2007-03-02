# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(FORMS_SRCS baseForm.cc laplaceInt.cc bdbInt.cc bdInt.cc massInt.cc
  baseoperator.cc linViscoElastInt.cc linElastInt.cc nLinElastInt.cc 
  gradfieldop.cc elecforceop.cc curlCurlEdgeInt.cc massEdgeInt.cc
  linearForm.cc curlfieldop.cc smoothInt.cc curlCurlNodeInt.cc 
  nLincurlCurlNodeInt.cc mechStressStrain.cc linStrainOp.cc elecchargeop.cc
  magforceop.cc adbInt.cc linPiezoCoupling.cc nLinPiezoCoupling.cc
  linPressureInt.cc baseForceOp.cc acouMechInt.cc abcInt.cc linSurfForm.cc
  linNeumannInt.cc PMLInt.cc laplaceXYZInt.cc 
  PMLTransXYZInt.cc acouPowerDensityOp.cc FlatShellInt.cc stokesFluidInt.cc
  stokesFluid3DInt.cc stokesFluidPlaneInt.cc stokesFluidAxiInt.cc
  nLinStokesFluidInt.cc linElecInt.cc FlatShellStiffInt.cc
  FlatShellMassInt.cc FlatShellPiezoInt.cc FlatShellElecInt.cc 
  singleEntryInt.cc nLinElecHystInt.cc bubbleDampInt.cc bubbleStiffInt.cc
  pierceInt.cc nLinElecMaterial.cc linSurfStressInt.cc nlinLaplaceInt.cc
  nLinMassInt.cc)

# nrbcInt.cc

ADD_LIBRARY(forms STATIC ${FORMS_SRCS})
