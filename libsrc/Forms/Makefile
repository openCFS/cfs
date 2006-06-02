# Makefile for discretization
C_DIR=../..
#
src = baseForm.cc laplaceInt.cc bdbInt.cc bdInt.cc massInt.cc baseoperator.cc\
      linViscoElastInt.cc linElastInt.cc nLinElastInt.cc\
      gradfieldop.cc elecforceop.cc curlCurlEdgeInt.cc massEdgeInt.cc\
      linearForm.cc curlfieldop.cc\
      smoothInt.cc curlCurlNodeInt.cc nLincurlCurlNodeInt.cc\
      mechStressStrain.cc linStrainOp.cc\
      elecchargeop.cc magforceop.cc adbInt.cc\
      linPiezoCoupling.cc linPressureInt.cc\
      baseForceOp.cc acouMechInt.cc abcInt.cc nrbcInt.cc linSurfForm.cc \
      linNeumannInt.cc PMLInt.cc laplaceXYZInt.cc PMLTransXYZInt.cc \
      acouPowerDensityOp.cc FlatShellInt.cc stokesFluidInt.cc stokesFluid3DInt.cc \
      stokesFluidPlaneInt.cc stokesFluidAxiInt.cc nLinStokesFluidInt.cc \
      linElecInt.cc
#
lib = forms
libpath = libsrc/Forms
CFLAGS3 =
#
include ../makefile.inc
-include $(DEPEND_FILE)
# DO NOT DELETE


