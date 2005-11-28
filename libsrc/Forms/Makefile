# Makefile for discretization
C_DIR=../..
#
src = baseForm.cc laplaceInt.cc bdbInt.cc linElastInt.cc nLinElastInt.cc linViscoElastInt.cc bdInt.cc\
      massInt.cc baseoperator.cc\
      gradfieldop.cc elecforceop.cc curlCurlEdgeInt.cc massEdgeInt.cc\
      linearForm.cc curlfieldop.cc\
      smoothInt.cc curlCurlNodeInt.cc linPiezoInt.cc nLincurlCurlNodeInt.cc\
      mechStressStrain.cc linStrainOp.cc\
      elecchargeop.cc magforceop.cc adbInt.cc\
      linPiezoCoupling3D.cc linPiezoCoupling2DAxi.cc linPressureInt.cc\
      baseForceOp.cc acouMechInt.cc abcInt.cc linSurfForm.cc acouNeumannInt.cc \
      laplacePMLInt.cc PMLInt.cc
#
lib = forms
libpath = libsrc/Forms
CFLAGS3 =
#
include ../makefile.inc
-include $(DEPEND_FILE)
# DO NOT DELETE


