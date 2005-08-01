##
# Makefile for physical system
include ../../main/Makefile.option

C_DIR=../..
#
src =  	piezoParamIdent.cc piezoParamIdent_methods.cc \
	newtonCG.cc newtonLandweber.cc nuMethods.cc optimalExpDesign.cc \
	optimalExpDesign.cc optimalExpDesignVarNrFreqs.cc \
#
ifeq ($(Adaptivity),yes)

  src = piezoParamIdent.cc piezoParamIdent_methods.cc \
	newtonCG.cc newtonLandweber.cc nuMethods.cc \
	optimalExpDesign.cc optimalExpDesignVarNrFreqs.cc \
	
endif
#
lib = paramIdent
libpath = libsrc/ParamIdent
CFLAGS3 =
#
include ../makefile.inc
-include $(DEPEND_FILE)
# DO NOT DELETE
