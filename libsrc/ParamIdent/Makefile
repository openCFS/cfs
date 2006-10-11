##
# Makefile for physical system
include ../../main/Makefile.option

C_DIR=../..
#
src =  	piezoParamIdent.cc piezoParamIdent_methods.cc \
		nuMethods.cc nonlinLandweber.cc\
	    optimalExpDesign.cc optimalExpDesignVarNrFreqs.cc leastSquare.cc \
#
ifeq ($(Adaptivity),yes)

  src = piezoParamIdent.cc piezoParamIdent_methods.cc \
	 nuMethods.cc nonlinLandweber.cc\
	optimalExpDesign.cc optimalExpDesignVarNrFreqs.cc leastSquare.cc\
#
endif
#
CFLAGS3 =
#
include ../makefile.inc
-include $(DEPEND_FILE)
# DO NOT DELETE
