#
C_DIR= ../..
include ../../main/Makefile.option
#
ifeq ($(USE_NEWBASEPDE),yes)
src =   newBasePDE.cc\
	newElecPDE.cc\
	newMechPDE.cc\
	newAcousticPDE.cc\
	newmark.cc\
	timestepping.cc\
	BaseEQN.cc\
	ScalarNodeEQN.cc\
	newsmoothPDE.cc\
	magneticPDE.cc\
	trapezoidal.cc\
	acouFlowNoise.cc\
	piezoPDE.cc
else
src =   basepde.cc\
	elecPDE.cc\
	acousticPDE.cc\
	acouFlowNoise.cc\
	mechPDE.cc\
	magEdgePDE.cc\
	smoothPDE.cc\
	smoothlaplacePDE.cc\
	timestepping.cc\
	newmark.cc\
	newmarkdamp.cc
endif
#
lib = pde
CFLAGS3 =
#
include ../makefile.inc

# DO NOT DELETE


