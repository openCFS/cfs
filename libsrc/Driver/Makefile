##
# Makefile for physical system
include ../../main/Makefile.option

C_DIR=../..
#
src =  	solveStepMech.cc solveStepMag.cc solveStepAcoustic.cc baseSolveStep.cc \
	solveStepAcouFlowNoise.cc solveStepMpCCI.cc solveStepSmooth.cc \
	solveStepElec.cc solveStepPiezo.cc basedriver.cc transientdriver.cc \
	staticdriver.cc harmonicDriver.cc multiSequenceDriver.cc\
	assemble.cc singleDriver.cc bubbledriver.cc \
#
ifeq ($(Adaptivity),yes)

  src = solveStepMech.cc solveStepMag.cc solveStepAcoustic.cc baseSolveStep.cc \
	solveStepAcouFlowNoise.cc solveStepMpCCI.cc solveStepSmooth.cc \
	solveStepElec.cc basedriver.cc transientdriver.cc staticdriver.cc \
	multiSequenceDriver.cc harmonicDriver.cc harmonicAdapt.cc staticAdaptSpaceDriver.cc \
	singleDriver.cc \


endif
#
lib = driver
libpath = libsrc/Driver
CFLAGS3 =
#
include ../makefile.inc
-include $(DEPEND_FILE)
# DO NOT DELETE


