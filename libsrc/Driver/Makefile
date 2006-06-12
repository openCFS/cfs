##
# Makefile for physical system
include ../../main/Makefile.option

C_DIR=../..
#
src =  	solveStepAcoustic.cc stdSolveStep.cc \
	solveStepMpCCI.cc solveStepSmooth.cc \
	baseSolveStep.cc solveStepElec.cc basedriver.cc transientdriver.cc \
	staticdriver.cc harmonicDriver.cc multiSequenceDriver.cc\
	assemble.cc singleDriver.cc transientHarmonicDriver.cc\
	iterSolveStep.cc bubbledriver.cc multiHarmonicDriver.cc \
	solveStepAcousticBubble.cc solveStepPiezo.cc MHassemble.cc \
	solveStepAcousticMechBubble.cc \
	solveStepStokesFluid.cc eigenFrequencyDriver.cc solveStepODE.cc
#
ifeq ($(Adaptivity),yes)

  src = solveStepAcoustic.cc stdSolveStep.cc \
	solveStepMpCCI.cc solveStepSmooth.cc \
	baseSolveStep.cc solveStepElec.cc basedriver.cc transientdriver.cc staticdriver.cc \
	multiSequenceDriver.cc harmonicDriver.cc harmonicAdapt.cc staticAdaptSpaceDriver.cc \
	piezoParamIdent.cc singleDriver.cc \
	iterSolveStep.cc multiHarmonicDriver.cc solveStepAcousticBubble.cc MHassemble.cc \
	solveStepStokesFluid.cc eigenFrequencyDriver.cc solveStepODE.cc


endif
#
lib = driver
libpath = libsrc/Driver
CFLAGS3 =
#
include ../makefile.inc
-include $(DEPEND_FILE)
# DO NOT DELETE


