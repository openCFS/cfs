##
# Makefile for physical system
include ../../main/Makefile.option

C_DIR=../..
#
src =  	solveStepAcoustic.cc stdSolveStep.cc \
	solveStepMpCCI.cc \
	baseSolveStep.cc solveStepElec.cc basedriver.cc transientdriver.cc \
	staticdriver.cc harmonicDriver.cc multiSequenceDriver.cc\
	assemble.cc singleDriver.cc transientHarmonicDriver.cc\
	iterSolveStep.cc \
	solveStepPiezo.cc \
	solveStepStokesFluid.cc eigenFrequencyDriver.cc formsContext.cc \
	solveStepODE.cc
#
ifeq ($(Adaptivity),yes)

  src = solveStepAcoustic.cc stdSolveStep.cc \
	solveStepMpCCI.cc\
	baseSolveStep.cc solveStepElec.cc basedriver.cc transientdriver.cc staticdriver.cc \
	multiSequenceDriver.cc harmonicDriver.cc harmonicAdapt.cc staticAdaptSpaceDriver.cc \
	piezoParamIdent.cc singleDriver.cc \
	iterSolveStep.cc  solveStepAcousticBubble.cc\
	solveStepStokesFluid.cc eigenFrequencyDriver.cc solveStepODE.cc


endif
#
CFLAGS3 =
#
include ../makefile.inc
-include $(DEPEND_FILE)
# DO NOT DELETE


