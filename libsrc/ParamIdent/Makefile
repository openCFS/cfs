#
# Makefile for physical system
include ../../main/Makefile.option

C_DIR=../..
#
src = basedriver.cc transientdriver.cc staticdriver.cc harmonicDriver.cc multiSequenceDriver.cc assemble.cc piezoParamIdent.cc
#
ifeq ($(Adaptivity),yes)
  src = basedriver.cc transientdriver.cc staticdriver.cc multiSequenceDriver.cc harmonicDriver.cc harmonicAdapt.cc staticAdaptSpaceDriver.cc piezoParamIdent.cc
endif
#
lib = driver
libpath = libsrc/Driver
CFLAGS3 =
#
include ../makefile.inc

# DO NOT DELETE


