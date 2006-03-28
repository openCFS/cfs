#
# Makefile for discretization
#
C_DIR=../..
include ../../main/Makefile.option
#
src = filetype.cc timefunc.cc freqfunc.cc writeresults.cc \
	XMLMaterialHandler.cc PlainMaterialHandler.cc WriteInfo.cc\
	MHMaterialDataFile.cc \

#
ifeq ($(DATABASE),yes)
  src := $(src) LoadMaterialDataDatabase.cc
endif
#
#
lib = datainout
libpath = libsrc/DataInOut
CFLAGS3 =
#
include ../makefile.inc
-include $(DEPEND_FILE)
# DO NOT DELETE


