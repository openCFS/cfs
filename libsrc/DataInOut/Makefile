#
# Makefile for discretization
#
C_DIR=../..
include ../../main/Makefile.option
#
src = filetype.cc timefunc.cc writeresults.cc MaterialData.cc\
	LoadMaterialData.cc LoadMaterialDataFile.cc WriteInfo.cc
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


