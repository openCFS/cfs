#
C_DIR = ../../..
include ../../../main/Makefile.option
#
ifeq ($(DATABASE),yes)
src = database.cc outDB.cc dbColumn.cc dbIntColumn.cc dbDoubleColumn.cc dbStringColumn.cc dbMatrix.cc dbLineData.cc
else
src = 
endif
#
lib = database
libpath = libsrc/DataInOut/Database
CFLAGS3 =
#
include ../../makefile.inc
ifeq ($(DATABASE),yes)
  -include $(Dependencies)
endif
# DO NOT DELETE


