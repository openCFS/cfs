#
C_DIR= ../..
#
include ../../main/Makefile.option

ifeq ($(MpCCI),yes)
  src = MpCCIexch.cc
else
  src = 
endif
#
lib = mpccicpl
CFLAGS3 =
#
include ../makefile.inc
ifeq ($(MpCCI),yes)
  -include $(Dependencies)
else
endif
# DO NOT DELETE


