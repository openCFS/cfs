#
# Makefile for discretization
#
C_DIR=../../..
#
include ../../../main/Makefile.option

ifeq ($(Adaptivity),yes)
   src = interface_adgrid.cc 
else
   src = 
endif
#
lib = toadaptgrid
libpath = libsrc/Domain/AdaptGrid
CFLAGS3 =
#
include ../../makefile.inc
ifeq ($(Adaptivity),yes)
  -include $(Dependencies)
else
endif
# DO NOT DELETE


