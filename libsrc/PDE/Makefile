#
C_DIR= ../..
include ../../main/Makefile.option
#
ifeq ($(CFLAGS2),-DNEWBASEPDE)
src = newBasePDE.cc newElecPDE.cc newMechPDE.cc newmark.cc timestepping.cc
else
src = basepde.cc elecPDE.cc acousticPDE.cc acouFlowNoise.cc mechPDE.cc magEdgePDE.cc smoothPDE.cc smoothlaplacePDE.cc \
      timestepping.cc newmark.cc newmarkdamp.cc
endif
#
lib = pde
CFLAGS3 =
#
include ../makefile.inc

# DO NOT DELETE


