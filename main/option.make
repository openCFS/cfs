program = cfs
appl =
applpath = main
Dependencies = dependencies.make
#
############ Directories
#
LIBSRC_DIR=$(C_DIR)/libsrc
#
LIB_DIR=$(C_DIR)/lib/$(MACHINE)
#
########### Flags
#
INCL =  -I$(LIBSRC_DIR)/General -I$(LIBSRC_DIR)/DataInOut -I$(LIBSRC_DIR)/DataInOut/DatFile \
          -I$(LIBSRC_DIR)/DataInOut/Unverg -I$(LIBSRC_DIR)/AlgebraicSystem/LinAlg \
          -I$(LIBSRC_DIR)/Domain \
    -I$(LIBSRC_DIR)/Driver -I$(LIBSRC_DIR)/Elements -I$(LIBSRC_DIR)/Forms \
          -I$(LIBSRC_DIR)/PDE -I$(LIBSRC_DIR)/Utils -I$(LIBSRC_DIR)/Matrix
CFLAGS1 = -c $(INCL)
#
FFLAGS1 = -c
#
LINKFLAGS3 =
#
FlagMakeDepend = -a -v  
#
