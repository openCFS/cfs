#
C_DIR = ../../..
include $(C_DIR)/main/Makefile.option
ifeq "$(TCL_INTERFACE)" "yes"
  src = cfsmessenger.cc tcl-messenger.cc scriptable.cc
else
  src = scriptable.cc cfsmessenger.cc
endif
#
lib = scripting
libpath = libsrc/DataInOut/Scripting
CFLAGS3 =
#
include ../../makefile.inc
-include $(DEPEND_FILE)
# DO NOT DELETE


