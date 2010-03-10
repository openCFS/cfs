# Make sure the compiler can find include files from TCL.
INCLUDE_DIRECTORIES(${TCL_INCLUDE_PATH})

SET(TCLMESS_SRCS tcl-messenger.cc)

ADD_LIBRARY(tclmess STATIC ${TCLMESS_SRCS})

TARGET_LINK_LIBRARIES(tclmess ${TCL_LIBRARY})
