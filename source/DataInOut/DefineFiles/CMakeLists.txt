# Depending on USE_SCRIPTING_TCL/PYTHON, we have to include the
# tcl/python headers, as they get referenced already in the
# tcl/py-messenger header
IF(USE_TCL)
  INCLUDE_DIRECTORIES(${TCL_INCLUDE_PATH})
ENDIF(USE_TCL)

IF(USE_PYTHON)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
ENDIF(USE_PYTHON)


SET(DEFINEFILES_SRCS definefiles.cc)

ADD_LIBRARY(definefiles STATIC ${DEFINEFILES_SRCS})

