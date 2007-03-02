# Depending on USE_SCRIPTING_PYTHON, we have to include the
# python headers, as they get referenced already in the py-messenger
# header
IF(USE_SCRIPTING_PYTHON)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
ENDIF(USE_SCRIPTING_PYTHON)

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(DEFINEFILES_SRCS definefiles.cc)

ADD_LIBRARY(definefiles STATIC ${DEFINEFILES_SRCS})

