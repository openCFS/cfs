# Find executables of a few required programs

#-----------------------------------------------------------------------------
# Find Perl interpreter and configure a Perl source file with the current
# settings for compilers and tools for building dependencies
#-----------------------------------------------------------------------------
FIND_PACKAGE(Perl)

FIND_PACKAGE(PythonInterp)
FIND_PACKAGE(PythonLibs)

IF(NOT PYTHONINTERP_FOUND)
  MESSAGE(FATAL_ERROR "No Python interpreter was found! Please make sure a Python interpreter is on the PATH.")
ENDIF(NOT PYTHONINTERP_FOUND)


# Since DOXYGEN and DOT are not cache variables in newer CMake versions
# we have to set them from the new cache variables before testing them.
IF(NOT DOXYGEN OR NOT DOT)
  SET(DOT ${DOXYGEN_DOT_EXECUTABLE}) 
  SET(DOXYGEN ${DOXYGEN_EXECUTABLE})
ENDIF(NOT DOXYGEN OR NOT DOT)

IF(NOT DOXYGEN OR NOT DOT)
 FIND_PACKAGE(Doxygen)
ENDIF(NOT DOXYGEN OR NOT DOT)


IF(NOT DOXYGEN_EXECUTABLE)
  SET(DOXYGEN_EXECUTABLE DOXYGEN)
ENDIF(NOT DOXYGEN_EXECUTABLE)

FIND_PACKAGE(LATEX)

IF(NOT LATEX_COMPILER)
  MESSAGE("Warning: LaTex could not be found! Note that you cannot build documentation!")
ENDIF(NOT LATEX_COMPILER)

MARK_AS_ADVANCED(
  SVN
  SVNVERSION
  )

FIND_PROGRAM(MATLAB_EXECUTABLE matlab)
MARK_AS_ADVANCED(MATLAB_EXECUTABLE)
