#=============================================================================
# Find executables of a few required programs
#=============================================================================

# preferres python3 over python2: https://cmake.org/cmake/help/latest/module/FindPython.html#module:FindPython
# Finds Python_* stuff. Python_EXECUTABLE, Python_VERSION, Python_LIBRARIES, Python_INCLUDE_DIRS
# find_package(Python COMPONENTS Interpreter Development)

# identifies "python" and sets PYTHON_EXECUTABLE which can be different from Python_EXECUTABLE
find_package(PythonInterp)

# find_package(PythonLibs) is in FindCFSDEPS


#-----------------------------------------------------------------------------
# Since DOXYGEN  and DOT are  not cache variables  in newer CMake  versions we
# have to set them from the new cache variables before testing them.
#-----------------------------------------------------------------------------
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

#-----------------------------------------------------------------------------
# Find LaTeX.
#-----------------------------------------------------------------------------
FIND_PACKAGE(LATEX)

IF(NOT LATEX_COMPILER)
  SET(MSG "LaTex could not be found! Note that you cannot build")
  SET(MSG "${MSG} documentation!")
  MESSAGE(WARNING "${MSG}")
ENDIF(NOT LATEX_COMPILER)

#-----------------------------------------------------------------------------
# Find Pygmentize syntax highlighter.
#-----------------------------------------------------------------------------
FIND_PROGRAM(PYGMENTIZE_EXECUTABLE
  NAMES pygmentize.cmd pygmentize)
MARK_AS_ADVANCED(PYGMENTIZE_EXECUTABLE)

#-----------------------------------------------------------------------------
# Find patch command.
#-----------------------------------------------------------------------------
FIND_PROGRAM(PATCH_EXECUTABLE patch)
MARK_AS_ADVANCED(PATCH_EXECUTABLE)


