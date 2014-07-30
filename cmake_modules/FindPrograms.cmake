#=============================================================================
# Find executables of a few required programs
#=============================================================================

#-----------------------------------------------------------------------------
# Find Python interpreter and libs
#-----------------------------------------------------------------------------
IF(MINGW)
  IF(NOT CMAKE_CROSSCOMPILING)
    FIND_PROGRAM(PYTHON_EXECUTABLE
      NAMES python)

    FIND_PROGRAM(PYTHON_CONFIG
      NAMES python-config.cmd python-config)
    MARK_AS_ADVANCED(PYTHON_CONFIG)

    IF(NOT PYTHON_INCLUDE_DIR AND PYTHON_CONFIG)
      EXECUTE_PROCESS(
	COMMAND ${PYTHON_CONFIG} --includes
	OUTPUT_VARIABLE PYINCDIRS
      )
      STRING(REPLACE "-I" "" PYINCDIRS "${PYINCDIRS}")
      STRING(REPLACE "\\" "/" PYINCDIRS "${PYINCDIRS}")
      STRING(REPLACE " " ";" PYINCDIRS "${PYINCDIRS}")
      STRING(STRIP "${PYINCDIRS}" PYINCDIRS)
      EXECUTE_PROCESS(
	COMMAND ${PYTHON_CONFIG} --prefix
	OUTPUT_VARIABLE PYPREFIX
      )
      STRING(REPLACE "\\" "/" PYPREFIX "${PYPREFIX}")
      STRING(STRIP "${PYPREFIX}" PYPREFIX)
      file(GLOB_RECURSE PYLIB "${PYPREFIX}/lib/*.dll.a")

      set(PYTHON_INCLUDE_DIR "${PYINCDIRS}" CACHE PATH
	"Path to where Python.h is found")
      set(PYTHON_LIBRARY "${PYLIB}" CACHE PATH
	"Path to Python library")
    ENDIF()
  ENDIF()
ENDIF()

FIND_PACKAGE(PythonInterp)
FIND_PACKAGE(PythonLibs)

IF(BUILD_PARAVIEW AND NOT PYTHONINTERP_FOUND)
  SET(MSG "No Python interpreter was found! Please make sure a Python")
  SET(MSG "${MSG} interpreter is on the PATH.")
  MESSAGE(FATAL_ERROR "${MSG}")
ENDIF(BUILD_PARAVIEW AND NOT PYTHONINTERP_FOUND)

#-----------------------------------------------------------------------------
# This  code  has  been  taken  from  CMake  2.8.8  FindPythonInterp.cmake  to
# determine the Python version string.
#-----------------------------------------------------------------------------
if(PYTHON_EXECUTABLE)
  set(PYPROG "import sys;
sys.stdout.write(';'.join([str(x) for x in sys.version_info[:3]]))")

  execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c "${PYPROG}"
    OUTPUT_VARIABLE _VERSION
    RESULT_VARIABLE _PYTHON_VERSION_RESULT
    ERROR_QUIET)
  if(NOT _PYTHON_VERSION_RESULT)
    string(REPLACE ";" "." PYTHON_VERSION_STRING "${_VERSION}")
    list(GET _VERSION 0 PYTHON_VERSION_MAJOR)
    list(GET _VERSION 1 PYTHON_VERSION_MINOR)
    list(GET _VERSION 2 PYTHON_VERSION_PATCH)
    if(PYTHON_VERSION_PATCH EQUAL 0)
      # it's called "Python 2.7", not "2.7.0"
      string(REGEX REPLACE "\\.0$" "" 
	PYTHON_VERSION_STRING "${PYTHON_VERSION_STRING}")
    endif()
  else()
    # sys.version predates sys.version_info, so use that
  set(PYPROG "import sys; sys.stdout.write(sys.version)")

  execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c "${PYPROG}"
    OUTPUT_VARIABLE _VERSION
    RESULT_VARIABLE _PYTHON_VERSION_RESULT
    ERROR_QUIET)
  if(NOT _PYTHON_VERSION_RESULT)
    string(REGEX REPLACE " .*" "" PYTHON_VERSION_STRING "${_VERSION}")
    string(REGEX REPLACE "^([0-9]+)\\.[0-9]+.*" "\\1"
      PYTHON_VERSION_MAJOR "${PYTHON_VERSION_STRING}")
    string(REGEX REPLACE "^[0-9]+\\.([0-9])+.*" "\\1"
      PYTHON_VERSION_MINOR "${PYTHON_VERSION_STRING}")
    if(PYTHON_VERSION_STRING MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+.*")
      string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1"
	PYTHON_VERSION_PATCH "${PYTHON_VERSION_STRING}")
    else()
      set(PYTHON_VERSION_PATCH "0")
    endif()
  else()
    # sys.version was first documented for Python 1.5, so assume
    # this is older.
    set(PYTHON_VERSION_STRING "1.4")
    set(PYTHON_VERSION_MAJOR "1")
    set(PYTHON_VERSION_MAJOR "4")
    set(PYTHON_VERSION_MAJOR "0")
  endif()
endif()
unset(_PYTHON_VERSION_RESULT)
unset(_VERSION)
endif(PYTHON_EXECUTABLE)

# MESSAGE("PYTHON_VERSION_STRING ${PYTHON_VERSION_STRING}")

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
