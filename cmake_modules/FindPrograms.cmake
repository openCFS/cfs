# Find executables of a few required programs

#-----------------------------------------------------------------------------
# Find Python interpreter and libs
#-----------------------------------------------------------------------------
FIND_PACKAGE(PythonInterp)
FIND_PACKAGE(PythonLibs)

IF(NOT PYTHONINTERP_FOUND)
  MESSAGE(FATAL_ERROR "No Python interpreter was found! Please make sure a Python interpreter is on the PATH.")
ENDIF(NOT PYTHONINTERP_FOUND)

# This code has been taken from CMake 2.8.8 FindPythonInterp.cmake.
# determine python version string
if(PYTHON_EXECUTABLE)
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                            "import sys; sys.stdout.write(';'.join([str(x) for x in sys.version_info[:3]]))"
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
            string(REGEX REPLACE "\\.0$" "" PYTHON_VERSION_STRING "${PYTHON_VERSION_STRING}")
        endif()
    else()
        # sys.version predates sys.version_info, so use that
        execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c "import sys; sys.stdout.write(sys.version)"
                        OUTPUT_VARIABLE _VERSION
                        RESULT_VARIABLE _PYTHON_VERSION_RESULT
                        ERROR_QUIET)
        if(NOT _PYTHON_VERSION_RESULT)
            string(REGEX REPLACE " .*" "" PYTHON_VERSION_STRING "${_VERSION}")
            string(REGEX REPLACE "^([0-9]+)\\.[0-9]+.*" "\\1" PYTHON_VERSION_MAJOR "${PYTHON_VERSION_STRING}")
            string(REGEX REPLACE "^[0-9]+\\.([0-9])+.*" "\\1" PYTHON_VERSION_MINOR "${PYTHON_VERSION_STRING}")
            if(PYTHON_VERSION_STRING MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+.*")
                string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" PYTHON_VERSION_PATCH "${PYTHON_VERSION_STRING}")
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
