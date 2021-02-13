#=============================================================================
# Find executables of a few required programs
#=============================================================================

# The pyhon situation is currently (10.2020) unpleasant. There is a new cmake module to find python
# find_package(Python COMPONENTS Interpreter Development)
# but this does not allow to guide the python to be found beside enforced version number
# see https://gitlab.kitware.com/cmake/cmake/-/issues/19492
# Note that it finds Python_* stuff. Python_EXECUTABLE, Python_VERSION, Python_LIBRARIES, Python_INCLUDE_DIRS 
# and also note, that the variables are case sensitive and different to the old PYTHON_EXECUTABLE, ... stuff

# this is here for TU-Wien, if it does not apply, it has no effect
# look for anaconda in mdmt environment
find_program(PYTHON_EXECUTABLE NAMES python python3 PATHS "/share/programs/anaconda3/2020.02/bin" "/share/programs/anaconda/latest/bin" "/share/programs/bin"  NO_DEFAULT_PATH)
message(STATUS "find_program found PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE} (CFS config)")

# identifies default "python" and sets PYTHON_EXECUTABLE. You can change this via -DPYTHON_EXECUTABLE=...
find_package(PythonInterp)
find_package(PythonInterp 3)

# sets PYTHON_LIBRARY and PYTHON_INCLUDE_DIR, can both the set via -DPYTHON_...
find_package(PythonLibs)

# Since DOXYGEN  and DOT are  not cache variables  in newer CMake  versions we
# have to set them from the new cache variables before testing them.
if(NOT DOXYGEN OR NOT DOT)
  set(DOT ${DOXYGEN_DOT_EXECUTABLE}) 
  set(DOXYGEN ${DOXYGEN_EXECUTABLE})
endif()

if(NOT DOXYGEN OR NOT DOT)
  find_package(Doxygen)
endif()


if(NOT DOXYGEN_EXECUTABLE)
  set(DOXYGEN_EXECUTABLE DOXYGEN)
endif()

# Find LaTeX.
find_package(LATEX)
if(NOT LATEX_COMPILER)
  set(MSG "LaTex could not be found! Note that you cannot build")
  set(MSG "${MSG} documentation!")
  message(WARNING "${MSG}")
endif()

# Find Pygmentize syntax highlighter.
find_program(PYGMENTIZE_EXECUTABLE NAMES pygmentize.cmd pygmentize)
mark_as_advanced(PYGMENTIZE_EXECUTABLE)

# check if we have relevant tools installed
find_program(PATCH_EXECUTABLE patch)
find_program(FIND_EXECUTABLE find)
find_program(DIFF_EXECUTABLE diff)

if(PATCH_EXECUTABLE MATCHES "NOTFOUND")
  message(WARNING "command patch not found, some cfsdeps will fail")
endif()
if(FIND_EXECUTABLE MATCHES "NOTFOUND")
  message(WARNING "command find not found, some cfsdeps (e.g. lis) will fail")
endif()
if(DIFF_EXECUTABLE MATCHES "NOTFOUND")
  message(WARNING "command diff not found, some cfsdeps (e.g. lis) will fail")
endif()

mark_as_advanced(PATCH_EXECUTABLE)
mark_as_advanced(FIND_EXECUTABLE)
mark_as_advanced(DIFF_EXECUTABLE)
