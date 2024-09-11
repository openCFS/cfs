#=============================================================================
# Find executables of a few required programs
#=============================================================================

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

# Find LaTeX to be used for doxygen formuls
find_package(LATEX)
if(NOT LATEX_COMPILER)
  message(STATUS "LaTex could not be found! Note that you cannot build documentation!")
endif()

# Find Pygmentize syntax highlighter.
find_program(PYGMENTIZE_EXECUTABLE NAMES pygmentize.cmd pygmentize)
mark_as_advanced(PYGMENTIZE_EXECUTABLE)
if(PYGMENTIZE_EXECUTABLE MATCHES "NOTFOUND")
  message(STATUS "command pygmentize not found, make doc might fail!")
endif()

# check if we have relevant tools installed
find_program(PATCH_EXECUTABLE patch)
find_program(FIND_EXECUTABLE find)
find_program(DIFF_EXECUTABLE diff)
find_program(UNZIP_EXECUTABLE unzip)

if(PATCH_EXECUTABLE MATCHES "NOTFOUND")
  message(WARNING "command patch not found, some cfsdeps will fail")
endif()

if(FIND_EXECUTABLE MATCHES "NOTFOUND")
  message(WARNING "command find not found, some cfsdeps (e.g. lis) will fail")
endif()

if(DIFF_EXECUTABLE MATCHES "NOTFOUND")
  message(WARNING "command diff not found, some cfsdeps (e.g. lis) will fail")
endif()

if(UNZIP_EXECUTABLE MATCHES "NOTFOUND" AND (USE_SNOPT OR USE_SCPIP))
  message(WARNING "command unzip not found, required for snopt and scpip")
endif()


mark_as_advanced(PATCH_EXECUTABLE)
mark_as_advanced(FIND_EXECUTABLE)
mark_as_advanced(DIFF_EXECUTABLE)
mark_as_advanced(UNZIP_EXECUTABLE)
