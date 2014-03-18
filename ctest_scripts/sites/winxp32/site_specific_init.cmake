MESSAGE(
"
=============================================================================
 Making a local copy of WC since we use SlikSVN 1.7 on Windows...
=============================================================================
"
)

SET(OLD_WC_DIR "z:/CFS_FESPACE_NIGHTLY")
SET(NEW_WC_DIR "c:/CFS_FESPACE_NIGHTLY")

IF(EXISTS ${NEW_WC_DIR})
  FILE(REMOVE_RECURSE ${NEW_WC_DIR})
  FILE(MAKE_DIRECTORY ${NEW_WC_DIR})
ENDIF()

EXECUTE_PROCESS(
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${OLD_WC_DIR} ${NEW_WC_DIR}
  RESULT_VARIABLE RETVAL
  )

MESSAGE(
"
=============================================================================
 Upgrading local SVN working copy...
=============================================================================
"
)
EXECUTE_PROCESS(
  COMMAND svn upgrade
  WORKING_DIRECTORY ${NEW_WC_DIR}
  RESULT_VARIABLE RETVAL
  )

MESSAGE(
"
=============================================================================
 Creating batch scripts for calling GCC 32-bit...
=============================================================================
"
)
FILE(REMOVE_RECURSE "c:/gcc32")
FILE(WRITE
  "c:/gcc32/gcc.bat"
"
C:\\MinGW64\\bin\\gcc.exe -m32 %*
"
)

FILE(WRITE
  "c:/gcc32/g++.bat"
"
C:\\MinGW64\\bin\\g++.exe -m32 %*
"
)

FILE(WRITE
  "c:/gcc32/gfortran.bat"
"
C:\\MinGW64\\bin\\gfortran.exe -m32 %*
"
)

