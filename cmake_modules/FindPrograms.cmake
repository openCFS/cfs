# Find executables of a few required programs

# Find the Subversion executable
FIND_PROGRAM(SVN svn)

# Find the Subversion executable
FIND_PROGRAM(SVNVERSION svnversion)

IF(NOT SVN OR
   NOT SVNVERSION)
  MESSAGE(SEND_ERROR "Subversion executable (svn) could not be found!")
  SUGGEST_INSTALL_PCKG()
ENDIF(NOT SVN OR
      NOT SVNVERSION)

FIND_PACKAGE(Doxygen)

IF(NOT DOXYGEN OR
   NOT DOT)
  MESSAGE("Doxygen could not be found! Note that you cannot build documentation!")
  SUGGEST_INSTALL_PCKG()
ENDIF(NOT DOXYGEN OR
      NOT DOT)

FIND_PACKAGE(LATEX)

IF(NOT LATEX_COMPILER)
  MESSAGE("LaTex could not be found! Note that you cannot build documentation!")
  SUGGEST_INSTALL_PCKG()
ENDIF(NOT LATEX_COMPILER)

MARK_AS_ADVANCED(
  SVN
  SVNVERSION
  )

IF(USE_ANSYSRST)
  # Set a sane default value for the ANSYS 10.0 path
  SET(ANSYS_ROOT_DIR "/home/data/programs/ansys/v100" CACHE STRING
    "Path to ANSYS installation.")
  
  IF(CFS_FORTRAN_COMPILER_NAME STREQUAL "GNU")
    IF(CFS_FORTRAN_COMPILER_VER LESS 4.1.0)
      MESSAGE(SEND_ERROR
	"The ANSYS RST writer requires a Fortran 95 compiler!")
    ENDIF(CFS_FORTRAN_COMPILER_VER LESS 4.1.0)
  ENDIF(CFS_FORTRAN_COMPILER_NAME STREQUAL "GNU")

  # Determine version of ANSYS.
  EXEC_PROGRAM("ls ${ANSYS_ROOT_DIR}/ansys/bin/ansys[0-9]*"
    ARGS
    OUTPUT_VARIABLE ANSYS_VERSION
    RETURN_VALUE RETVAL)
  
  STRING(REGEX REPLACE ".*ansys"
    "" ANSYS_VERSION
    ${ANSYS_VERSION})

  STRING(REGEX REPLACE "([0-9]$)"
    ".\\1" ANSYS_VERSION
    ${ANSYS_VERSION})

ENDIF(USE_ANSYSRST)
