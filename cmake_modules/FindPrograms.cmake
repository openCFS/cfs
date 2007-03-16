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
