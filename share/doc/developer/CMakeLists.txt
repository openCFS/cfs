# FILE OBTAINED FROM JAN WOETZEL AND ADJUSTED BY
# SIMON TRIEBENBACHER FOR PROJECT CFS++
#
#  Dipl.-Ing. Jan Woetzel
#--------------------------------------------------
#  University of Kiel
#  Institute of Computer Science and Applied Mathematics
#  Hermann-Rodewald-Str. 3 [room 310]
#  24098 Kiel/Germany
#--------------------------------------------------
#  Phone +49-431-880-4477
#  Fax   +49-431-880-4845
#  Mob.  +49-179-2937346
#--------------------------------------------------
#  Url   www.mip.informatik.uni-kiel.de/~jw
#  Email jw at mip.informatik.uni-kiel.de

# -helper macro to add a "doc" target with CMake build system. 
# and configure doxy.config.in to doxy.config
#
# target "doc" allows building the documentation with doxygen/dot on WIN32 and Linux
# Creates .chm windows help file if MS HTML help workshop 
# (available from http://msdn.microsoft.com/workshop/author/htmlhelp)
# is installed with its DLLs in PATH.
#
#
# Please note, that the tools, e.g.:
# doxygen, dot, latex, dvips, makeindex, gswin32, etc.
# must be in path.
#
# Note about Visual Studio Projects: 
# MSVS hast its own path environment which may differ from the shell.
# See "Menu Tools/Options/Projects/VC++ Directories" in VS 7.1
#
# author Jan Woetzel 2004-2006
# www.mip.informatik.uni-kiel.de/~jw
cmake_minimum_required(VERSION 2.4)

IF (DOXYGEN)

  # click+jump in Emacs and Visual Studio (for doxy.config) (jw)
  IF    (CMAKE_BUILD_TOOL MATCHES "(msdev|devenv)")
    SET(DOXY_WARN_FORMAT "\"$file($line) : $text \"")
  ELSE  (CMAKE_BUILD_TOOL MATCHES "(msdev|devenv)")
    SET(DOXY_WARN_FORMAT "\"$file:$line: $text \"")
  ENDIF (CMAKE_BUILD_TOOL MATCHES "(msdev|devenv)")

  # we need latex for doxygen because of the formulas
  FIND_PACKAGE(LATEX)
  IF    (NOT LATEX_COMPILER)
    MESSAGE(STATUS "latex command LATEX_COMPILER not found but usually required. You will probably get warnings and user inetraction on doxy run.")
  ENDIF (NOT LATEX_COMPILER)
  IF    (NOT MAKEINDEX_COMPILER)
    MESSAGE(STATUS "makeindex command MAKEINDEX_COMPILER not found but usually required.")
  ENDIF (NOT MAKEINDEX_COMPILER)
  IF    (NOT DVIPS_CONVERTER)
    MESSAGE(STATUS "dvips command DVIPS_CONVERTER not found but usually required.")
  ENDIF (NOT DVIPS_CONVERTER)


  SET(DOXY_CONFIG_IN "${CMAKE_CURRENT_SOURCE_DIR}/DoxygenConfig/doxy-config.in")
  SET(DOXY_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/DoxygenConfig/doxy-config")

  
  IF (EXISTS ${DOXY_CONFIG_IN})
    MESSAGE(STATUS "configured ${DOXY_CONFIG_IN} --> ${DOXY_CONFIG}")
    CONFIGURE_FILE(${DOXY_CONFIG_IN} 
      ${DOXY_CONFIG}
      @ONLY )
    # use (configured) doxy.config from (out of place) BUILD tree:
    SET(DOXY_CONFIG ${DOXY_CONFIG})
  ELSE (EXISTS ${DOXY_CONFIG_IN})
    # use static hand-edited doxy.config from SOURCE tree:
    SET(DOXY_CONFIG ${DOXY_CONFIG})
    IF   (EXISTS ${DOXY_CONFIG})
      MESSAGE(STATUS "WARNING: using existing ${DOXY_CONFIG} instead of configuring from ${DOXY_CONFIG_IN} file.")
    ELSE (EXISTS ${DOXY_CONFIG})
      MESSAGE(SEND_ERROR "FATAL_ERROR: no existing ${DOXY_CONFIG} and no ${DOXY_CONFIG_IN} found.")
    ENDIF(EXISTS ${DOXY_CONFIG})
  ENDIF(EXISTS ${DOXY_CONFIG_IN})
  
  IF(WIN32)
    SET(GEN_DOCU_IN "${CMAKE_CURRENT_SOURCE_DIR}/makedoc.bat.in")
    SET(GEN_DOCU "${CMAKE_CURRENT_SOURCE_DIR}/makedoc.bat")
  ELSE(WIN32)
    SET(GEN_DOCU_IN "${CMAKE_CURRENT_SOURCE_DIR}/makedoc.sh.in")
    SET(GEN_DOCU "${CMAKE_CURRENT_SOURCE_DIR}/makedoc.sh")
  ENDIF(WIN32)

  CONFIGURE_FILE(${GEN_DOCU_IN}
    ${GEN_DOCU}
    @ONLY )
  MESSAGE(STATUS "configured ${GEN_DOCU_IN} --> ${GEN_DOCU}")

  ADD_CUSTOM_TARGET(doc-devel
                    ${GEN_DOCU} "build"
		    )

  ADD_CUSTOM_TARGET(doc-devel-rebuild
                    ${GEN_DOCU} "rebuild"
		    )

  ADD_CUSTOM_TARGET(doc-devel-clean
                    ${GEN_DOCU} "clean"
		    )

  # create a windows help .chm file using hhc.exe
  # HTMLHelp DLL must be in path!
  # fallback: use hhw.exe interactively
  IF    (WIN32)
    FIND_PACKAGE(HTMLHelp)
    IF   (HTML_HELP_COMPILER)      
      SET (TMP "${CMAKE_CURRENT_BINARY_DIR}\\html\\index.hhp")
      STRING(REGEX REPLACE "[/]" "\\\\" HHP_FILE ${TMP} )
      # MESSAGE(SEND_ERROR "DBG  HHP_FILE=${HHP_FILE}")
      ADD_CUSTOM_TARGET(winhelp ${HTML_HELP_COMPILER} ${HHP_FILE})
      ADD_DEPENDENCIES (winhelp docu)
     
      IF (NOT TARGET_DOC_SKIP_INSTALL)
      # install windows help?
      # determine useful name for output file 
      # should be project and version unique to allow installing 
      # multiple projects into one global directory      
      IF   (EXISTS "${PROJECT_BINARY_DIR}/html/index.chm")
        IF   (PROJECT_NAME)
          SET(OUT "${PROJECT_NAME}")
        ELSE (PROJECT_NAME)
          SET(OUT "Documentation") # default
        ENDIF(PROJECT_NAME)
        IF   (${PROJECT_NAME}_VERSION_MAJOR)
          SET(OUT "${OUT}-${${PROJECT_NAME}_VERSION_MAJOR}")
          IF   (${PROJECT_NAME}_VERSION_MINOR)
            SET(OUT  "${OUT}.${${PROJECT_NAME}_VERSION_MINOR}")
            IF   (${PROJECT_NAME}_VERSION_PATCH)
              SET(OUT "${OUT}.${${PROJECT_NAME}_VERSION_PATCH}")      
            ENDIF(${PROJECT_NAME}_VERSION_PATCH)
          ENDIF(${PROJECT_NAME}_VERSION_MINOR)
        ENDIF(${PROJECT_NAME}_VERSION_MAJOR)
        # keep suffix
        SET(OUT  "${OUT}.chm")
        
        #MESSAGE("DBG ${PROJECT_BINARY_DIR}/html/index.chm \n${OUT}")
        # create target used by install and package commands 
        INSTALL(FILES "${PROJECT_BINARY_DIR}/html/index.chm"
          DESTINATION "doc"
          RENAME "${OUT}"
        )
      ENDIF(EXISTS "${PROJECT_BINARY_DIR}/html/index.chm")
      ENDIF(NOT TARGET_DOC_SKIP_INSTALL)

    ENDIF(HTML_HELP_COMPILER)
    # MESSAGE(SEND_ERROR "HTML_HELP_COMPILER=${HTML_HELP_COMPILER}")
  ENDIF (WIN32) 
ELSE (DOXYGEN)
  MESSAGE(STATUS "Doxygen not found")
ENDIF(DOXYGEN)
