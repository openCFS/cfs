# - Find SCPIP installation 
# This module finds if SCPIP is installed and determines where 
# the libraries are.
#
#  SCPIP_FOUND     = system has SCPIP 
#  SCPIP_LIBRARY = path to the SCPIP libraries
#
# Fabian Wein (11/07)
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Please note that one needs an academic use license from 
# PD. Christian Zillober (Uni-Würzburg) before using SCPIP!!

SET(SCPIP_FOUND 0)

FIND_LIBRARY(SCPIP_LIBRARY
  NAMES scpip
  PATHS ${CFSDEPS_LIBRARY_DIR}
)

IF(SCPIP_LIBRARY)
  SET(SCPIP_FOUND 1)
ENDIF(SCPIP_LIBRARY)

MARK_AS_ADVANCED(SCPIP_LIBRARY)

IF(NOT SCPIP_FOUND)
  MESSAGE("Warning: SCPIP could not be found! Please specify proper paths.")
ENDIF(NOT SCPIP_FOUND)
