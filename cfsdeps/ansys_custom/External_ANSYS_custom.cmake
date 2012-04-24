#-------------------------------------------------------------------------------
# Set paths to ansys_custom sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(ansys_custom_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/ansys_custom")
set(ansys_custom_source  "${ansys_custom_prefix}/src/ansys_custom")
set(ansys_custom_install  "${CMAKE_CURRENT_BINARY_DIR}")

#-------------------------------------------------------------------------------
# Set names of build file and template file.
#-------------------------------------------------------------------------------
SET(BUILD_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/ansys_custom/ansys_custom-build.cmake.in")
SET(BUILD "${ansys_custom_prefix}/ansys_custom-build.cmake")
CONFIGURE_FILE("${BUILD_TEMPL}" "${BUILD}" @ONLY) 

#-------------------------------------------------------------------------------
# Set names of install file and template file.
#-------------------------------------------------------------------------------
SET(INST_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/ansys_custom/ansys_custom-install.cmake.in")
SET(INST "${ansys_custom_prefix}/ansys_custom-install.cmake")
CONFIGURE_FILE("${INST_TEMPL}" "${INST}" @ONLY) 

#-------------------------------------------------------------------------------
# The ansys_custom external project
#-------------------------------------------------------------------------------
ExternalProject_Add(ansys_custom
  PREFIX "${ansys_custom_prefix}"
  SOURCE_DIR "${ansys_custom_source}"
  BUILD_IN_SOURCE 1
  DOWNLOAD_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND  ${CMAKE_COMMAND} -P ${BUILD}
  INSTALL_COMMAND  ${CMAKE_COMMAND} -P ${INST}
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  ansys_custom
)

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
#SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/ansys_custom/ansys_custom-patch.cmake.in")
#SET(PFN "${ansys_custom_prefix}/ansys_custom-patch.cmake")
#CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# We do not use the PATCH_COMMAND  of ExternalProject_Add since we do not only
# want to apply the patch script  during configuration time but also if it has
# changed.  Therefore,   we  need  a   dependency  on  the   configured  patch
# script. This can be achieved by  adding an additional build step between the
# download and configure steps.
#
# NOTE: The  patch script should  be designed  in such a  way, that it  can be
# applied to  an already patched  source tree. This  is due to the  fact, that
# ExternalProject_Add only extracts the source if the MD5 sum has has changed.
#-------------------------------------------------------------------------------
#ExternalProject_Add_Step(ansys_custom custom_patch
#   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
#   DEPENDEES download
#   DEPENDERS configure
#   DEPENDS "${PFN}"
#   WORKING_DIRECTORY ${ansys_custom_source}
#)

#-------------------------------------------------------------------------------
# Determine paths of ANSYS Customization libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")

IF(CFS_ARCH STREQUAL "I386")
  SET(ANSYS_V100_LIBBIN
    "${LD}/v100/linia32/libbin.so"
    CACHE FILEPATH "ANSYS binary interface library.")
  SET(ANSYS_V110_LIBBIN
    "${LD}/v110/linia32/libbin.so"
    CACHE FILEPATH "ANSYS binary interface library.")
ENDIF(CFS_ARCH STREQUAL "I386")

IF(CFS_ARCH STREQUAL "X86_64")
  IF(CFS_SUBARCH STREQUAL "OPTERON")
    SET(ANSYS_V100_LIBBIN
      "${LD}/v100/linop64/libbin.so"
      CACHE FILEPATH "ANSYS binary interface library.")
    SET(ANSYS_V110_LIBBIN
      "${LD}/v110/linop64/libbin.so"
      CACHE FILEPATH "ANSYS binary interface library.")
  ENDIF(CFS_SUBARCH STREQUAL "OPTERON")

  IF(CFS_SUBARCH STREQUAL "EM64T")
    SET(ANSYS_V100_LIBBIN
      "${LD}/v100/linem64t/libbin.so"
      CACHE FILEPATH "ANSYS binary interface library.")
    SET(ANSYS_V110_LIBBIN
      "${LD}/v110/linem64t/libbin.so;-lguide;-lpthread"
      CACHE FILEPATH "ANSYS binary interface library.")
  ENDIF(CFS_SUBARCH STREQUAL "EM64T")
  
  SET(ANSYS_V120_LIBBIN
    "${LD}/v120/libbin.so;-lpthread"
    CACHE FILEPATH "ANSYS binary interface library.")

  SET(ANSYS_V130_LIBBIN
    "${LD}/v130/libbin.so;-lpthread"
    CACHE FILEPATH "ANSYS binary interface library.")
  
ENDIF(CFS_ARCH STREQUAL "X86_64")

SET(ANSYS_V100_INCLUDE_DIR "${CFS_BINARY_DIR}/include/v100")
SET(ANSYS_V110_INCLUDE_DIR "${CFS_BINARY_DIR}/include/v110")
SET(ANSYS_V120_INCLUDE_DIR "${CFS_BINARY_DIR}/include/v120")
SET(ANSYS_V130_INCLUDE_DIR "${CFS_BINARY_DIR}/include/v130")

MARK_AS_ADVANCED(ANSYS_V100_LIBBIN)
MARK_AS_ADVANCED(ANSYS_V110_LIBBIN)
MARK_AS_ADVANCED(ANSYS_V120_LIBBIN)
MARK_AS_ADVANCED(ANSYS_V130_LIBBIN)

