SET(UTILS_SRCS
  piezoMicroModelBK.cc
  mathParser/mathParser.cc
  mathfunctions.cc
  biotSavart.cc
  nodestoresol.cc
  elemstoresol.cc
  piezoMicroModel.cc
  Coil.cc
  interpolate.cc
  tools.cc
  cylCoordSys.cc
  polCoordSys.cc
  trivialCartesianCoordSys.cc
  dampLayer.cc
  cartesianCoordSys.cc
  SmoothSpline.cc
  preisach.cc
  coordSystem.cc
  simplePreisachInv.cc
  ApproxData.cc
  Timer.cc
  jiles.cc
  result.cc
  LinInterpolate.cc
  defaultCoordSys.cc
  Point.cc
  hysteresis.cc
  readCFDbinary.cc
  )


ADD_LIBRARY(utils STATIC ${UTILS_SRCS})

TARGET_LINK_LIBRARIES(utils
  cfsgeneral
  datainout
  domain
  elements
  logging
  materials
  matvec
  odesolver
  paramh
  pde
  ${BOOST_SERIALIZATION_LIB}
  ${BOOST_SIGNALS_LIB}
  ${MUPARSER_LIBRARY}
)

IF(USE_ANSYSRST)
  ADD_LIBRARY(CFSDynLd SHARED DynamicLibrary.cc)
  TARGET_LINK_LIBRARIES(CFSDynLd
    -ldl
  )
ENDIF(USE_ANSYSRST)

IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
set_source_files_properties(piezoMicroModelBK.cc
                              PROPERTIES COMPILE_FLAGS
                              -wd279)
ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")

IF(CFS_CXX_COMPILER_NAME STREQUAL "GCC")
set_source_files_properties(mathParser/mathParser.cc
                            PROPERTIES COMPILE_FLAGS
                            -fno-strict-aliasing)
ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "GCC")

