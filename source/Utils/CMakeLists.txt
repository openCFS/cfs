SET(UTILS_SRCS
  ApproxData.cc
  cartesianCoordSys.cc
  Coil.cc
  coordSystem.cc
  cylCoordSys.cc
  dampLayer.cc
  defaultCoordSys.cc
  elemstoresol.cc
  hysteresis.cc
  interpolate.cc
  jiles.cc
  LinInterpolate.cc
  mathfunctions.cc
  Timer.cc
  nodestoresol.cc
  polCoordSys.cc
  preisach.cc
  readCFDbinary.cc
  result.cc
  SmoothSpline.cc
  tools.cc
  Point.cc
  trivialCartesianCoordSys.cc
  piezoMicroModel.cc
  piezoMicroModelBK.cc
  )


ADD_LIBRARY(utils STATIC ${UTILS_SRCS})

SET(TARGET_LL 
  datainout
  matvec
  ${BOOST_SERIALIZATION_LIB}
)

TARGET_LINK_LIBRARIES(utils
  ${TARGET_LL}
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
