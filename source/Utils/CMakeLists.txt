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
  myclock.cc
  nodestoresol.cc
  polCoordSys.cc
  preisach.cc
  profiler.cc
  readCFDbinary.cc
  result.cc
  SmoothSpline.cc
  tools.cc
  trivialCartesianCoordSys.cc
  vector.cc
  )


ADD_LIBRARY(utils STATIC ${UTILS_SRCS})


SET(TARGET_LL 
  datainout
  matrix
  ${BOOST_SERIALIZATION_LIB}
)

TARGET_LINK_LIBRARIES(utils
  ${TARGET_LL}
)

ADD_LIBRARY(CFSDynLd SHARED DynamicLibrary.cc)
TARGET_LINK_LIBRARIES(CFSDynLd
  -ldl
)
  