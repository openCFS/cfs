SET(UTILS_SRCS myclock.cc nodestoresol.cc elemstoresol.cc tools.cc vector.cc
  Coil.cc ApproxData.cc SmoothSpline.cc mathfunctions.cc
  hysteresis.cc preisach.cc LinInterpolate.cc jiles.cc profiler.cc
  coordSystem.cc cylCoordSys.cc defaultCoordSys.cc dampLayer.cc
  interpolate.cc result.cc polCoordSys.cc readCFDbinary.cc
  cartesianCoordSys.cc trivialCartesianCoordSys.cc)

ADD_LIBRARY(utils STATIC ${UTILS_SRCS})


SET(TARGET_LL 
  datainout
  matrix
  ${BOOST_SERIALIZATION_LIB}
)

TARGET_LINK_LIBRARIES(utils
  ${TARGET_LL}
)