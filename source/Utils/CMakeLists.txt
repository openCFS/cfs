SET(UTILS_SRCS
  ApproxData.cc
  LinInterpolate.cc
  BiLinInterpolate.cc
  TriLinInterpolate.cc
  CubicInterpolate.cc
  BiCubicInterpolate.cc
  TriCubicInterpolate.cc
  BSpline.cc
  Coil.cc
  Interpolate1D.cc
  mathfunctions.cc
  Point.cc
  Polynomial.cc
  SmoothSpline.cc
  Timer.cc
  tools.cc)

ADD_LIBRARY(utils STATIC ${UTILS_SRCS})

SET(TARGET_LL 
  datainout
  matvec
  coordsystems
  ${BOOST_SERIALIZATION_LIB}
  ${BOOST_REGEX_LIB}
)

TARGET_LINK_LIBRARIES(utils 
  ${TARGET_LL}
)

IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
  set_source_files_properties(piezoMicroModelBK.cc PROPERTIES COMPILE_FLAGS -wd279)
ENDIF()
