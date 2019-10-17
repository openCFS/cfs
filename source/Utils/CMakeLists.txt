SET(UTILS_SRCS
  ApproxData.cc
  Coil.cc
  Interpolate1D.cc
  LinInterpolate.cc
  BiLinInterpolate.cc
  TriLinInterpolate.cc
  mathfunctions.cc
  Polynomial.cc
  Timer.cc
  SmoothSpline.cc
  tools.cc
  Point.cc
  )

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
