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

if(USE_EMBEDDED_PYTHON)
  list(APPEND TARGET_LL ${Python_LIBRARIES})
  set(UTILS_SRCS ${UTILS_SRCS} PythonKernel.cc PythonKernelFunctions.cc)
endif()

ADD_LIBRARY(utils STATIC ${UTILS_SRCS})

SET(TARGET_LL 
  datainout
  matvec
  coordsystems
  ${BOOST_LIBRARY} 
  ${ZLIB_LIBRARY} )

TARGET_LINK_LIBRARIES(utils 
  ${TARGET_LL}
)

IF(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set_source_files_properties(piezoMicroModelBK.cc PROPERTIES COMPILE_FLAGS -wd279)
ENDIF()
