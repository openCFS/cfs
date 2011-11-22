SET(UTILS_SRCS
  ApproxData.cc
  Coil.cc
  DampLayer.cc
  Interpolate1D.cc
  LinInterpolate.cc
  mathfunctions.cc
  Timer.cc
  SmoothSpline.cc
  tools.cc
  )

ADD_LIBRARY(utils STATIC ${UTILS_SRCS})

SET(TARGET_LL 
  datainout
  matvec
  coordsystems
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
