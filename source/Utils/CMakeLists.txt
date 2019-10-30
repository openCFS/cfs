SET(UTILS_SRCS
  ApproxData.cc
  BiLinInterpolate.cc
  BiCubicInterpolate.cc
  BSpline.cc
  Coil.cc
  Interpolate1D.cc
  LinInterpolate.cc
  mathfunctions.cc
  Point.cc
  Polynomial.cc
  SmoothSpline.cc
  Timer.cc
  tools.cc
  TriLinInterpolate.cc
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
  # In order to make core CFS++ independent from the ANSYS version,
  # we do not want to link against ANSYS DLLs directly anymore. Instead
  # we provide interface DLLs for each ANSYS release, which will be
  # loaded on demand at runtime. This can even be triggered from the XML file:
  #
  #       <rst id="rst" ansysRevision="14.0"/>
  #
  # The interface for loading C++ objects from DLLs / shared object files
  # has been implemented following:
  #
  # http://www.drdobbs.com/dynamically-loaded-c-objects/184401900
  # http://www.codeproject.com/Articles/28969/HowTo-Export-C-classes-from-a-DLL#CppMatureApproach
  ADD_LIBRARY(CFSDynLd SHARED DynamicLibrary.cc)

  IF(NOT CMAKE_SYSTEM_NAME STREQUAL "Windows")
    TARGET_LINK_LIBRARIES(CFSDynLd -ldl )
  ENDIF(NOT CMAKE_SYSTEM_NAME STREQUAL "Windows")

  ADD_LIBRARY(DynamicLibTestClassIface SHARED DynamicLibTestClassIface.cc)
  TARGET_LINK_LIBRARIES(DynamicLibTestClassIface CFSDynLd )

  ADD_LIBRARY(DynamicLibTestClass SHARED DynamicLibTestClass.cc)
  TARGET_LINK_LIBRARIES(DynamicLibTestClass DynamicLibTestClassIface )

  ADD_EXECUTABLE(DynamicLibTest DynamicLibTest.cc)
  TARGET_LINK_LIBRARIES(DynamicLibTest DynamicLibTestClassIface )
ENDIF(USE_ANSYSRST)

IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
  set_source_files_properties(piezoMicroModelBK.cc PROPERTIES COMPILE_FLAGS -wd279)
ENDIF()
