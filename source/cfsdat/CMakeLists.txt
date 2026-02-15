INCLUDE_DIRECTORIES(
  ${HDF5_INCLUDE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR})

IF(MKL_INCLUDE_DIR)
  INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR})
ENDIF(MKL_INCLUDE_DIR)

# List of source codes for current target.
SET(CFSDAT_SRCS
  cfsdat.cc
  Filters/BaseFilter.cc
  Filters/BaseMeshFilterType.cc
  Filters/MeshFilter.cc
  Filters/Arithmetic/BinOpFilter.cc
  Filters/Arithmetic/TensorFilter.cc
  Filters/Input/InputFilter.cc
  Filters/Output/OutputFilter.cc
  Filters/Derivatives/TimeDerivFilter.cc
  Filters/Derivatives/RotatingSubstDt.cc
  Filters/Derivatives/GradientDifferentiator.cc
  Filters/Derivatives/DivergenceDifferentiator.cc
  Filters/Derivatives/CurlDifferentiator.cc
  Filters/Derivatives/Lighthill.cc
  Filters/Derivatives/AeroacousticBase.cc
  Filters/Derivatives/PostLighthillSource.cc
  Filters/Derivatives/GaussDerivative.cc
  Filters/Interpolators/AbstractInterpolator.cc
  Filters/Interpolators/CentroidInterpolator.cc
  Filters/Interpolators/NearestNeighbourInterpolator.cc
  Filters/Interpolators/Cell2NodeInterpolator.cc
  Filters/Interpolators/FEBasedInterpolator.cc
  Filters/Interpolators/GridIntersectionFilter.cc
  Filters/Interpolators/RBFInterpolator.cc
  Filters/Interpolators/Node2CellInterpolator.cc
  Filters/Interpolators/L2norm.cc
  Filters/SignalProcessing/FftFilter.cc
  Filters/SignalProcessing/FIRFilter.cc
  Filters/SignalProcessing/TemporalBlendingFilter.cc
  Filters/SignalProcessing/TimeMeanFilter.cc
  Filters/SynteticSources/BaseSNGR.cc
  Filters/Utils/VolumeMultiplication.cc
  DatUtils/KNNIndexSearch.cc
  DatUtils/EqnNumberingSimple.cc
  DatUtils/CFSDatProgramOptions.cc
  DatUtils/ResultManager.cc
  DatUtils/ResultCache.cc
#  Utils/KNNSearch.cc
)

# Main cfsDAT executable target.
ADD_EXECUTABLE(cfsdat ${CFSDAT_SRCS})

SET(TARGET_LL 
  datainout
  matvec
)

IF(MKL_BLAS_LIB)
  LIST(APPEND TARGET_LL ${MKL_BLAS_LIB})
ENDIF(MKL_BLAS_LIB)

if(USE_BLAS_LAPACK STREQUAL "NETLIB")
  assert_set(LAPACK_LIBRARY)
  list(APPEND TARGET_LL ${LAPACK_LIBRARY})
  list(APPEND TARGET_LL ${BLAS_LIBRARY})
endif()

IF(DEPS_SEQUENTIAL)
  LIST(APPEND TARGET_LL ${DEPS_SEQUENTIAL})
ENDIF(DEPS_SEQUENTIAL)

TARGET_LINK_LIBRARIES(cfsdat ${TARGET_LL})

# IntelLLVM und RHEL8 chooses a too old lib for std::filesystem. CFS_LINKER_FLAGS does not work
if(UNIX AND CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  target_link_libraries(cfsdat stdc++fs)
endif()


SET_TARGET_PROPERTIES(cfsdat PROPERTIES LINK_FLAGS "${CFS_LINK_FLAGS}")
