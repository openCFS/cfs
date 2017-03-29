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
  Filters/Interpolators/CentroidInterpolator.cc
  Filters/Interpolators/NearestNeighbourInterpolator.cc
  Filters/Interpolators/Cell2NodeInterpolator.cc
  Filters/Interpolators/FEBasedInterpolator.cc
  Filters/Interpolators/GridIntersectionFilter.cc
  Filters/Interpolators/RBFInterpolator.cc
  Filters/Interpolators/Node2CellInterpolator.cc
  Filters/SignalProcessing/FftFilter.cc
  Filters/SignalProcessing/TimeMeanFilter.cc
  Utils/EqnNumberingSimple.cc
  Utils/CFSDatProgramOptions.cc
  Utils/ResultManager.cc
#  Utils/KNNSearch.cc
)

# Main CFS++ Dat executable target.
ADD_EXECUTABLE(cfsdatbin ${CFSDAT_SRCS})

SET(TARGET_LL 
  datainout
  matvec
)

IF(MKL_BLAS_LIB)
  LIST(APPEND TARGET_LL ${MKL_BLAS_LIB})
ENDIF(MKL_BLAS_LIB)

TARGET_LINK_LIBRARIES(cfsdatbin ${TARGET_LL})

SET(CFS_LINK_FLAGS "")

SET_TARGET_PROPERTIES(cfsdatbin PROPERTIES LINK_FLAGS "${CFS_LINK_FLAGS}")

# Install the cfstool executable target to ${CMAKE_INSTALL_PREFIX}/bin
#INSTALL_TARGETS(/bin cfsdat)
