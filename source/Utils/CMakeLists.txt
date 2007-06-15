# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(UTILS_SRCS myclock.cc nodestoresol.cc elemstoresol.cc tools.cc vector.cc
  tracing.cc Coil.cc ApproxData.cc SmoothSpline.cc mathfunctions.cc
  hysteresis.cc preisach.cc LinInterpolate.cc jiles.cc profiler.cc
  coordSystem.cc cylCoordSys.cc defaultCoordSys.cc dampLayer.cc
  interpolate.cc result.cc polCoordSys.cc readCFDbinary.cc VecStat.cc)

ADD_LIBRARY(utils STATIC ${UTILS_SRCS})


