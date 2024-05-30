set(OPTIMIZER_SRCS
  BaseOptimizer.cc
  ShapeOptimizer.cc
  OptimalityCondition.cc
  EvaluateOnly.cc
  GradientCheck.cc
  MMA.cc
  BFGS.cc)

set(TARGET_LL
  optimization
  cfsgeneral
  datainout
  design
  domain
  driver
  logging
  matvec
  paramh
  pde
  utils)

# we need it almost always ?!
set(DEPENDS
  "${LAPACK_LIBRARY}"
  "${BLAS_LIBRARY}"
  "${CFS_FORTRAN_LIBS}"
  )


if(USE_SCPIP)
   set(TARGET_LL ${TARGET_LL} ${SCPIP_LIBRARY} ${DEPENDS})
   # message("${CFS_FORTRAN_LIBS}")
   # SCPIPBase is also the open source project C++SCPIP.
   # With this definition the file lives in both worlds! 
   add_definitions(-DUSE_4_CFS)
   
   set(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} SCPIP.cc SCPIPBase.cc)
endif()


if(USE_SNOPT)
   set(TARGET_LL ${TARGET_LL}  ${SNOPT_LIBRARY} ${DEPENDS})
   set(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} SnOpt.cc)
   set(OPTIMIZER_FSRCS SnOptInterface.f)
   ADD_LIBRARY(snopt_ftnif STATIC ${OPTIMIZER_FSRCS})
endif()

if(USE_IPOPT)
   set(TARGET_LL ${TARGET_LL}  ${IPOPT_LIBRARY} ${DEPENDS})
   include_directories(${IPOPT_INCLUDE_DIR})
   set(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} IPOPT.cc FeasPP.cc FeasSubProblem.cc)
endif()

if(USE_EMBEDDED_PYTHON)
  # see FindPrograms.cmake on how the PYTHON_ stuff is set
  set(TARGET_LL ${TARGET_LL} ${PYTHON_LIBRARY})
  # required for numpy
  include_directories(${PYTHON_SITE_PACKAGES_DIR})
  set(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} PythonOptimizer.cc)
endif()

if(USE_SGP)
  # we always need MKL for building SGP
  set(TARGET_LL ${TARGET_LL} ${SGP_LIBRARY} ${DEPENDS} ${MKL_BLAS_LIB})
  include_directories(${MKL_INCLUDE_DIR})
  set(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} SGP.cc SGPHolder.cc)
endif()

add_library(optimizer STATIC ${OPTIMIZER_SRCS})

IF(USE_SNOPT)
  TARGET_LINK_LIBRARIES(optimizer snopt_ftnif ${TARGET_LL})
ELSE()
  TARGET_LINK_LIBRARIES(optimizer ${TARGET_LL})
ENDIF()
