if(USE_PHIST_EV AND USE_PHIST_CG)
SET(PHIST_SRCS PhistEigenSolver.cc PhistLinearSolver.cc PhistCore.cc)
elseif(USE_PHIST_EV)
SET(PHIST_SRCS PhistEigenSolver.cc PhistCore.cc)
elseif(USE_PHIST_CG)
SET(PHIST_SRCS PhistLinearSolver.cc PhistCore.cc)
endif()

ADD_LIBRARY(phist-olas STATIC ${PHIST_SRCS})

SET(TARGET_LL
  ${PHIST_LIBRARY} 
  ${GHOST_LIBRARY}
  ${HWLOC_LIBRARY} )

TARGET_LINK_LIBRARIES(phist-olas ${TARGET_LL})

