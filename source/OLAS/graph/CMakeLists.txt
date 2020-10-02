# Make sure the compiler can find include files
if(USE_METIS)
  INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR})
endif()

SET(GRAPH_SRCS
    BaseGraph.cc
    BaseOrdering.cc
    GraphManager.cc
    IDBC_Graph.cc
    Sloan.cc
	)

ADD_DEFINITIONS()

ADD_LIBRARY(graph-olas STATIC ${GRAPH_SRCS})

if(USE_METIS)
  TARGET_LINK_LIBRARIES(graph-olas ${METIS_LIBRARY})
endif()

ADD_DEPENDENCIES(graph-olas boost)
