# Make sure the compiler can find include files
INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR})

SET(GRAPH_SRCS
    BaseGraph.cc
    BaseOrdering.cc
    GraphManager.cc
    IDBC_Graph.cc
    Sloan.cc
	)

ADD_DEFINITIONS()

ADD_LIBRARY(graph-olas STATIC ${GRAPH_SRCS})

TARGET_LINK_LIBRARIES(graph-olas ${METIS_LIBRARY})

ADD_DEPENDENCIES(graph-olas boost)
