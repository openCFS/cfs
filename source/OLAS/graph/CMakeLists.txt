# Make sure the compiler can find include files
INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR})

SET(GRAPH_SRCS basegraph.cc baseordering.cc sloan.cc graphmanagersimple.cc
	       graphmanagerstdmat.cc graphmanagersbmmat.cc idbcgraph.cc)

ADD_DEFINITIONS()

ADD_LIBRARY(graph-olas STATIC ${GRAPH_SRCS})

TARGET_LINK_LIBRARIES(graph-olas ${METIS_LIBRARY})

