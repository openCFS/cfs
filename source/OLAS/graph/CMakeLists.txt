# Make sure the compiler can find include files
INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIR})

SET(GRAPH_SRCS
  graphmanagersbmmat.cc
  graphmanagerstdmat.cc
  basegraph.cc
  graphmanagersimple.cc
  baseordering.cc
  sloan.cc
  idbcgraph.cc
	)

ADD_DEFINITIONS()

ADD_LIBRARY(graph-olas STATIC ${GRAPH_SRCS})

TARGET_LINK_LIBRARIES(graph-olas 
  ${METIS_LIBRARY}
  cfsgeneral
  matvec
  utils
)

