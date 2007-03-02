# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(GRAPH_SRCS basegraph.cc baseordering.cc sloan.cc graphmanagersimple.cc
	       graphmanagerstdmat.cc graphmanagersbmmat.cc idbcgraph.cc)

ADD_DEFINITIONS()

ADD_LIBRARY(graph-olas STATIC ${GRAPH_SRCS})

TARGET_LINK_LIBRARIES(graph-olas metis)

