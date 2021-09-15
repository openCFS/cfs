CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

PROJECT(METIS C)

SET(METIS_SRCS
  Lib/balance.c
  Lib/bucketsort.c
  Lib/ccgraph.c
  Lib/coarsen.c
  Lib/compress.c
  Lib/debug.c
  Lib/estmem.c
  Lib/fm.c
  Lib/fortran.c
  Lib/frename.c
  Lib/graph.c
  Lib/initpart.c
  Lib/kmetis.c
  Lib/kvmetis.c
  Lib/kwayfm.c
  Lib/kwayrefine.c
  Lib/kwayvolfm.c
  Lib/kwayvolrefine.c
  Lib/match.c
  Lib/mbalance.c
  Lib/mbalance2.c
  Lib/mcoarsen.c
  Lib/memory.c
  Lib/mesh.c
  Lib/meshpart.c
  Lib/mfm.c
  Lib/mfm2.c
  Lib/mincover.c
  Lib/minitpart.c
  Lib/minitpart2.c
  Lib/mkmetis.c
  Lib/mkwayfmh.c
  Lib/mkwayrefine.c
  Lib/mmatch.c
  Lib/mmd.c
  Lib/mpmetis.c
  Lib/mrefine.c
  Lib/mrefine2.c
  Lib/mutil.c
  Lib/myqsort.c
  Lib/ometis.c
  Lib/parmetis.c
  Lib/pmetis.c
  Lib/pqueue.c
  Lib/refine.c
  Lib/separator.c
  Lib/sfm.c
  Lib/srefine.c
  Lib/stat.c
  Lib/subdomains.c
  Lib/timing.c
  Lib/util.c
  )

INCLUDE_DIRECTORIES(Lib)

SET(METIS_LIBNAME metis)

IF(WIN32 AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    SET(METIS_LIBNAME metisd)
ENDIF(WIN32 AND CMAKE_BUILD_TYPE STREQUAL "Debug")

ADD_LIBRARY(${METIS_LIBNAME} STATIC ${METIS_SRCS})

INSTALL(TARGETS ${METIS_LIBNAME}
  ARCHIVE DESTINATION "${LIB_SUFFIX}")

INSTALL(FILES
  Lib/defs.h
  Lib/macros.h
  Lib/metis.h
  Lib/proto.h
  Lib/rename.h
  Lib/struct.h DESTINATION include)