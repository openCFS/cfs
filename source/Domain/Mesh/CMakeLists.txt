if(USE_XERCES)
  include_directories(${XERCES_INCLUDE_DIR})
endif()

set(MESH_SRCS 
  Grid.cc
  NcInterfaces/BaseNcInterface.cc
  NcInterfaces/MortarInterface.cc
  GridCFS/GridCFS.cc
  MeshUtils/EntityAssociation.cc
  MeshUtils/Intersection/VolumeGridIntersection.cc
  MeshUtils/Intersection/IntersectAlgos/CoordTetra.cc
  MeshUtils/Intersection/IntersectAlgos/TriaIntersect.cc
  MeshUtils/Intersection/IntersectAlgos/ElemIntersect2D.cc)


add_library(mesh STATIC ${MESH_SRCS})

set(TARGET_LL
  utils
  elemmapping
  coordsystems
  datainout
  paramh
  matvec
  cfsgeneral
  integration)

if(USE_CGAL)
  if(UNIX)
    list(APPEND TARGET_LL ${MPFR_LIBRARY} ${GMP_LIBRARY})
  endif()
  list(APPEND TARGET_LL ${BOOST_THREAD_LIB})				
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel" AND USE_OPENMP)
  list(APPEND TARGET_LL ${BOOST_THREAD_LIB})
endif()

if(USE_XERCES)
  list(APPEND TARGET_LL ${XERCES_LIBRARY})
endif()

list(REMOVE_DUPLICATES TARGET_LL)

target_link_libraries(mesh ${TARGET_LL})
