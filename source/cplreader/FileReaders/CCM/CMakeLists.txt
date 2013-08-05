SET(CCM_SRCS
  FileReader_CCM.cc
  
  MeshSplitter.cc
  MeshPreparator.cc
  FieldDataAcceptor.cc
  DataFile.cc
  MeshFile.cc
  FieldDataAcceptor.cc
  Utilities.cc
  NameFieldDataAcceptor.cc
  PreAllocatedDataAcceptor.cc
  Meshes.cc
  DualMesher.cc
  ValueMapCatcher.cc
  FlowDataLoader.cc
  divUtil.cpp
  )

ADD_LIBRARY(cplreader_ccm STATIC ${CCM_SRCS})

TARGET_LINK_LIBRARIES(cplreader_ccm ${CCMIO_LIBRARY})

