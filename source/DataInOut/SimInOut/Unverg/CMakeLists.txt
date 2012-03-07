SET(SIMOUTUNV_SRCS
  SimOutputUnv.cc
  SimInputUnv.cc
  unv.cc
  unv_dat.cc
  unv_if.cc
  unv_reader.cc
)

ADD_LIBRARY(simoutunv STATIC ${SIMOUTUNV_SRCS})

IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
set_source_files_properties(SimOutputUnv.cc
                            PROPERTIES COMPILE_FLAGS
                            "-D '__sync_fetch_and_add(ptr,addend)=_InterlockedExchangeAdd(const_cast<void*>(reinterpret_cast<volatile void*>(ptr)), addend)'")
ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")

TARGET_LINK_LIBRARIES(simoutunv
  cfsgeneral
  datainout
  domain
  elements
  logging
  matvec
  paramh
  pde
  utils
)
