SET(SIMOUTUNV_SRCS
  SimOutputUnv.cc
  SimInputUnv.cc
  unv_reader.cc
  unv.cc
  unv_dat.cc
  unv_if.cc
)

ADD_LIBRARY(simoutunv STATIC ${SIMOUTUNV_SRCS})

ADD_DEPENDENCIES(simoutunv boost)

IF(TARGET cgal)
  ADD_DEPENDENCIES(simoutunv cgal)
ENDIF()

if(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
  if(NOT WIN32)
    set_source_files_properties(SimOutputUnv.cc PROPERTIES COMPILE_FLAGS
         "-D '__sync_fetch_and_add(ptr,addend)=_InterlockedExchangeAdd(const_cast<void*>(reinterpret_cast<volatile void*>(ptr)), addend)'")
  endif()
endif()

