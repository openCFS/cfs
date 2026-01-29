SET(GENERAL_SRCS Environment.cc Exception.cc)

ADD_LIBRARY(cfsgeneral STATIC ${GENERAL_SRCS})

if(USE_MKL)
  target_include_directories(cfsgeneral PUBLIC ${MKL_INCLUDE_DIR})
endif()

ADD_DEPENDENCIES(cfsgeneral boost)
