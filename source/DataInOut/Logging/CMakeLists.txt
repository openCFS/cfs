if(DEBUG)
  set(LOGGING_SRCS LogConfigurator.cc)
else()
  # release mode we need no logging but have dummy code to avoid linker warnings/errors
  set(LOGGING_SRCS LogConfiguratorDummy.cc)
endif()

SET(TARGET_LL
  paramh
  utils
  ${BOOST_LIBRARY} 
  ${ZLIB_LIBRARY} )

ADD_LIBRARY(logging STATIC ${LOGGING_SRCS})

ADD_DEPENDENCIES(logging boost)

TARGET_LINK_LIBRARIES(logging ${TARGET_LL} ${BOOST_LIBRARY} ${ZLIB_LIBRARY} )
