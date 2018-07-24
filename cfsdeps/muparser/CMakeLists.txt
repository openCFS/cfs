CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

PROJECT(muparser CXX)

SET(SRCS
  ./src/muParser.cpp               
  ./src/muParserBase.cpp           
  ./src/muParserBytecode.cpp       
  ./src/muParserCallback.cpp       
  ./src/muParserDLL.cpp            
  ./src/muParserError.cpp          
  ./src/muParserInt.cpp            
  ./src/muParserTest.cpp        
  ./src/muParserTokenReader.cpp 
)

IF(MSVC)
  SET(CMAKE_CXX_FLAGS "/EHsc")
ELSEIF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "-fno-common")
ENDIF()

INCLUDE_DIRECTORIES("./include")

ADD_LIBRARY(muparser STATIC ${SRCS})

# get CXX flags to fix issues with Intel ICC
SET(CMAKE_C_FLAGS "${CFS_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CFS_CXX_FLAGS}")

INSTALL(TARGETS muparser
  ARCHIVE DESTINATION "${LIB_SUFFIX}/${CFS_ARCH_STR}")

INSTALL(FILES
  include/muParser.h
  include/muParserBase.h
  include/muParserBytecode.h
  include/muParserCallback.h
  include/muParserDef.h
  include/muParserDLL.h
  include/muParserError.h
  include/muParserFixes.h
  include/muParserInt.h
  include/muParserStack.h
  include/muParserTemplateMagic.h
  include/muParserTest.h
  include/muParserToken.h
  include/muParserTokenReader.h
  DESTINATION include
)