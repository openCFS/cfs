INCLUDE_DIRECTORIES(${ANSYS_V110_INCLUDE_DIR})

# Add the ANSYS ResWrDemo program for test purposes
#ADD_EXECUTABLE(ResWrDemo110 ../ResWrDemo.F ../ResWr.F)
#TARGET_LINK_LIBRARIES(ResWrDemo110
#  ${CFSDEPS_LIBRARY_DIR}/ansys/v110/libbin.so
#  ${GFORTRAN_LIBRARY}
#  )

ADD_EXECUTABLE(ResWrHeadTestCC110
  ../ResWrHeadTest.cc
  ../ResWrTest.F
  ../ResWr.F)

SET(TARGET_LL
  ${ANSYS_V110_LIBBIN}
  ${CFS_FORTRAN_DYNRT_LIBS}
  ${BOOST_FILESYSTEM_LIB}
  ${BOOST_SYSTEM_LIB}  
  )

IF(NOT CMAKE_SYSTEM_NAME STREQUAL "Windows")
  LIST(APPEND TARGET_LL
    -ldl
    -lpthread
    )
ENDIF()

TARGET_LINK_LIBRARIES(ResWrHeadTestCC110 ${TARGET_LL})

SET_TARGET_PROPERTIES(ResWrHeadTestCC110 PROPERTIES LINKER_LANGUAGE CXX 
                                                    LINK_FLAGS -Wl,--allow-shlib-undefined)

#ADD_EXECUTABLE(bintst110 ../bintst.F ../bintrd.F ../bintwr.F)
#TARGET_LINK_LIBRARIES(bintst110 
#  ${CFSDEPS_LIBRARY_DIR}/ansys/v110/libbin.so
#  ${GFORTRAN_LIBRARY}
#  )


ADD_LIBRARY(CFSAnsysBinlibIfaceRev110 SHARED
  ../AnsysBinlibIfaceGeneral.cc
  ../ResWrTest.F
  ../ResWr.F)

SET(TARGET_LL
  CFSAnsysBinlibIface
  ${ANSYS_V110_LIBBIN}
  ${CFS_FORTRAN_DYNRT_LIBS}
)

IF(MINGW)
  LIST(APPEND TARGET_LL
   logging
   matvec
   febasis
   cfsgeneral
   
#  
#  datainout
#  cfsgeneral
#  utils
#  definefiles
#  paramh
#  utils
#  domain
#  elemmapping
#  utils
  )
ENDIF()

TARGET_LINK_LIBRARIES(CFSAnsysBinlibIfaceRev110
  ${TARGET_LL}
  )

SET_TARGET_PROPERTIES(CFSAnsysBinlibIfaceRev110 PROPERTIES LINKER_LANGUAGE CXX
                                                           LINK_FLAGS -Wl,--allow-shlib-undefined)
