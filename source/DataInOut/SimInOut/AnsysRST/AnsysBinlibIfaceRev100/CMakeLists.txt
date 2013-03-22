INCLUDE_DIRECTORIES(${ANSYS_V100_INCLUDE_DIR})

# Add the ANSYS ResWrDemo program for test purposes
#ADD_EXECUTABLE(ResWrDemo100 ../ResWrDemo.F ../ResWr.F)
#TARGET_LINK_LIBRARIES(ResWrDemo100
#  ${CFSDEPS_LIBRARY_DIR}/ansys/v100/libbin.so
#  ${GFORTRAN_LIBRARY}
#  )

ADD_EXECUTABLE(ResWrHeadTestCC100 ../ResWrHeadTest.cc ../ResWrTest.F ../ResWr.F)

SET(TARGET_LL
  ${ANSYS_V100_LIBBIN}
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

TARGET_LINK_LIBRARIES(ResWrHeadTestCC100 ${TARGET_LL})

SET_TARGET_PROPERTIES(ResWrHeadTestCC100 PROPERTIES LINKER_LANGUAGE CXX 
                                                    LINK_FLAGS -Wl,--allow-shlib-undefined)

#ADD_EXECUTABLE(bintst100 ../bintst.F ../bintrd.F ../bintwr.F)
#TARGET_LINK_LIBRARIES(bintst100 
#  ${CFSDEPS_LIBRARY_DIR}/ansys/v100/libbin.so
#  ${GFORTRAN_LIBRARY}
#  )


ADD_LIBRARY(CFSAnsysBinlibIfaceRev100 SHARED
  ../AnsysBinlibIfaceGeneral.cc
  ../ResWrTest.F
  ../ResWr.F)

TARGET_LINK_LIBRARIES(CFSAnsysBinlibIfaceRev100
  ${ANSYS_V100_LIBBIN}
  ${CFS_FORTRAN_DYNRT_LIBS}
  )

SET_TARGET_PROPERTIES(CFSAnsysBinlibIfaceRev100 PROPERTIES LINKER_LANGUAGE CXX
                                                           LINK_FLAGS -Wl,--allow-shlib-undefined)
