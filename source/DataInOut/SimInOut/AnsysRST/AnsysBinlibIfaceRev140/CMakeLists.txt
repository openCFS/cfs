INCLUDE_DIRECTORIES(${ANSYS_V140_INCLUDE_DIR})

# Add the ANSYS ResWrDemo program for test purposes
#ADD_EXECUTABLE(ResWrDemo140 ../ResWrDemo.F ../ResWr.F)
#TARGET_LINK_LIBRARIES(ResWrDemo140
#  ${CFSDEPS_LIBRARY_DIR}/ansys/v140/libbin.so
#  ${GFORTRAN_LIBRARY}
#  )

ADD_EXECUTABLE(ResWrHeadTestCC140
  ../ResWrHeadTest.cc
  ../ResWrTest.F
  ../ResWr.F)
TARGET_LINK_LIBRARIES(ResWrHeadTestCC140
  ${ANSYS_V140_LIBBIN}
  ${CFS_FORTRAN_DYNRT_LIBS}
  -ldl
  ${BOOST_FILESYSTEM_LIB}
  ${BOOST_SYSTEM_LIB}
  -lpthread
  )
SET_TARGET_PROPERTIES(ResWrHeadTestCC140 PROPERTIES LINKER_LANGUAGE CXX 
                                                    LINK_FLAGS -Wl,--allow-shlib-undefined)

#ADD_EXECUTABLE(bintst140 ../bintst.F ../bintrd.F ../bintwr.F)
#TARGET_LINK_LIBRARIES(bintst140 
#  ${CFSDEPS_LIBRARY_DIR}/ansys/v140/libbin.so
#  ${GFORTRAN_LIBRARY}
#  )

ADD_EXECUTABLE(ResRdDemo140
  ../ResRdDemo.F
  ../ResRd.F)
TARGET_LINK_LIBRARIES(ResRdDemo140
  ${ANSYS_V140_LIBBIN}
  ${CFS_FORTRAN_DYNRT_LIBS}
  -ldl
  ${BOOST_FILESYSTEM_LIB}
  ${BOOST_SYSTEM_LIB}
  -lpthread
  )
SET_TARGET_PROPERTIES(ResRdDemo140 PROPERTIES LINKER_LANGUAGE Fortran 
                                              LINK_FLAGS -Wl,--allow-shlib-undefined)


ADD_LIBRARY(CFSAnsysBinlibIfaceRev140 SHARED
  ../AnsysBinlibIfaceGeneral.cc
  ../ResWrTest.F
  ../ResWr.F)

TARGET_LINK_LIBRARIES(CFSAnsysBinlibIfaceRev140
  ${ANSYS_V140_LIBBIN}
  ${CFS_FORTRAN_DYNRT_LIBS}
  )

SET_TARGET_PROPERTIES(CFSAnsysBinlibIfaceRev140 PROPERTIES LINKER_LANGUAGE CXX
                                                           LINK_FLAGS -Wl,--allow-shlib-undefined)
