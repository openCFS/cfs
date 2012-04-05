INCLUDE_DIRECTORIES(${ANSYS_V130_INCLUDE_DIR})

# Add the ANSYS ResWrDemo program for test purposes
#ADD_EXECUTABLE(ResWrDemo130 ../ResWrDemo.F ../ResWr.F)
#TARGET_LINK_LIBRARIES(ResWrDemo130
#  ${CFSDEPS_LIBRARY_DIR}/ansys/v130/libbin.so
#  ${GFORTRAN_LIBRARY}
#  )

ADD_EXECUTABLE(ResWrHeadTestCC130
  ../ResWrHeadTest.cc
  ../ResWrTest.F
  ../ResWr.F)

SET(TARGET_LL
  ${ANSYS_V130_LIBBIN}
  ${CFS_FORTRAN_DYNRT_LIBS}
  -ldl
  ${BOOST_FILESYSTEM_LIB}
  ${BOOST_SYSTEM_LIB}
  -lpthread
)
IF(NOT CFS_FORTRAN_COMPILER_NAME STREQUAL "IFORT")
  SET(TARGET_LL
    ${TARGET_LL}
    -lirc
  )
ENDIF(NOT CFS_FORTRAN_COMPILER_NAME STREQUAL "IFORT")

TARGET_LINK_LIBRARIES(ResWrHeadTestCC130
  ${TARGET_LL}
  )
SET_TARGET_PROPERTIES(ResWrHeadTestCC130 PROPERTIES LINKER_LANGUAGE CXX 
                                                    LINK_FLAGS -Wl,--allow-shlib-undefined)

#ADD_EXECUTABLE(bintst130 ../bintst.F ../bintrd.F ../bintwr.F)
#TARGET_LINK_LIBRARIES(bintst130 
#  ${CFSDEPS_LIBRARY_DIR}/ansys/v130/libbin.so
#  ${GFORTRAN_LIBRARY}
#  )

ADD_EXECUTABLE(ResRdDemo130
  ../ResRdDemo.F
  ../ResRd.F)
TARGET_LINK_LIBRARIES(ResRdDemo130
  ${ANSYS_V130_LIBBIN}
  ${CFS_FORTRAN_DYNRT_LIBS}
  -ldl
  ${BOOST_FILESYSTEM_LIB}
  ${BOOST_SYSTEM_LIB}
  -lpthread
  )
SET_TARGET_PROPERTIES(ResRdDemo130 PROPERTIES LINKER_LANGUAGE Fortran 
                                              LINK_FLAGS -Wl,--allow-shlib-undefined)


ADD_LIBRARY(CFSAnsysBinlibIfaceRev130 SHARED
  ../AnsysBinlibIfaceGeneral.cc
  ../ResWrTest.F
  ../ResWr.F)

TARGET_LINK_LIBRARIES(CFSAnsysBinlibIfaceRev130
  ${ANSYS_V130_LIBBIN}
  ${CFS_FORTRAN_DYNRT_LIBS}
  )

SET_TARGET_PROPERTIES(CFSAnsysBinlibIfaceRev130 PROPERTIES LINKER_LANGUAGE CXX
                                                           LINK_FLAGS -Wl,--allow-shlib-undefined)
