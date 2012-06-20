INCLUDE_DIRECTORIES(${ANSYS_V120_INCLUDE_DIR})

# Add the ANSYS ResWrDemo program for test purposes
#ADD_EXECUTABLE(ResWrDemo120 ../ResWrDemo.F ../ResWr.F)
#TARGET_LINK_LIBRARIES(ResWrDemo120
#  ${CFSDEPS_LIBRARY_DIR}/ansys/v120/libbin.so
#  ${GFORTRAN_LIBRARY}
#  )

ADD_EXECUTABLE(ResWrHeadTestCC120
  ../ResWrHeadTest.cc
  ../ResWrTest.F
  ../ResWr.F)

SET(TARGET_LL
  ${ANSYS_V120_LIBBIN}
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

TARGET_LINK_LIBRARIES(ResWrHeadTestCC120
  ${TARGET_LL}
  )
SET_TARGET_PROPERTIES(ResWrHeadTestCC120 PROPERTIES LINKER_LANGUAGE CXX 
                                                    LINK_FLAGS -Wl,--allow-shlib-undefined)

#ADD_EXECUTABLE(bintst120 ../bintst.F ../bintrd.F ../bintwr.F)
#TARGET_LINK_LIBRARIES(bintst120 
#  ${CFSDEPS_LIBRARY_DIR}/ansys/v120/libbin.so
#  ${GFORTRAN_LIBRARY}
#  )


ADD_LIBRARY(CFSAnsysBinlibIfaceRev120 SHARED
  ../AnsysBinlibIfaceGeneral.cc
  ../ResWrTest.F
  ../ResWr.F)

TARGET_LINK_LIBRARIES(CFSAnsysBinlibIfaceRev120
  ${ANSYS_V120_LIBBIN}
  ${CFS_FORTRAN_DYNRT_LIBS}
  )

SET_TARGET_PROPERTIES(CFSAnsysBinlibIfaceRev120 PROPERTIES LINKER_LANGUAGE CXX
                                                           LINK_FLAGS -Wl,--allow-shlib-undefined)
