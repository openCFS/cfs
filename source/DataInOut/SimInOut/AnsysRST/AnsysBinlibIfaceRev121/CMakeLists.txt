INCLUDE_DIRECTORIES(${ANSYS_V121_INCLUDE_DIR})

# Add the ANSYS ResWrDemo program for test purposes
#ADD_EXECUTABLE(ResWrDemo121 ../ResWrDemo.F ../ResWr.F)
#TARGET_LINK_LIBRARIES(ResWrDemo121
#  ${CFSDEPS_LIBRARY_DIR}/ansys/v121/libbin.so
#  ${GFORTRAN_LIBRARY}
#  )

ADD_EXECUTABLE(ResWrHeadTestCC121
  ../ResWrHeadTest.cc
  ../ResWrTest.F
  ../ResWr.F)

SET(TARGET_LL
  ${ANSYS_V121_LIBBIN}
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

TARGET_LINK_LIBRARIES(ResWrHeadTestCC121 ${TARGET_LL})

SET_TARGET_PROPERTIES(ResWrHeadTestCC121 PROPERTIES LINKER_LANGUAGE CXX 
                                                    LINK_FLAGS -Wl,--allow-shlib-undefined)

#ADD_EXECUTABLE(bintst121 ../bintst.F ../bintrd.F ../bintwr.F)
#TARGET_LINK_LIBRARIES(bintst121 
#  ${CFSDEPS_LIBRARY_DIR}/ansys/v121/libbin.so
#  ${GFORTRAN_LIBRARY}
#  )


ADD_LIBRARY(CFSAnsysBinlibIfaceRev121 SHARED
  ../AnsysBinlibIfaceGeneral.cc
  ../ResWrTest.F
  ../ResWr.F)

SET(TARGET_LL
  CFSAnsysBinlibIface
  ${ANSYS_V121_LIBBIN}
  ${CFS_FORTRAN_DYNRT_LIBS}
)

IF(MINGW)
  LIST(APPEND TARGET_LL
   logging
   matvec
   febasis
   cfsgeneral
  )
ENDIF()

TARGET_LINK_LIBRARIES(CFSAnsysBinlibIfaceRev121
  ${TARGET_LL}
  )

SET_TARGET_PROPERTIES(CFSAnsysBinlibIfaceRev121 PROPERTIES LINKER_LANGUAGE CXX
                                                           LINK_FLAGS -Wl,--allow-shlib-undefined)
