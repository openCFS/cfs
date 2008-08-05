INCLUDE_DIRECTORIES(${ANSYS_ROOT_DIR}/ansys/customize/include)

IF(CFS_ARCH STREQUAL "I386")
  SET(CMAKE_Fortran_FLAGS "-safe-cray-ptr -fcray-pointer -DLINUXIA32_SYS")
ENDIF(CFS_ARCH STREQUAL "I386")

IF(CFS_ARCH STREQUAL "X86_64")
  SET(CMAKE_Fortran_FLAGS "-safe-cray-ptr -fcray-pointer -DLINUXOP64_SYS")
ENDIF(CFS_ARCH STREQUAL "X86_64")

IF(DEBUG)
  SET(CMAKE_Fortran_FLAGS "-g ${CMAKE_Fortran_FLAGS}")
ENDIF(DEBUG)

SET(SIMOUTANSYSRST_SRCS simOutputRST.cc ResWr.F)

ADD_LIBRARY(simoutansysrst STATIC ${SIMOUTANSYSRST_SRCS})

SET_TARGET_PROPERTIES(simoutansysrst PROPERTIES LINKER_LANGUAGE CXX)

IF(CFS_ARCH STREQUAL "I386")
  TARGET_LINK_LIBRARIES(simoutansysrst
    ${ANSYS_ROOT_DIR}/ansys/customize/misc/linia32/libbin.so
    ${ANSYS_ROOT_DIR}/ansys/syslib/linia32/libcxa.so.5
    ${ANSYS_ROOT_DIR}/ansys/syslib/linia32/libunwind.so.5
    -ldl
    ${GFORTRAN_LIBRARY})
ENDIF(CFS_ARCH STREQUAL "I386")

IF(CFS_ARCH STREQUAL "X86_64")
  TARGET_LINK_LIBRARIES(simoutansysrst
    ${ANSYS_ROOT_DIR}/ansys/customize/misc/linop64/libbin.so
    ${ANSYS_ROOT_DIR}/ansys/syslib/linop64/libpgf90.so
    ${GFORTRAN_LIBRARY})
ENDIF(CFS_ARCH STREQUAL "X86_64")

# Add the ANSYS ResWrDemo program for test purposes
ADD_EXECUTABLE(ResWrDemo ResWrDemo.F)
TARGET_LINK_LIBRARIES(ResWrDemo simoutansysrst)

ADD_EXECUTABLE(ResWrHeadTestCC ResWrHeadTest.cc ResWrTest.F)
TARGET_LINK_LIBRARIES(ResWrHeadTestCC
  simoutansysrst
  ${BOOST_FILESYSTEM_LIB}
  )

ADD_EXECUTABLE(bintst bintst.F bintrd.F bintwr.F)
TARGET_LINK_LIBRARIES(bintst simoutansysrst)
