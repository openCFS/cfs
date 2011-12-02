INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})

SET(SIMINOUTHDF5_SRCS
  hdf5io.cc
  SimInputHDF5.cc
  SimOutputHDF5.cc
)


ADD_LIBRARY(siminouthdf5 STATIC ${SIMINOUTHDF5_SRCS})

TARGET_LINK_LIBRARIES(siminouthdf5
  mesh
  ${HDF5_CPP_LIBRARY}
  ${HDF5_LIBRARY}
  ${ZLIB_LIBRARY}
  ${BOOST_DATE_TIME_LIB}
  )

ADD_EXECUTABLE(h5tool h5tool.cc)

SET(TARGET_LL siminouthdf5 cfsgeneral)

# Fix link problems with Intel 9 compiler
IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC" AND
  CFS_CXX_COMPILER_VER MATCHES "9\\.")

SET(TARGET_LL ${TARGET_LL} "-Wl,-Bstatic -lstdc++")

ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "ICC" AND
     CFS_CXX_COMPILER_VER MATCHES "9\\.")

TARGET_LINK_LIBRARIES(h5tool ${TARGET_LL})

# Only build the Boost HDF5 serialization test executable with GCC compiler
# since it is not absolutely necessary anyway and some Intel compilers complain
# about it.
IF(FALSE AND CFS_CXX_COMPILER_NAME STREQUAL "GCC")
  ADD_EXECUTABLE(h5SerializeTest serializeTest.cc)

  TARGET_LINK_LIBRARIES(h5SerializeTest
    ${BOOST_SERIALIZATION_HDF5_LIB}
    ${BOOST_SERIALIZATION_LIB}
    ${HDF5_CPP_LIBRARY}
    ${HDF5_LIBRARY}
  )
ENDIF(FALSE AND CFS_CXX_COMPILER_NAME STREQUAL "GCC")

