
IF(USE_IPOPT)
  INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
  LINK_DIRECTORIES(${IPOPT_LINK_DIR})
ENDIF(USE_IPOPT)

IF(USE_ADOLC)
  INCLUDE_DIRECTORIES(${ADOLC_INCLUDE_DIR})
  #LINK_DIRECTORIES(${ADOLC_LIBRARY})
ENDIF(USE_ADOLC)




# List of source codes for current target.
SET(CFS_SRCS CFS.cc)

# Main CFS++ executable target.
ADD_EXECUTABLE(cfsbin ${CFS_SRCS})

# Define a list of libraries CFS++ depends upon. At the moment there
# are some nasty circular dependencies which should be resolved in
# the future (see e.g. pde <-> domain, domain <-> gridcfs !!!)
  
IF(USE_ADOLC)
  SET(TARGET_LL
    cfsgeneral
    utils
    datainout
    paramh
    utils
    domain
    elemmapping
    ${CFS_FORTRAN_DYNRT_LIBS}
    #${ADOLC_LIBRARY}
    )
ELSE(USE_ADOLC)
  SET(TARGET_LL
    cfsgeneral
    utils
    datainout
    paramh
    utils
    domain
    elemmapping
    ${CFS_FORTRAN_DYNRT_LIBS}
    )
ENDIF(USE_ADOLC)



IF(USE_XERCES)
  LIST(APPEND TARGET_LL ${XERCES_LIBRARY})
ENDIF(USE_XERCES)

IF(USE_LIBXML2)
  LIST(APPEND TARGET_LL ${LIBXML2_LIBRARY})
ENDIF(USE_LIBXML2)  

IF(USE_LIBFBI)
  LIST(APPEND TARGET_LL
    ${BOOST_THREAD_LIB}
  )
ENDIF(USE_LIBFBI)

IF(MINGW)
  LIST(APPEND TARGET_LL
    -lws2_32
    -lwsock32
    )
ELSEIF(DIST STREQUAL "FEDORA" AND REV STREQUAL "22")
  LIST(APPEND TARGET_LL
    -lpthread
    -llapack
    -lblas
    -lquadmath
    )
ELSEIF(DIST STREQUAL "CENTOS" AND REV STREQUAL "6.7")
  LIST(APPEND TARGET_LL
    -lpthread
########## RHEL is detected as CENTOS on Lima of RRZE: blas and lapack were not found, but we have mkl instead #############
#    -llapack            
#    -lblas
    )
ELSE()
  LIST(APPEND TARGET_LL
    -lpthread
    -lquadmath
     )
ENDIF()

# Target link libraries for CFS++ executable target.
TARGET_LINK_LIBRARIES(cfsbin 
  ${TARGET_LL}
  )


#IF(MINGW)
#  IF(CFS_BLAS_LAPACK STREQUAL MKL)
#    SET_TARGET_PROPERTIES(cfsbin PROPERTIES
#      LINK_FLAGS -Wl,--wrap,__chkstk
#      LINKER_LANGUAGE CXX)
#  ENDIF()
#ENDIF()

SET(CFS_LINK_FLAGS "")

# Make sure CFS++ gets linked with system libgcc.so instead of ANSYS libgcc.so.
IF(UNIX AND USE_ANSYSRST)
  SET(CFS_LINK_FLAGS "${CFS_LINK_FLAGS} -Wl,-rpath,/${LIB_SUFFIX}")
ENDIF(UNIX AND USE_ANSYSRST)

# Even if we did not compile CFS++ itself with OpenMP enabled, we may have
# compiled external libraries like LIS with OpenMP support.
IF(OPENMP_FOUND)
  SET(CFS_LINK_FLAGS "${CFS_LINK_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

SET_TARGET_PROPERTIES(cfsbin PROPERTIES LINK_FLAGS "${CFS_LINK_FLAGS}")

# Install the CFS++ executable target to ${CMAKE_INSTALL_PREFIX}/bin
#INSTALL_TARGETS(/bin cfs)
