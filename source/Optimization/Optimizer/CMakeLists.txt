SET(OPTIMIZER_SRCS
  BaseOptimizer.cc
  ShapeOptimizer.cc
  OptimalityCondition.cc
  EvaluateOnly.cc
  GradientCheck.cc
)

SET(TARGET_LL
  cfsgeneral
  datainout
  design
  domain
  driver
  logging
  matvec
  optimization
  paramh
  pde
  utils
)

# we need it almost always ?!
SET(DEPENDS
  "${LAPACK_LIBRARY}"
  "${BLAS_LIBRARY}"
  "${CFS_FORTRAN_LIBS}"
)

IF(USE_IPOPT)
  # gfortranbegin gcc_s
   #SET(TARGET_LL ${TARGET_LL} -Wl,-Bstatic ipopt coinhsl )
   SET(COMBINED ${IPOPT_LIBRARY} ${IPOPT_HSL_LIBRARY} )
   ADD_DEPENDENT_LIBRARY(ipopt-lib "${COMBINED}" "${DEPENDS}")
   
   SET(TARGET_LL  ${TARGET_LL} ipopt-lib1)

   INCLUDE_DIRECTORIES(${IPOPT_INCLUDE_DIR})
   
   SET(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} IPOPT.cc)
ENDIF(USE_IPOPT)

IF(USE_KNITRO)

   #SET(TARGET_LL ${TARGET_LL} ${KNITRO_LIBRARY})
   SET(TARGET_LL ${TARGET_LL} knitro700)
   
   INCLUDE_DIRECTORIES(${KNITRO_INCLUDE_DIR})
   
   SET(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} KNITRO.cc)
ENDIF(USE_KNITRO)

IF(USE_SCPIP)
   IF(NOT USE_FEAS_SCP)
      # For FeasSCP we need to compile SCPIP and SCPIPBase but must not link
      # the lib as most functions have the same name
      SET(TARGET_LL ${TARGET_LL} ${SCPIP_LIBRARY} ${CFS_FORTRAN_LIBS})
   ENDIF(NOT USE_FEAS_SCP)
   # message("${CFS_FORTRAN_LIBS}")
   # SCPIPBase is also the open source project C++SCPIP.
   # With this definition the file lives in both worlds! 
   ADD_DEFINITIONS(-DUSE_4_CFS)
   
   SET(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} SCPIP.cc SCPIPBase.cc)
ENDIF(USE_SCPIP)

IF(USE_FEAS_SCP)
   SET(FEAS_SCP_LIBRARY /home/fwein/code/cfs/debug/lib64/OPENSUSE_11.4_X86_64/libscpip40i.a)
   SET(IPOPT23_LIBRARY /home/fwein/code/cfs/debug/lib64/OPENSUSE_11.4_X86_64/libipopt23.a) 

   SET(COMBINED ${FEAS_SCP_LIBRARY}  ${IPOPT23_LIBRARY})
   ADD_DEPENDENT_LIBRARY(feas_scp-lib "${COMBINED}" "${DEPENDS}")
   SET(TARGET_LL ${TARGET_LL} feas_scp-lib1)
   
   # SET(TARGET_LL ${TARGET_LL} libdmumps.a  libmpiseq.a  libmumps_common.a  libpord.a  libscp40.a ${CFS_FORTRAN_LIBS})
   
   SET(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} FeasSCP.cc)
ENDIF(USE_FEAS_SCP)


IF(USE_SNOPT)
   SET(TARGET_LL ${TARGET_LL} ${SNOPT_LIBRARY} ${CFS_FORTRAN_LIBS})
   
   SET(OPTIMIZER_SRCS ${OPTIMIZER_SRCS} SnOptInterface.f SnOpt.cc)
ENDIF(USE_SNOPT)

ADD_DEFINITIONS(-DUSE_4_CFS)

ADD_LIBRARY(optimizer STATIC ${OPTIMIZER_SRCS})

TARGET_LINK_LIBRARIES(optimizer optimization ${TARGET_LL})
