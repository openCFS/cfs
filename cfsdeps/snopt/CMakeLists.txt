cmake_minimum_required(VERSION 2.8.6)

#IF(NOT DEFINED CMAKE_BUILD_TYPE)
#  SET (CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
#ENDIF()

PROJECT(SNOPT Fortran)
SET(SNOPT_DIR snopt7)

SET(SNOPT_SRCS 
     ${SNOPT_DIR}/src/sqopt.f
     ${SNOPT_DIR}/src/snopta.f
     ${SNOPT_DIR}/src/snoptb.f
     ${SNOPT_DIR}/src/snoptc.f
     ${SNOPT_DIR}/src/snoptq.f
     ${SNOPT_DIR}/src/npopt.f
     ${SNOPT_DIR}/src/sq02lib.f
     ${SNOPT_DIR}/src/sn02lib.f
     ${SNOPT_DIR}/src/np02lib.f
     ${SNOPT_DIR}/src/sn05wrpa.f
     ${SNOPT_DIR}/src/sn05wrpb.f
     ${SNOPT_DIR}/src/sn05wrpc.f
     ${SNOPT_DIR}/src/sn05wrpn.f
     ${SNOPT_DIR}/src/sn10mach.f
     ${SNOPT_DIR}/src/sn12ampl.f
     ${SNOPT_DIR}/src/sn17util.f
     ${SNOPT_DIR}/src/sn20amat.f
     ${SNOPT_DIR}/src/sn25bfac.f
     ${SNOPT_DIR}/src/sn27lu.f
     ${SNOPT_DIR}/src/sn30spec.f
     ${SNOPT_DIR}/src/sn35mps.f
     ${SNOPT_DIR}/src/sn37wrap.f
     ${SNOPT_DIR}/src/sn40bfil.f
     ${SNOPT_DIR}/src/sn50lp.f
     ${SNOPT_DIR}/src/sn55qp.f
     ${SNOPT_DIR}/src/sn56qncg.f
     ${SNOPT_DIR}/src/sn57qopt.f
     ${SNOPT_DIR}/src/sn60srch.f
     ${SNOPT_DIR}/src/sn65rmod.f
     ${SNOPT_DIR}/src/sn70nobj.f
     ${SNOPT_DIR}/src/sn80ncon.f
     ${SNOPT_DIR}/src/sn85hess.f
     ${SNOPT_DIR}/src/sn87sopt.f
     ${SNOPT_DIR}/src/sn90lmqn.f
     ${SNOPT_DIR}/src/sn95fmqn.f
   )
SET(BLAS_SRCS ${SNOPT_DIR}/src/sn15blas.f)
SET(SNPRNT_SRCS ${SNOPT_DIR}/src/sn03prnt.f)

SET(SNOPTAD_SRCS
     ${SNOPT_DIR}/src/snopta.f
     ${SNOPT_DIR}/src/sq02lib.f
     ${SNOPT_DIR}/src/sn02lib.f
     ${SNOPT_DIR}/src/sn10mach.f
     ${SNOPT_DIR}/src/sn12ampl.f
     ${SNOPT_DIR}/src/sn17util.f
     ${SNOPT_DIR}/src/sn20amat.f
     ${SNOPT_DIR}/src/sn25bfac.f
     ${SNOPT_DIR}/src/sn27lu.f
     ${SNOPT_DIR}/src/sn30spec.f
     ${SNOPT_DIR}/src/sn37wrap.f
     ${SNOPT_DIR}/src/sn40bfil.f
     ${SNOPT_DIR}/src/sn50lp.f
     ${SNOPT_DIR}/src/sn55qp.f
     ${SNOPT_DIR}/src/sn56qncg.f
     ${SNOPT_DIR}/src/sn57qopt.f
     ${SNOPT_DIR}/src/sn60srch.f
     ${SNOPT_DIR}/src/sn65rmod.f
     ${SNOPT_DIR}/src/sn70nobj.f
     ${SNOPT_DIR}/src/sn80ncon.f
     ${SNOPT_DIR}/src/sn85hess.f
     ${SNOPT_DIR}/src/sn87sopt.f
     ${SNOPT_DIR}/src/sn90lmqn.f
     ${SNOPT_DIR}/src/sn95fmqn.f
   )

SET(SNADIOPT_GEN_SRCS
     ${SNOPT_DIR}/snadiopt/common_lib_src/ad.f
     ${SNOPT_DIR}/snadiopt/common_lib_src/mtx.f
     ${SNOPT_DIR}/snadiopt/lib_src/nlp.f
   )

SET(SNADIOPT_DENSE_SRCS
     ${SNOPT_DIR}/snadiopt/common_lib_src/dad.f 
     ${SNOPT_DIR}/snadiopt/lib_src/daduserfg.f
   )

SET(SNADIOPT_SPARSE_SRCS
     ${SNOPT_DIR}/snadiopt/common_lib_src/sad.f 
     ${SNOPT_DIR}/snadiopt/lib_src/saduserfg.f
   )

SET(SNADIOPT_SNADFUN_SRCS ${SNOPT_DIR}/snadiopt/common_lib_src/sn05adwrp.f)

ADD_LIBRARY(snopt STATIC ${SNOPT_SRCS})
ADD_LIBRARY(snoptad STATIC ${SNOPTAD_SRCS})
ADD_LIBRARY(snprint STATIC ${SNPRNT_SRCS})
ADD_LIBRARY(blas_snopt STATIC ${BLAS_SRCS})
ADD_LIBRARY(snadiopt STATIC ${SNADIOPT_GEN_SRCS} ${SNADIOPT_DENSE_SRCS} ${SNADIOPT_SPARSE_SRCS}
	${SNADIOPT_SNADFUN_SRCS})

INSTALL(TARGETS snopt snoptad snadiopt snprint ARCHIVE DESTINATION "${LIB_SUFFIX}")
