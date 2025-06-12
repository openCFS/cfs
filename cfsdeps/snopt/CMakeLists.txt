# this file generates a static and on request a shared snopt lib (the later needs blas on the system) 
cmake_minimum_required(VERSION 3.16)

include(CMakePrintHelpers)

project(snopt Fortran)

# we don't need the shared lib for openCFS but probably for own use
option(BUILD_SHARED_LIBS "build also a dynamic snopt lib" OFF)
if(BUILD_SHARED_LIBS)
  find_package(BLAS REQUIRED) # for the shared lib!
endif() 

set(SNOPT_SRCS 
     src/sqopt.f
     src/snopta.f
     src/snoptb.f
     src/snoptc.f
     src/snoptq.f
     src/npopt.f
     src/sq02lib.f
     src/sn02lib.f
     src/np02lib.f
     src/sn05wrpa.f
     src/sn05wrpb.f
     src/sn05wrpc.f
     src/sn05wrpn.f
     src/sn10mach.f
     src/sn12ampl.f
     src/sn17util.f
     src/sn20amat.f
     src/sn25bfac.f
     src/sn27lu.f
     src/sn30spec.f
     src/sn35mps.f
     src/sn37wrap.f
     src/sn40bfil.f
     src/sn50lp.f
     src/sn55qp.f
     src/sn56qncg.f
     src/sn57qopt.f
     src/sn60srch.f
     src/sn65rmod.f
     src/sn70nobj.f
     src/sn80ncon.f
     src/sn85hess.f
     src/sn87sopt.f
     src/sn90lmqn.f
     src/sn95fmqn.f
     src/sn03prnt.f) # print stuff

set(BLAS_SRCS src/sn15blas.f)

set(SNOPTAD_SRCS
     src/snopta.f
     src/sq02lib.f
     src/sn02lib.f
     src/sn10mach.f
     src/sn12ampl.f
     src/sn17util.f
     src/sn20amat.f
     src/sn25bfac.f
     src/sn27lu.f
     src/sn30spec.f
     src/sn37wrap.f
     src/sn40bfil.f
     src/sn50lp.f
     src/sn55qp.f
     src/sn56qncg.f
     src/sn57qopt.f
     src/sn60srch.f
     src/sn65rmod.f
     src/sn70nobj.f
     src/sn80ncon.f
     src/sn85hess.f
     src/sn87sopt.f
     src/sn90lmqn.f
     src/sn95fmqn.f)

set(SNADIOPT_GEN_SRCS
     snadiopt/common_lib_src/ad.f
     snadiopt/common_lib_src/mtx.f
     snadiopt/lib_src/nlp.f)

set(SNADIOPT_DENSE_SRCS
     snadiopt/common_lib_src/dad.f 
     snadiopt/lib_src/daduserfg.f)

set(SNADIOPT_SPARSE_SRCS
     snadiopt/common_lib_src/sad.f 
     snadiopt/lib_src/saduserfg.f)

set(SNADIOPT_SNADFUN_SRCS snadiopt/common_lib_src/sn05adwrp.f)

add_library(snopt STATIC ${SNOPT_SRCS})
# for some very strange reason, the precompiled package containing only lib/libsnopt.a is unpacked flat to CMAKE_BINARY_DIR/libsnopt.a
# but only when executed by execute_process "${CMAKE_COMMAND}" -E tar xzf
# with another file all works?! cmake version 3.22.2 on macOS
install(TARGETS snopt ARCHIVE DESTINATION "${LIB_SUFFIX}")

if(BUILD_SHARED_LIBS)
  add_library(snopt_shared SHARED ${SNOPT_SRCS}) # shared lib
  set_target_properties(snopt_shared PROPERTIES OUTPUT_NAME "snopt") # to get a snopt.so, snopt.dylib
  target_link_libraries(snopt_shared ${BLAS_LIBRARIES})
  install(TARGETS snopt_shared ARCHIVE DESTINATION "${LIB_SUFFIX}")
endif()
