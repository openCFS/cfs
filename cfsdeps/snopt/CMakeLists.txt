# this file generates a static and on request a shared snopt lib (the later needs blas on the system)
cmake_minimum_required(VERSION 3.16)

include(CMakePrintHelpers)

project(snopt Fortran)

# we don't need the shared lib for openCFS but probably for own use
option(BUILD_SHARED_LIBS "build also a dynamic snopt lib" OFF)
if(BUILD_SHARED_LIBS)
  find_package(BLAS REQUIRED) # for the shared lib!
endif()

# The file list follows the variable SNOPT of the original src/Submakefile.
# Deliberately not included:
# - snblas.f: the bundled BLAS. openCFS links its own BLAS (see below for the shared lib).
# - snopth.f: defines snOptCH/snKerCH a second time and would clash with snoptch.f.
# - sn27lu77.f, sn27lu90.f: alternatives to sn27lu.f.
# - snopt_cwrap.c, sqopt_cwrap.c, snopt_wrapper.f90, sqopt_wrapper.f90, snoptProblem.cpp,
#   sqoptProblem.cpp: the official C/C++ interface. openCFS calls the Fortran entry points
#   directly, so we build Fortran only.
set(SNOPT_SRCS
     src/sqopt.f
     src/snopta.f
     src/snoptb.f
     src/snoptc.f
     src/snoptch.f
     src/snoptq.f
     src/npopt.f
     src/snctrl.f
     src/sq02slib.f
     src/sn02slib.f
     src/np02slib.f
     src/sq02tlib.f
     src/sn02tlib.f
     src/np02tlib.f
     src/sq02lib.f
     src/sn02lib.f
     src/np02lib.f
     src/sn04wrap.f
     src/sn04stats.f
     src/sn05wrpa.f
     src/sn05wrpb.f
     src/sn05wrpc.f
     src/sn05wrpch.f
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
     src/sn82qn.f
     src/sn83bfgs.f
     src/sn87sopt.f
     src/snfilewrapper.f
     src/sn01dum.f # dummy snCheck, called by all interfaces
     src/sn03prnt.f) # print stuff

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
