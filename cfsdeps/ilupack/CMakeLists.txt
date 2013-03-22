CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

PROJECT(ilupack C Fortran)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(INSTALL_BIN_DIR
  "${CMAKE_INSTALL_PREFIX}/bin/${CFS_ARCH_STR}"
  CACHE PATH "Installation directory for executables")
set(INSTALL_LIB_DIR
  "${CMAKE_INSTALL_PREFIX}/${LIB_SUFFIX}/${CFS_ARCH_STR}"
  CACHE PATH "Installation directory for libraries")
set(INSTALL_INC_DIR
  "${CMAKE_INSTALL_PREFIX}/include"
  CACHE PATH "Installation directory for headers")


IF(APPLE)
#  ADD_DEFINITIONS("-D_DECLARE_ILUPACK_GLOBALS_")
  ADD_DEFINITIONS("-fno-common")
ENDIF(APPLE)

# Decide which type of matching should be used
SET(MATCHING_TYPE "MUMPS"
  CACHE STRING
  "Type matching to use. May be one of MUMPS, PARDISO or MC64")

# Should we build samples also?
OPTION(BUILD_SAMPLES "OFF")

# Let the user choose the type of floating point numbers for ilupack
SET(FLOAT_TYPE "DOUBLE_REAL"
  CACHE STRING
  "Type of floating points numbers to use. May be one of SINGLE_REAL, DOUBLE_REAL, SINGLE_COMPLEX, DOUBLE_COMPLEX")

#-----------------------------------------------------------------------------
# Include all the necessary files for macros
#-----------------------------------------------------------------------------
INCLUDE (${CMAKE_ROOT}/Modules/CheckIncludeFile.cmake)

CHECK_INCLUDE_FILE ("sys/times.h"   HAVE_SYS_TIMES_H)
CHECK_INCLUDE_FILE ("sys/time.h"    HAVE_SYS_TIME_H)

CONFIGURE_FILE( 
  "include/ilupack_config.h.in"
  "${CMAKE_SOURCE_DIR}/include/ilupack_config.h"
)

# Set internal MATCHING variable according to MATCHING_TYPE
IF(MATCHING_TYPE STREQUAL "MUMPS")
  SET(MATCHING "_MUMPS_MATCHING_")
ENDIF(MATCHING_TYPE STREQUAL "MUMPS")

IF(MATCHING_TYPE STREQUAL "PARDISO")
  SET(MATCHING "_PARDISO_MATCHING_")
ENDIF(MATCHING_TYPE STREQUAL "PARDISO")

IF(MATCHING_TYPE STREQUAL "MC64")
  SET(MATCHING "_MC64_MATCHING_")
ENDIF(MATCHING_TYPE STREQUAL "MC64")


# Set internal ARITHMETIC and FLOAT variables according to FLOAT_TYPE
IF(FLOAT_TYPE STREQUAL "SINGLE_COMPLEX")
  SET(ARITHMETIC "_SINGLE_COMPLEX_")
  SET(FLOAT "C")
ENDIF(FLOAT_TYPE STREQUAL "SINGLE_COMPLEX")

IF(FLOAT_TYPE STREQUAL "DOUBLE_COMPLEX")
  SET(ARITHMETIC "_DEFAULT_DOUBLE_COMPLEX_")
  SET(FLOAT "Z")
ENDIF(FLOAT_TYPE STREQUAL "DOUBLE_COMPLEX")

IF(FLOAT_TYPE STREQUAL "SINGLE_REAL")
  SET(ARITHMETIC "_SINGLE_REAL_")
  SET(FLOAT "S")
ENDIF(FLOAT_TYPE STREQUAL "SINGLE_REAL")

IF(FLOAT_TYPE STREQUAL "DOUBLE_REAL")
  SET(ARITHMETIC "_DOUBLE_REAL_")
  SET(FLOAT "D")
ENDIF(FLOAT_TYPE STREQUAL "DOUBLE_REAL")

# We want to build ILUPACK with the LP64 model.
# CFS is compiled with 32 integers (check via cfs --version) on all plattforms
# otherwise check CFS_ARCH and set for "X86_64" "_LONG_INTEGER_" and for gfortran
# SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-integer-8")
SET(LONGINTEGER "_NOT_USING_LONG_INTEGER_")

# switch for names alias of fortran routines
# -D__UNDERSCORE__   use this, if a fortran routine say "daxpy" has to be
#                    called from C using the name "daxpy_"
# -D__CAPS__         use this, if a fortran routine say "daxpy" has to be
#                    called from C using the name "DAXPY"
# -D__2UNDERSCORES__ use this, if a fortran routine say "daxpy" has to be
#                    called from C using the name "daxpy__"
# You can combine __CAPS__ with either __UNDERSCORE__ or __2UNDERSCORES__
# to obtain names like "DAXPY_" or "DAXPY__"
# SOLARIS OS:      -D__UNDERSCORE__
# Red Hat Linux:   -D__UNDERSCORE__
# AIX:             none
SET(FORTRANNAMES "__UNDERSCORE__")

IF(WIN32)
  IF(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    SET(FORTRANNAMES "__CAPS__")
  ENDIF()

  if(MSVC)
    set(CMAKE_DEBUG_POSTFIX "d")
  endif()
ELSE(WIN32)

  # If compiler is gfortran add a flag that implicit variable declarations
  # produce an error
  IF(CMAKE_COMPILER_IS_GNUFortran)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fimplicit-none")
  ENDIF()

ENDIF(WIN32)

SET(BLASLIKE_SRCS 
  #blaslike/cdotc2.F
  #blaslike/cdotu2.F
  #blaslike/cprivatehptrs.F
  #blaslike/crot.F
  blaslike/ddot2.F
  blaslike/ddot2t.F
  blaslike/dprivatesptrs.F
  #blaslike/sdot2.F
  #blaslike/sprivatesptrs.F
  blaslike/zdotc2.F
  blaslike/zdotu2.F
  blaslike/zprivatehptrs.F)
  #blaslike/zrot.F)

SET_SOURCE_FILES_PROPERTIES(${BLASLIKE_SRCS}
  PROPERTIES COMPILE_FLAGS "-D${ARITHMETIC} -D${LONGINTEGER}")

INCLUDE_DIRECTORIES(include)

SET(BLASLIKE_LIBNAME blaslike)

ADD_LIBRARY(${BLASLIKE_LIBNAME} STATIC ${BLASLIKE_SRCS})

INSTALL(TARGETS ${BLASLIKE_LIBNAME}
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}")

# SPARSPAK
SET(SPARSPAK_SRCS 
  sparspak/rcm.F
  sparspak/rootls.F
  sparspak/degree.F
  sparspak/fndsep.F
  sparspak/fnroot.F
  sparspak/gennd.F
  sparspak/genqmd.F
  sparspak/genrcm.F
  sparspak/qmdmrg.F
  sparspak/qmdqt.F
  sparspak/qmdrch.F
  sparspak/qmdupd.F
  sparspak/revrse.c)

SET_SOURCE_FILES_PROPERTIES(${SPARPAK_SRCS}
  PROPERTIES COMPILE_FLAGS "-D${MATCHING} -D${FORTRANNAMES} -D${ARITHMETIC} -D${LONGINTEGER}")

SET(SPARSPAK_LIBNAME sparspak)

ADD_LIBRARY(${SPARSPAK_LIBNAME} STATIC ${SPARSPAK_SRCS})

INSTALL(TARGETS ${SPARSPAK_LIBNAME}
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}")


SET(FORTRAN_SRCS 
  fortran/ilupackdelete.c
  fortran/ilupackfactor.c
  fortran/ilupackinfo.c
  fortran/ilupackinit.c
  fortran/ilupacknnz.c
  fortran/ilupacksol.c
  fortran/ilupacksolver.c
  fortran/spdilupackdelete.c
  fortran/spdilupackfactor.c
  fortran/spdilupackinfo.c
  fortran/spdilupackinit.c
  fortran/spdilupacknnz.c
  fortran/spdilupacksol.c
  fortran/spdilupacksolver.c
  fortran/symilupackdelete.c
  fortran/symilupackdeletes.c
  fortran/symilupackfactor.c
  fortran/symilupackfactors.c
  fortran/symilupackinfo.c
  fortran/symilupackinfos.c
  fortran/symilupackinit.c
  fortran/symilupackinits.c
  fortran/symilupacknnz.c
  fortran/symilupacknnzs.c
  fortran/symilupacksol.c
  fortran/symilupacksols.c
  fortran/symilupacksolver.c
  fortran/symilupacksolvers.c
  fortran/symspdilupackconvert.c
  fortran/unscale.c
)

SET_SOURCE_FILES_PROPERTIES(${FORTRAN_SRCS}
  PROPERTIES COMPILE_FLAGS "-D${MATCHING} -D${FORTRANNAMES} -D${ARITHMETIC} -D${LONGINTEGER}")

SET(NOTDISTRIBUTED_SRCS
  notdistributed/AMF4.F
)

SET_SOURCE_FILES_PROPERTIES(${FORTRAN_SRCS}
  PROPERTIES COMPILE_FLAGS "-D${ARITHMETIC} -D${LONGINTEGER}")

SET(SOLVERS_SRCS
  solvers/AMGdelete.c
  solvers/AMGgetparams.c
  solvers/AMGinit.c
  solvers/AMGsetparams.c
  solvers/AMGsetupparameters.c
  solvers/AMGsolver.c
  solvers/AMGsolver_real_prec.c
  solvers/AMGsolver_real_spd_prec.c
  solvers/AMGsolver_real_sym_prec.c
  solvers/AMGsolver_spd_prec.c
  solvers/AMGsolver_sym_prec.c
  solvers/AMGsolver_sym_precs.c
  solvers/spdAMGdelete.c
  solvers/spdAMGgetparams.c
  solvers/spdAMGinit.c
  solvers/spdAMGsetparams.c
  solvers/spdAMGsetupparameters.c
  solvers/spdAMGsolver.c
  solvers/spdAMGsolver_real_prec.c
  solvers/symAMGdelete.c
  solvers/symAMGdeletes.c
  solvers/symAMGgetparams.c
  solvers/symAMGgetparamss.c
  solvers/symAMGinit.c
  solvers/symAMGinits.c
  solvers/symAMGsetparams.c
  solvers/symAMGsetparamss.c
  solvers/symAMGsetupparameters.c
  solvers/symAMGsetupparameterss.c
  solvers/symAMGsolver.c
  solvers/symAMGsolver_real_prec.c
  solvers/symAMGsolver_real_precs.c
  solvers/symAMGsolver_real_spd_prec.c
  solvers/symAMGsolver_real_spd_precs.c
  solvers/symAMGsolvers.c
  solvers/symAMGsolver_spd_prec.c
  solvers/symAMGsolver_spd_precs.c
  solvers/symspd.c)

IF(MATCHING STREQUAL "_PARDISO_MATCHING_")
  SET(SOLVERS_SRCS
    ${SOLVERS_SRCS}
    solvers/symAMGinit_pardiso.c)
ENDIF(MATCHING STREQUAL "_PARDISO_MATCHING_")

SET_SOURCE_FILES_PROPERTIES(${SOLVERS_SRCS}
  PROPERTIES COMPILE_FLAGS "-D${MATCHING} -D${FORTRANNAMES} -D${ARITHMETIC} -D${LONGINTEGER}")

SET(TOOLS_SRCS
  tools/cc_etimes.c
  tools/clear.c
  tools/cpermC.c
  tools/CSRmathvec.c
  tools/CSRmattvec.c
  tools/CSRmatvec.c
  tools/CSRSetupGraph.c
  tools/CSRSetupGraph_epsilon.c
  tools/CSRSparTran.c
  tools/csymmatvec.c
  tools/etree.c
  tools/etree_postorder.c
  tools/geteps.c
  tools/globals.c
  tools/GNLSYM.c
  tools/GNLSYMs.c
  tools/Malloc.c
  tools/mps_arms.c
  tools/qqsort2.F
  tools/qqsort.F
  tools/qqsorti.F
  tools/qqsorts2.F
  tools/qqsorts.F
  tools/qsort2.F
  tools/qsort.F
  tools/Realloc.c
  tools/rpermR.c
  tools/skewmatvec.c
  tools/swapj.c
  tools/swapm.c
  tools/symmatvec.c)

SET_SOURCE_FILES_PROPERTIES(${TOOLS_SRCS}
  PROPERTIES COMPILE_FLAGS "-D${MATCHING} -D${FORTRANNAMES} -D${ARITHMETIC} -D${LONGINTEGER}")

SET_SOURCE_FILES_PROPERTIES(tools/globals.c 
    PROPERTIES COMPILE_FLAGS "-DILUPACK_LIB_FLOATNESS=${FLOAT} -D${MATCHING} -D${FORTRANNAMES} -D${ARITHMETIC} -D${LONGINTEGER}")

INCLUDE_DIRECTORIES(sparskit)

SET(SPARSKIT_SRCS
  sparskit/gmres.F
  sparskit/pcg.F
  sparskit/fpcg.F
  sparskit/readmtc.F
  sparskit/writemtc.F
  sparskit/lusol.F
  sparskit/lutsol.F
  sparskit/lulsol.F
  sparskit/lutlsol.F
  sparskit/luusol.F
  sparskit/lutusol.F
  sparskit/ludlsol.F
  sparskit/lutdlsol.F
  sparskit/ludusol.F
  sparskit/lutdusol.F
  sparskit/ilutp2.F
  sparskit/bisinit.F
  sparskit/brkdn.F
  sparskit/mgsro.F
  sparskit/stopbis.F
  sparskit/tidycg.F
  sparskit/qsplit.F
  sparskit/qsplit2.F
  sparskit/qsplit3.F
  sparskit/scale.F
  sparskit/spdscale.F
  sparskit/symscale.F
  sparskit/csrcsc.F
  sparskit/fgmres.F
  sparskit/ilut2.F
  sparskit/bcg.F
  sparskit/readvectors.F
  sparskit/writevectors.F
  sparskit/sbcg.F
  sparskit/sbcgs.F
  sparskit/sqmr.F
  sparskit/sqmrs.F)



SET_SOURCE_FILES_PROPERTIES(${SPARSKIT_SRCS}
  PROPERTIES COMPILE_FLAGS "-D${ARITHMETIC} -D${LONGINTEGER}")

SET(ORDERINGS_SRCS
  orderings/qsortR2I.c
  orderings/indset.c
  orderings/permindset.c
  orderings/permrcm.c
  orderings/perm_rcm_FC.c
  orderings/perm_rcm_FCv.c
  orderings/symperm_rcm_FC.c
  orderings/symperm_rcm_FCv.c
  orderings/symperm_rcm_FCs.c
  orderings/symperm_rcm_FCvs.c
  orderings/permnd.c
  orderings/permnd_FC.c
  orderings/permnd_FCv.c
  orderings/symperm_nd_FC.c
  orderings/symperm_nd_FCv.c
  orderings/symperm_nd_FCs.c
  orderings/symperm_nd_FCvs.c
  orderings/permamf.c
  orderings/perm_amf_FC.c
  orderings/perm_amf_FCv.c
  orderings/symperm_amf_FC.c
  orderings/symperm_amf_FCv.c
  orderings/symperm_amf_FCs.c
  orderings/symperm_amf_FCvs.c
  orderings/permmmd.c
  orderings/perm_mmd_FC.c
  orderings/perm_mmd_FCv.c
  orderings/symperm_mmd_FC.c
  orderings/symperm_mmd_FCv.c
  orderings/symperm_mmd_FCs.c
  orderings/symperm_mmd_FCvs.c
  orderings/permnull.c
  orderings/permPQ.c
  orderings/spdpermrcm.c
  orderings/spdpermPP.c
  orderings/indPPF3.c
  orderings/indPQF5.c
  orderings/indFC.c
  orderings/permFC.c
  orderings/symindFC.c
  orderings/sympindFC.c
  orderings/sympindFCs.c
  orderings/sympindFCv.c
  orderings/sympindFCvs.c
  orderings/symbuildblock.c
  orderings/symbuildblocks.c
  orderings/sympermFC.c
  orderings/indFCv.c
  orderings/permFCv.c
  orderings/symindFCv.c
  orderings/symindFCvs.c
  orderings/sympermFCv.c
  orderings/perm_metis_E.c
  orderings/perm_metis_E_FC.c
  orderings/perm_metis_E_FCv.c
  orderings/symperm_metis_e_FC.c
  orderings/symperm_metis_e_FCv.c
  orderings/symperm_metis_e_FCs.c
  orderings/symperm_metis_e_FCvs.c
  orderings/perm_metis_N.c
  orderings/perm_metis_N_FC.c
  orderings/perm_metis_N_FCv.c
  orderings/symperm_metis_n_FC.c
  orderings/symperm_metis_n_FCv.c
  orderings/symperm_metis_n_FCs.c
  orderings/symperm_metis_n_FCvs.c
  orderings/permamd.c
  orderings/perm_amd_FC.c
  orderings/perm_amd_FCv.c
  orderings/symperm_amd_FC.c
  orderings/symperm_amd_FCv.c
  orderings/symperm_amd_FCs.c
  orderings/symperm_amd_FCvs.c
  orderings/permmwmrcm.c
  orderings/permmwmmmd.c
  orderings/permmwmamf.c
  orderings/permmwm_metis_e.c
  orderings/permmwm_metis_n.c
  orderings/permmwmamd.c
  orderings/sympermmwm_metis_n.c
  orderings/sympermmwm_metis_ns.c
  orderings/symperm_mwm_metis_n_FC.c
  orderings/symperm_mwm_metis_n_FCv.c
  orderings/symperm_mwm_metis_n_FCs.c
  orderings/symperm_mwm_metis_n_FCvs.c
  orderings/sympermmwm_metis_e.c
  orderings/sympermmwm_metis_es.c
  orderings/symperm_mwm_metis_e_FC.c
  orderings/symperm_mwm_metis_e_FCv.c
  orderings/symperm_mwm_metis_e_FCs.c
  orderings/symperm_mwm_metis_e_FCvs.c
  orderings/sympermmwm_amd.c
  orderings/sympermmwm_amds.c
  orderings/symperm_mwm_amd_FC.c
  orderings/symperm_mwm_amd_FCv.c
  orderings/symperm_mwm_amd_FCs.c
  orderings/symperm_mwm_amd_FCvs.c
  orderings/sympermmwm_mmd.c
  orderings/sympermmwm_mmds.c
  orderings/symperm_mwm_mmd_FC.c
  orderings/symperm_mwm_mmd_FCv.c
  orderings/symperm_mwm_mmd_FCs.c
  orderings/symperm_mwm_mmd_FCvs.c
  orderings/sympermmwm_amf.c
  orderings/sympermmwm_amfs.c
  orderings/symperm_mwm_amf_FC.c
  orderings/symperm_mwm_amf_FCv.c
  orderings/symperm_mwm_amf_FCs.c
  orderings/symperm_mwm_amf_FCvs.c
  orderings/sympermmwm_rcm.c
  orderings/sympermmwm_rcms.c
  orderings/symperm_mwm_rcm_FC.c
  orderings/symperm_mwm_rcm_FCv.c
  orderings/symperm_mwm_rcm_FCs.c
  orderings/symperm_mwm_rcm_FCvs.c
  orderings/permMC64rcm.c
  orderings/permMC64mmd.c
  orderings/permMC64amf.c
  orderings/permMC64_metis_e.c
  orderings/permMC64_metis_n.c
  orderings/permMC64amd.c
  orderings/permMWA_rcm.c
  orderings/permMWA_mmd.c
  orderings/permMWA_amf.c
  orderings/permMWA_amd.c
  orderings/permMWA_metis_e.c
  orderings/permMWA_metis_n.c
  orderings/sympermMC64rcm.c
  orderings/sympermMC64rcms.c
  orderings/symperm_MC64_rcm_FC.c
  orderings/symperm_MC64_rcm_FCv.c
  orderings/symperm_MC64_rcm_FCs.c
  orderings/symperm_MC64_rcm_FCvs.c
  orderings/sympermMC64_metis_n.c
  orderings/sympermMC64_metis_ns.c
  orderings/symperm_MC64_metis_n_sp.c
  orderings/symperm_MC64_metis_n_sps.c
  orderings/symperm_MC64_metis_n_FC.c
  orderings/symperm_MC64_metis_n_FCv.c
  orderings/symperm_MC64_metis_n_FCs.c
  orderings/symperm_MC64_metis_n_FCvs.c
  orderings/sympermMC64_metis_e.c
  orderings/sympermMC64_metis_es.c
  orderings/symperm_MC64_metis_e_FC.c
  orderings/symperm_MC64_metis_e_FCv.c
  orderings/symperm_MC64_metis_e_FCs.c
  orderings/symperm_MC64_metis_e_FCvs.c
  orderings/sympermMC64amd.c
  orderings/sympermMC64amds.c
  orderings/symperm_MC64_amd_FC.c
  orderings/symperm_MC64_amd_FCv.c
  orderings/symperm_MC64_amd_FCs.c
  orderings/symperm_MC64_amd_FCvs.c
  orderings/sympermMC64mmd.c
  orderings/sympermMC64mmds.c
  orderings/symperm_MC64_mmd_FC.c
  orderings/symperm_MC64_mmd_FCv.c
  orderings/symperm_MC64_mmd_FCs.c
  orderings/symperm_MC64_mmd_FCvs.c
  orderings/sympermMC64amf.c
  orderings/sympermMC64amfs.c
  orderings/symperm_MC64_amf_FC.c
  orderings/symperm_MC64_amf_FCv.c
  orderings/symperm_MC64_amf_FCs.c
  orderings/symperm_MC64_amf_FCvs.c
  orderings/symwms.c
  orderings/symwm.c
  orderings/swm.c
  orderings/perm_matching_metis_n.c
  orderings/symperm_matching_metis_n.c
  orderings/symperm_matching_metis_ns.c
  orderings/symperm_matching_metis_n_FC.c
  orderings/symperm_matching_metis_n_FCs.c
  orderings/symperm_matching_metis_n_FCv.c
  orderings/symperm_matching_metis_n_FCvs.c
#  The following files have been excluded from build due to problems with
#  saddle point code.
#  orderings/symperm_matching_metis_n_sp.c
#  orderings/symperm_matching_metis_n_sps.c
  orderings/perm_matching_metis_e.c
  orderings/symperm_matching_metis_e.c
  orderings/symperm_matching_metis_es.c
  orderings/symperm_matching_metis_e_FC.c
  orderings/symperm_matching_metis_e_FCs.c
  orderings/symperm_matching_metis_e_FCv.c
  orderings/symperm_matching_metis_e_FCvs.c
  orderings/perm_matching_amf.c
  orderings/symperm_matching_amf.c
  orderings/symperm_matching_amd.c
  orderings/symperm_matching_amfs.c
  orderings/symperm_matching_amf_FC.c
  orderings/symperm_matching_amf_FCs.c
  orderings/symperm_matching_amf_FCv.c
  orderings/symperm_matching_amf_FCvs.c
  orderings/perm_matching_amd.c
  orderings/symperm_matching_amds.c
  orderings/symperm_matching_amd_FC.c
  orderings/symperm_matching_amd_FCs.c
  orderings/symperm_matching_amd_FCv.c
  orderings/symperm_matching_amd_FCvs.c
  orderings/perm_matching_mmd.c
  orderings/symperm_matching_mmd.c
  orderings/symperm_matching_mmds.c
  orderings/symperm_matching_mmd_FC.c
  orderings/symperm_matching_mmd_FCs.c
  orderings/symperm_matching_mmd_FCv.c
  orderings/symperm_matching_mmd_FCvs.c
  orderings/perm_matching_rcm.c
  orderings/symperm_matching_rcm.c
  orderings/symperm_matching_rcms.c
  orderings/symperm_matching_rcm_FC.c
  orderings/symperm_matching_rcm_FCs.c
  orderings/symperm_matching_rcm_FCv.c
  orderings/symperm_matching_rcm_FCvs.c

#  amd/amd_1.c
#  amd/amd_2.c
#  amd/amd_aat.c
#  amd/amd_control.c
#  amd/amd_defaults.c
#  amd/amd_dump.c
#  amd/amd_info.c
#  amd/amd_order.c
#  amd/amd_postorder.c
#  amd/amd_post_tree.c
  amd/ilupack_amd_preprocess.c
#  amd/amd_valid.c
)

# Sources for saddle point problems which do not work at the moment
#  orderings/symperm_matching_amd_sp.c
#  orderings/symperm_matching_amf_sp.c
#  orderings/symperm_matching_metis_e_sp.c
#  orderings/symperm_matching_metis_n_sp.c
#  orderings/symperm_matching_mmd_sp.c
#  orderings/symperm_matching_rcm_sp.c
#  orderings/symperm_MC64_amd_sp.c
#  orderings/symperm_MC64_amf_sp.c
#  orderings/symperm_MC64_metis_e_sp.c
#  orderings/symperm_MC64_metis_n_sp.c
#  orderings/symperm_MC64_mmd_sp.c
#  orderings/symperm_MC64_rcm_sp.c
#  orderings/symperm_mwm_amd_sp.c
#  orderings/symperm_mwm_amf_sp.c
#  orderings/symperm_mwm_metis_e_sp.c
#  orderings/symperm_mwm_metis_n_sp.c
#  orderings/symperm_mwm_mmd_sp.c
#  orderings/symperm_mwm_rcm_sp.c)


SET_SOURCE_FILES_PROPERTIES(${ORDERINGS_SRCS}
  PROPERTIES COMPILE_FLAGS "-D${MATCHING} -D${FORTRANNAMES} -D${ARITHMETIC} -D${LONGINTEGER}")

INCLUDE_DIRECTORIES(preconditioners)

SET(PRECONDITIONERS_SRCS
  preconditioners/AMGdlsol.c
  preconditioners/AMGdusol.c
  preconditioners/AMGextract.c
  preconditioners/AMGlsol.c
  preconditioners/AMGsetup.c
  preconditioners/AMGsol1.c
  preconditioners/AMGsol2.c
  preconditioners/AMGsol.c
  preconditioners/AMGtdlsol.c
  preconditioners/AMGtdusol.c
  preconditioners/AMGtlsol.c
  preconditioners/AMGtsol.c
  preconditioners/AMGtusol.c
  preconditioners/AMGusol.c
  preconditioners/csymAMGextract.c
  preconditioners/ildlc.F
  preconditioners/ildlcs.F
  preconditioners/ildlcsol.F
  preconditioners/ildlcsols.F
  preconditioners/ildlcsolz.F
  preconditioners/ildlcsolzs.F
  preconditioners/ildlcz.F
  preconditioners/ildlczs.F
  preconditioners/ilucdlsol.F
  preconditioners/ilucdusol.F
  preconditioners/iluc.F
  preconditioners/iluclsol.F
  preconditioners/ilucsol.F
  preconditioners/iluctdlsol.F
  preconditioners/iluctdusol.F
  preconditioners/iluctlsol.F
  preconditioners/iluctsol.F
  preconditioners/iluctusol.F
  preconditioners/ilucusol.F
  preconditioners/LDLP.F
  preconditioners/LDLPsol.F
  preconditioners/LUPQdlsol.F
  preconditioners/LUPQdusol.F
  preconditioners/LUPQ.F
  preconditioners/LUPQlsol.F
  preconditioners/LUPQsol.F
  preconditioners/LUPQtdlsol.F
  preconditioners/LUPQtdusol.F
  preconditioners/LUPQtlsol.F
  preconditioners/LUPQtsol.F
  preconditioners/LUPQtusol.F
  preconditioners/LUPQusol.F
  preconditioners/mpildlc.F
  preconditioners/mpiluc.F
  preconditioners/pildlcdlsol.F
  preconditioners/pildlcdusol.F
  preconditioners/pildlc.F
  preconditioners/pildlclsol.F
  preconditioners/pildlcusol.F
  preconditioners/pilucdlsol.F
  preconditioners/pilucdusol.F
  preconditioners/piluc.F
  preconditioners/piluclsol.F
  preconditioners/piluctdlsol.F
  preconditioners/piluctdusol.F
  preconditioners/piluctlsol.F
  preconditioners/piluctusol.F
  preconditioners/pilucusol.F
  preconditioners/psympiluclsol.F
  preconditioners/psympiluclsols.F
  preconditioners/psympilucusol.F
  preconditioners/psympilucusols.F
  preconditioners/shrAMGextract.c
  preconditioners/spdAMGsetup.c
  preconditioners/spdAMGsol1.c
  preconditioners/spdAMGsol2.c
  preconditioners/spdAMGsol.c
  preconditioners/spiluc.F
  preconditioners/ssmAMGextract.c
  preconditioners/symAMGextract.c
  preconditioners/symAMGsetup.c
  preconditioners/symAMGsetups.c
  preconditioners/symAMGsol1.c
  preconditioners/symAMGsol1s.c
  preconditioners/symAMGsol2.c
  preconditioners/symAMGsol2s.c
  preconditioners/symAMGsol.c
  preconditioners/symAMGsols.c
  preconditioners/symiluc.F
  preconditioners/symilucs.F
  preconditioners/sympiluc.F
  preconditioners/sympiluclsol.F
  preconditioners/sympiluclsols.F
  preconditioners/sympilucs.F
  preconditioners/sympilucsol.F
  preconditioners/sympilucsols.F
  preconditioners/sympilucusol.F
  preconditioners/sympilucusols.F)

SET_SOURCE_FILES_PROPERTIES(${PRECONDITIONERS_SRCS}
  PROPERTIES COMPILE_FLAGS "-D${MATCHING} -D${FORTRANNAMES} -D${ARITHMETIC} -D${LONGINTEGER}")


IF(MATCHING_TYPE STREQUAL "MUMPS")

  # Define a list of source files for MUMPS
  SET(MUMPS_SRCS
    mumps/Dmumps_struc_def.F
    mumps/Dmumps_match.F
    mumps/Dmumps_203.F
    mumps/Dmumps_444.F
    mumps/Dmumps_445.F
    mumps/Dmumps_446.F
    mumps/Dmumps_447.F
    mumps/Dmumps_448.F
    mumps/Dmumps_450.F
    mumps/Dmumps_451.F
    mumps/Dmumps_452.F
    mumps/Dmumps_453.F
    mumps/Dmumps_454.F
    mumps/Dmumps_455.F
    mumps/Dmumps_457.F
    mumps/Dmumps_551.F
    mumps/Dmumps_559.F
    mumps/Dmumps_562.F
    mumps/Dmumps_563.F
    mumps/Dmumps_abort.F
    mumps/Dmumps_metric2x2.F
    mumps/Dmumps_updatescore.F
    mumps/Dmumps_update_inverse.F)
  
  # Rename source files according to ${FLOAT}
  STRING(REPLACE "Dmumps" "${FLOAT}mumps" MUMPS_SRCS ${MUMPS_SRCS})
  
  # Make sure that MUMPS_SRCS will still get interpreted as list
  # of strings instead of a single string by adding semicolons
  # after each item
  STRING(REPLACE "F" "F;" MUMPS_SRCS ${MUMPS_SRCS})
  
  SET_SOURCE_FILES_PROPERTIES(${MUMPS_SRCS}
    PROPERTIES COMPILE_FLAGS "-D${ARITHMETIC} -D${LONGINTEGER}")
  
ENDIF(MATCHING_TYPE STREQUAL "MUMPS")


# Set name for ilupack library
SET(ILUPACK_LIBNAME "${FLOAT}ilupack")

IF(MSVC AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    SET(ILUPACK_LIBNAME ${ILUPACK_LIBNAME}d)
ENDIF(MSVC AND CMAKE_BUILD_TYPE STREQUAL "Debug")

# Build an overall list of ilupack sources
SET(ILUPACK_SRCS
  ${FORTRAN_SRCS}
  ${NOTDISTRIBUTED_SRCS}
  ${SOLVERS_SRCS}
  ${SPARSKIT_SRCS}
  ${ORDERINGS_SRCS}
  ${PRECONDITIONERS_SRCS}
  ${TOOLS_SRCS})

# If matching type is equal to MUMPS add sources to ILUPACK_SRCS
IF(MATCHING_TYPE STREQUAL "MUMPS")
  SET(ILUPACK_SRCS
    ${ILUPACK_SRCS}
    ${MUMPS_SRCS})
ENDIF(MATCHING_TYPE STREQUAL "MUMPS")

ADD_LIBRARY(${ILUPACK_LIBNAME} STATIC
  ${ILUPACK_SRCS}
)

# Put ilupack sources into a lib directory when installing
INSTALL(TARGETS ${ILUPACK_LIBNAME}
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}")

INSTALL(FILES
  "include/ilupack.h"
  "include/f2c.h"
  "include/long_integer.h"
  "include/namesilupack.h"
  DESTINATION "${INSTALL_INC_DIR}")

IF(BUILD_SAMPLES)
  # MESSAGE("Query ${PREFIX}/${LIB_SUFFIX}/${CFS_ARCH_STR}")

  # The libraries AMD, BLAS, LAPACK and Metis are built
  # as standalone libraries independently of ilupack.
  # Therefore they are searched for in a given directory
  # only if we need them to build the samples.
  FIND_LIBRARY(AMD_LIB amd
    "${PREFIX}/${LIB_SUFFIX}/${CFS_ARCH_STR}"
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH)
  
  FIND_LIBRARY(BLAS_LIB blas
    "${PREFIX}/${LIB_SUFFIX}/${CFS_ARCH_STR}"
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH)
  
  FIND_LIBRARY(LAPACK_LIB lapack
    "${PREFIX}/${LIB_SUFFIX}/${CFS_ARCH_STR}"
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH)
  
  FIND_LIBRARY(METIS_LIB metis
    "${PREFIX}/${LIB_SUFFIX}/${CFS_ARCH_STR}"
    NO_DEFAULT_PATH
    NO_CAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH)

  SET(MAIN_SRCS samples/main.c)
  SET(MAINSPD_SRCS samples/mainspd.c)
  SET(MAINSYM_SRCS samples/mainsym.c)
  SET(MAINSYMS_SRCS samples/mainsyms.c)
  SET(MAINGNL_SRCS samples/maingnl.c)
#  SET(FMAINSYMS_SRCS samples/fmainsyms.F)
#  SET(MAINILUC_SRCS samples/mainiluc.c)
  SET(MAINSYM3S_SRCS samples/mainsym3s.c)

  SET_SOURCE_FILES_PROPERTIES(
    ${MAIN_SRCS}
    ${MAINSPD_SRCS}
    ${MAINSYM_SRCS}
    ${MAINSYMS_SRCS}
    ${MAINGNL_SRCS}
 #   ${FMAINSYMS_SRCS}
 #   ${MAINILUC_SRCS}
    ${MAINSYM3S_SRCS}
     PROPERTIES COMPILE_FLAGS "-D${MATCHING} -D${FORTRANNAMES} -D${ARITHMETIC} -D${LONGINTEGER}")

  SET(TARGET_LL
    ${ILUPACK_LIBNAME})

  IF(FLOAT STREQUAL "Z")
    LIST(APPEND TARGET_LL
      ${INSTALL_LIB_DIR}/libDilupack.a
      )
  ENDIF(FLOAT STREQUAL "Z")

  IF(FLOAT STREQUAL "C")
    LIST(APPEND TARGET_LL
      ${INSTALL_LIB_DIR}/libSilupack.a
      )
  ENDIF(FLOAT STREQUAL "C")

  LIST(APPEND TARGET_LL
    ${BLASLIKE_LIBNAME}
    ${SPARSPAK_LIBNAME}
    ${AMD_LIB}
    ${LAPACK_LIB}
    ${BLAS_LIB}
    ${METIS_LIB}
  )

IF(CMAKE_COMPILER_IS_GNUC)
  LIST(APPEND TARGET_LL m)
ELSEIF(CMAKE_C_COMPILER_ID STREQUAL "Intel")
  IF(UNIX)
    LIST(APPEND TARGET_LL m)
  ENDIF()
ENDIF()


  # Make sure that C executables get linked against Fortran runtime libs.
  # SET(CMAKE_Fortran_LINKER_PREFERENCE_PROPAGATES TRUE)

  # Search explicitely for implicitely defined Fortran libs
  foreach(lib IN LISTS CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)
    FIND_LIBRARY(${lib}_LIBRARY
      NAMES ${lib}
      PATHS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}
      NO_DEFAULT_PATH
      NO_CMAKE_ENVIRONMENT_PATH
      NO_CMAKE_PATH
      NO_SYSTEM_ENVIRONMENT_PATH
      NO_CMAKE_SYSTEM_PATH 
      )

    LIST(APPEND TARGET_LL "${${lib}_LIBRARY}")
  endforeach()

  UNSET(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)

#  IF(CMAKE_Fortran_COMPILER MATCHES "ifort")
#    STRING(REGEX REPLACE "bin/ifort" "lib" IFORT_LIB_PATH "${CMAKE_Fortran_COMPILER}")
#    LINK_DIRECTORIES(${IFORT_LIB_PATH})
#    SET(TARGET_LL
#      ${TARGET_LL}
#      ifcoremt_pic
#      irc
#      imf
#      dl
#       pthread)
#  ENDIF(CMAKE_Fortran_COMPILER MATCHES "ifort")

  ADD_EXECUTABLE("${FLOAT}main" ${MAIN_SRCS})
  TARGET_LINK_LIBRARIES("${FLOAT}main" ${TARGET_LL})

  ADD_EXECUTABLE("${FLOAT}mainspd" ${MAINSPD_SRCS})
  TARGET_LINK_LIBRARIES("${FLOAT}mainspd" ${TARGET_LL})

  ADD_EXECUTABLE("${FLOAT}mainsym" ${MAINSYM_SRCS})
  TARGET_LINK_LIBRARIES("${FLOAT}mainsym" ${TARGET_LL})

  ADD_EXECUTABLE("${FLOAT}mainsyms" ${MAINSYMS_SRCS})
  TARGET_LINK_LIBRARIES("${FLOAT}mainsyms" ${TARGET_LL})

  ADD_EXECUTABLE("${FLOAT}maingnl" ${MAINGNL_SRCS})
  TARGET_LINK_LIBRARIES("${FLOAT}maingnl" ${TARGET_LL})

 # ADD_EXECUTABLE("${FLOAT}fmainsyms" ${FMAINSYMS_SRCS})
 # TARGET_LINK_LIBRARIES("${FLOAT}fmainsyms" ${TARGET_LL})

 # ADD_EXECUTABLE("${FLOAT}mainiluc" ${MAINILUC_SRCS})
 # TARGET_LINK_LIBRARIES("${FLOAT}mainiluc" ${TARGET_LL})

 ADD_EXECUTABLE("${FLOAT}mainsym3s" ${MAINSYM3S_SRCS})
 TARGET_LINK_LIBRARIES("${FLOAT}mainsym3s" ${TARGET_LL})


  INSTALL(TARGETS
    ${FLOAT}main
    ${FLOAT}mainspd
    ${FLOAT}mainsym
    ${FLOAT}mainsyms
    ${FLOAT}maingnl
#    ${FLOAT}fmainsyms
#    ${FLOAT}mainiluc
    ${FLOAT}mainsym3s
    RUNTIME DESTINATION "${INSTALL_BIN_DIR}")

ENDIF(BUILD_SAMPLES)