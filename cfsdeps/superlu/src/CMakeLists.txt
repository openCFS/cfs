#######################################################################
#  CMakeLists.txt file for sparse supernodal LU, implemented in ANSI C
#
#  The files are organized as follows:
#
#       ALLAUX -- Auxiliary routines called from all precisions of SuperLU
#       LAAUX  -- LAPACK auxiliary routines called from all precisions
#       SLASRC -- LAPACK single precision real routines
#       DLASRC -- LAPACK double precision real routines
#       CLASRC -- LAPACK single precision complex routines
#       ZLASRC -- LAPACK double precision complex routines
#       SCLAUX -- LAPACK Auxiliary routines called from both real and complex
#       DZLAUX -- LAPACK Auxiliary routines called from both double precision
#                 and complex*16
#       SLUSRC -- Single precision real SuperLU routines
#       DLUSRC -- Double precision real SuperLU routines
#       CLUSRC -- Single precision complex SuperLU routines
#       ZLUSRC -- Double precision complex SuperLU routines
#
#######################################################################

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)

if(NOT SUPERLU_PRECISIONS)
  set(SUPERLU_PRECISIONS single double complex complex16)
endif(NOT SUPERLU_PRECISIONS)
string(TOUPPER "${SUPERLU_PRECISIONS}" SUPERLU_PRECISIONS)

set(LAAUX_SRCS lsame.c xerbla.c)

set(ALLAUX_LU_SRCS
  superlu_timer.c util.c memory.c get_perm_c.c mmd.c
  sp_coletree.c sp_preorder.c sp_ienv.c relax_snode.c
  heap_relax_snode.c colamd.c ilu_relax_snode.c 
  ilu_heap_relax_snode.c mark_relax.c mc64ad.c 
  qselect.c
  )

# Per original makefiles, override target level settings and
# don't optimize this file's compilation (The "last" optimization
# flag on the gcc command line is the one gcc uses, so append the
# default optimizmation setting to the end of the compile flags.)
set_source_files_properties(superlu_timer.c PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS_DISABLE_OPTIMIZATION}")

set(SINGLE_LA_SRCS slacon.c)

set(SINGLE_AUX_SRCS slamch.c)
# Per original makefiles, override target level settings and
# don't optimize this file's compilation (The "last" optimization
# flag on the gcc command line is the one gcc uses, so append the
# default optimizmation setting to the end of the compile flags.)
set_source_files_properties(slamch.c PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS_DISABLE_OPTIMIZATION}")

set(SINGLE_LU_SRCS 
  sgssv.c sgssvx.c ssp_blas2.c ssp_blas3.c sgscon.c  
  slangs.c sgsequ.c slaqgs.c spivotgrowth.c 
  sgsrfs.c sgstrf.c sgstrs.c scopy_to_ucol.c 
  ssnode_dfs.c ssnode_bmod.c spanel_dfs.c spanel_bmod.c 
  sreadhb.c sreadrb.c sreadtriple.c scolumn_dfs.c
  scolumn_bmod.c spivotL.c spruneL.c smemory.c
  sutil.c smyblas2.c sgsisx.c sgsitrf.c sldperm.c
  ilu_sdrop_row.c ilu_ssnode_dfs.c ilu_scolumn_dfs.c
  ilu_spanel_dfs.c ilu_scopy_to_ucol.c ilu_spivotL.c
  sdiagonal.c
  )

set(DOUBLE_LA_SRCS dlacon.c)

set(DOUBLE_AUX_SRCS dlamch.c)

# Per original makefiles, override target level settings and
# don't optimize this file's compilation (The "last" optimization
# flag on the gcc command line is the one gcc uses, so append the
# default optimizmation setting to the end of the compile flags.)
set_source_files_properties(dlamch.c PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS_DISABLE_OPTIMIZATION}")

set(DOUBLE_LU_SRCS 
  dgssv.c dgssvx.c dsp_blas2.c dsp_blas3.c dgscon.c 
  dlangs.c dgsequ.c dlaqgs.c dpivotgrowth.c  
  dgsrfs.c dgstrf.c dgstrs.c dcopy_to_ucol.c 
  dsnode_dfs.c dsnode_bmod.c dpanel_dfs.c dpanel_bmod.c 
  dreadhb.c dreadrb.c dreadtriple.c dcolumn_dfs.c
  dcolumn_bmod.c dpivotL.c dpruneL.c dmemory.c
  dutil.c dmyblas2.c dgsisx.c dgsitrf.c dldperm.c
  ilu_ddrop_row.c ilu_dsnode_dfs.c ilu_dcolumn_dfs.c
  ilu_dpanel_dfs.c ilu_dcopy_to_ucol.c ilu_dpivotL.c
  ddiagonal.c
  )


set(COMPLEX_LA_SRCS clacon.c scsum1.c icmax1.c)
set(COMPLEX_AUX_SRCS slamch.c)
set(COMPLEX_LU_SRCS
  scomplex.c cgssv.c cgssvx.c csp_blas2.c csp_blas3.c cgscon.c 
  clangs.c cgsequ.c claqgs.c cpivotgrowth.c cgsrfs.c cgstrf.c 
  cgstrs.c ccopy_to_ucol.c csnode_dfs.c csnode_bmod.c 
  cpanel_dfs.c cpanel_bmod.c creadhb.c creadrb.c creadtriple.c 
  ccolumn_dfs.c ccolumn_bmod.c cpivotL.c cpruneL.c cmemory.c
  cutil.c cmyblas2.c cgsisx.c cgsitrf.c cldperm.c
  ilu_cdrop_row.c ilu_csnode_dfs.c ilu_ccolumn_dfs.c
  ilu_cpanel_dfs.c ilu_ccopy_to_ucol.c ilu_cpivotL.c
  cdiagonal.c
  )

set(COMPLEX16_LA_SRCS zlacon.c dzsum1.c izmax1.c)
set(COMPLEX16_AUX_SRCS dlamch.c)
set(COMPLEX16_LU_SRCS 
  dcomplex.c zgssv.c zgssvx.c zsp_blas2.c zsp_blas3.c zgscon.c 
  zlangs.c zgsequ.c zlaqgs.c zpivotgrowth.c zgsrfs.c zgstrf.c 
  zgstrs.c zcopy_to_ucol.c zsnode_dfs.c zsnode_bmod.c 
  zpanel_dfs.c zpanel_bmod.c zreadhb.c zreadrb.c zreadtriple.c 
  zcolumn_dfs.c zcolumn_bmod.c zpivotL.c zpruneL.c zmemory.c
  zutil.c zmyblas2.c zgsisx.c zgsitrf.c zldperm.c
  ilu_zdrop_row.c ilu_zsnode_dfs.c ilu_zcolumn_dfs.c
  ilu_zpanel_dfs.c ilu_zcopy_to_ucol.c ilu_zpivotL.c
  zdiagonal.c
  )

#set(SUPERLU_SRCS ${ALLAUX_LU_SRCS} ${LAAUX_SRCS})
set(SUPERLU_SRCS ${ALLAUX_LU_SRCS})
foreach(precision ${SUPERLU_PRECISIONS})
  set(SUPERLU_SRCS ${SUPERLU_SRCS} ${${precision}_LA_SRCS} ${${precision}_AUX_SRCS} ${${precision}_LU_SRCS})
endforeach(precision ${SUPERLU_PRECISIONS})
list(REMOVE_DUPLICATES SUPERLU_SRCS)

set(SUPERLU_INCLUDE_FILES
  colamd.h
  html_mainpage.h
  slu_cdefs.h
  slu_Cnames.h
  slu_dcomplex.h
  slu_ddefs.h
  slu_scomplex.h
  slu_sdefs.h
  slu_util.h
  slu_zdefs.h
  superlu_enum_consts.h
  supermatrix.h
  )

add_library(superlu STATIC ${SUPERLU_SRCS})
set_target_properties(superlu PROPERTIES COMPILE_FLAGS "-DPRNTlevel=0")
install(TARGETS superlu
  RUNTIME DESTINATION ${INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR})

install(FILES ${SUPERLU_INCLUDE_FILES} DESTINATION ${INSTALL_INC_DIR}/superlu)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
