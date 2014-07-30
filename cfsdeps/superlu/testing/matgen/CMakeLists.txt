#######################################################################
#  CMakeLists.txt file used to create a library of the test matrix
#  generators used in LAPACK.  The files are organized as follows:
#
#     SCATGEN  -- Auxiliary routines called from both REAL and COMPLEX
#     DZATGEN  -- Auxiliary routines called from both DOUBLE PRECISION
#                 and COMPLEX*16
#     SMATGEN  -- Single precision real matrix generation routines
#     CMATGEN  -- Single precision complex matrix generation routines
#     DMATGEN  -- Double precision real matrix generation routines
#     ZMATGEN  -- Double precision complex matrix generation routines

#######################################################################

set(ALLAUX lsamen.c)
set(SCATGEN slatm1.c slaran.c slarnd.c slaruv.c slabad.c slarnv.c)
set(SLASRC  slatb4.c slaset.c slartg.c)
set(SMATGEN slatms.c slatme.c slatmr.c slagge.c slagsy.c 
  slarge.c slaror.c slarot.c slatm2.c slatm3.c)
set(SINTRINSIC  r_lg10.c r_sign.c pow_dd.c)

set(DZATGEN dlatm1.c dlaran.c dlarnd.c dlaruv.c dlabad.c dlarnv.c)
set(DLASRC  dlatb4.c dlaset.c dlartg.c)
set(DMATGEN dlatms.c dlatme.c dlatmr.c dlagge.c dlagsy.c dlarge.c
  dlaror.c dlarot.c dlatm2.c dlatm3.c)
set(DINTRINSIC  d_lg10.c d_sign.c pow_dd.c)

set(CLASRC  clatb4.c claset.c clartg.c clarnv.c clacgv.c csymv.c)
set(CMATGEN clatms.c clatme.c clatmr.c clagge.c clagsy.c clarge.c
  claror.c clarot.c clatm2.c clatm3.c claghe.c clarnd.c
  cdotc.c)

set(ZLASRC  zlatb4.c zlaset.c zlartg.c zlarnv.c zlacgv.c zsymv.c)
set(ZMATGEN zlatms.c zlatme.c zlatmr.c zlagge.c zlagsy.c 
  zlarge.c zlaror.c zlarot.c zlatm2.c zlatm3.c zlaghe.c 
  zlarnd.c zdotc.c)

set(matgen_SRCS ${ALLAUX} ${SMATGEN} ${SCATGEN} ${SLASRC} ${SINTRINSIC} 
  ${DMATGEN} ${DZATGEN} ${DLASRC} ${DINTRINSIC} ${CMATGEN} ${CLASRC} 
  ${ZMATGEN} ${ZLASRC}
  )

add_library(tmglib STATIC ${matgen_SRCS})
target_link_libraries(tmglib ${BLAS_LIBRARIES})
set_target_properties(tmglib PROPERTIES COMPILE_FLAGS "-DPRNTlevel=0")

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
