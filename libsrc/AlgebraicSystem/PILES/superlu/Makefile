#
CCLIB = cc
C_DIR=../../../..
#
src = dgssv.c dgssvx.c \
        dsp_blas2.c dsp_blas3.c dgscon.c dlacon.c \
        dlangs.c dgsequ.c dlaqgs.c dpivotgrowth.c  \
        dgsrfs.c dgstrf.c dgstrs.c dcopy_to_ucol.c \
        dsnode_dfs.c dsnode_bmod.c \
        dpanel_dfs.c dpanel_bmod.c dreadhb.c \
        dcolumn_dfs.c dcolumn_bmod.c dpivotL.c dpruneL.c \
        dmemory.c dutil.c \
	superlu_timer.c lsame.c util.c memory.c get_perm_c.c mmd.c \
	sp_coletree.c sp_preorder.c sp_ienv.c relax_snode.c \
        xerbla.c colamd.c \
	dlamch.c \
	idamax.c dasum.c daxpy.c dcopy.c ddot.c dnrm2.c \
        drot.c dscal.c \
	dgemv.c dsymv.c dtrsv.c dger.c dsyr2.c \
	dmyblas2.c

#
lib = superlu-piles

libpath = libsrc/AlgebraicSystem/PILES/superlu
CFLAGS3 = -DAdd_
#
#
# make an include file for library
#
include $(C_DIR)/main/Makefile.option
#
include $(LIBSRC_DIR)/makefilemach/make.mach.$(MACHINE)
#
ARFLAGS = r 
#
LIBB = $(LIB_DIR)/lib$(lib).so
#
.SUFFIXES: .c .o

.c.o:
	$(CCLIB) $(CFLAGS1) $(CFLAGS2) $(CFLAGS3) $<
#
libgoal: $(LIB_DIR) $(LIBB)
#
$(LIBB): $(src:.c=.o)
	$(CCLIB) $(FLAG_MKSHLIB) -o $@ $(?:.c=.o) 

#
########## Make clean for each directory
#
clean:
	-rm *.o
#
######### Make depend 
#
depend: # $(src)
	makedepend $(INCL)  -f$(C_DIR)/main/$(Dependencies) $(FlagMakeDepend) $(src) 

# DO NOT DELETE THIS LINE -- make depend depends on it


# DO NOT DELETE
