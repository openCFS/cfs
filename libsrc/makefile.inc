#
# make an include file for library
#
include $(C_DIR)/main/Makefile.option
#
ifeq "$(AdaptivLib)" "yes"
	ADAPT = -DADAPTGRID
endif

include $(LIBSRC_DIR)/makefilemach/make.mach.$(MACHINE)
#
ARFLAGS = r 
#
ifeq "$(LIBRARY_TYPE)" "shared"
	LIBB := $(LIB_DIR)/lib$(lib).so
else
	LIBB := $(LIB_DIR)/lib$(lib).a 
endif

ifeq "$(PETSc)" "yes"
	CFLAGSPETSc = -D USE_PETSC
endif

#
.SUFFIXES: .f .cc .cpp .o
.cc.o:
	$(CC) $(CFLAGS1) $(CFLAGS2) $(CFLAGS3) $(CFLAGSPETS) $(INCL) $(ADAPT) $< 

.cpp.o:
	$(CC) $(CFLAGS1) $(CFLAGS2) $(CFLAGS3) $(CFLAGSPETS) $(INCL) $<


#
libgoal: $(LIB_DIR) $(LIBB)
#
#################### to make library
#
ifeq "$(LIBRARY_TYPE)" "shared"

$(LIBB): $(src:.cc=.o)
	$(CC) $(FLAG_MKSHLIB) -o $@ $(?:.cc=.o) 
else
$(LIBB): $(src:.cc=.o)
	@$(AR) $(ARFLAGS) $@ $(?:.cc=.o); 
#	-@rm $(?:.cc=.o);  
	-@$(RANLIB) $@
endif
#
$(LIB_DIR) :
	-@mkdir $(C_DIR)/lib
	@mkdir $(LIB_DIR)
#
########## Make clean for each directory
#
clean:
	-rm -f *.o && rm -f $(LIBB);
allclean:
	-rm *.o && -rm *~; 
#
########## Make link to include
#
make-include: $(src:.cc=.hh) $(src_h) 
	@for n in $(src:.cc=.hh); do \
        (ln -sf $$n $(Include)/$$n) ;\
        done
#
######### Make depend 
#
depend: 
	makedepend $(INCL)  -f$(C_DIR)/main/$(Dependencies) $(FlagMakeDepend) $(src) 

# DO NOT DELETE THIS LINE -- make depend depends on it
