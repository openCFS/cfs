#
# make an include file for library
#
include $(C_DIR)/main/option.make
#
include $(LIBSRC_DIR)/makefilemach/make.mach.$(MACHINE)
#
ARFLAGS = r 
#
LIBB = $(LIB_DIR)/lib$(lib).a
#
.SUFFIXES: .f .cc .o

.cc.o:
	$(CC) $(CFLAGS1) $(CFLAGS2) $(CFLAGS3) $< 
#
libgoal: $(LIB_DIR) $(LIBB)
#
# compile each file seperately
#	
$(LIBB): $(src:.cc=.o)
	$(AR) $(ARFLAGS) $@ $(?:.cc=.o) 
	-$(RANLIB) $@ 
#
#
$(LIB_DIR) :
	-@mkdir $(C_DIR)/lib
	@mkdir $(LIB_DIR)
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
