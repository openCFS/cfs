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
.SUFFIXES: .f .cc .o

.cc.o:
	$(CC) $(CFLAGS1) $(CFLAGS2) $(CFLAGS3) $< 
#
libgoal: $(LIB_DIR) $(LIBB)
#
# compile each file seperately
#	
#LIBB = $(LIB_DIR)/lib$(lib).so
#
#ifeq "$(TYPEOFLIB)" "STATIC"
#
#$(LIBB): $(LIBB).a
#
#$(LIBB).a: $(src:.cc=.o)
#	$(AR) $(ARFLAGS) $@ $(?:.cc=.o) 
#	i$(RANLIB) $@ 
#
#else

#$(LIBB): $(LIBB).so

$(LIBB): $(src:.cc=.o)
	$(CC) $(FLAG_MKSHLIB) -o $@ $(?:.cc=.o) 

#endif
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
depend: 
	makedepend $(INCL)  -f$(C_DIR)/main/$(Dependencies) $(FlagMakeDepend) $(src) 

# DO NOT DELETE THIS LINE -- make depend depends on it
