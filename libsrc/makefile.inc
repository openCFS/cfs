#
# make an include file for library
#
include $(C_DIR)/main/Makefile.option
#


include $(LIBSRC_DIR)/makefilemach/make.mach.$(MACHINE)

#
LIBB := $(LIB_DIR)/lib$(CFS_LIB).a 

#
.SUFFIXES: .f .cc .cpp .o
.cc.o:
	$(CC) -c -I$(LIBSRC_DIR) $(CFLAGS_GLOBAL) $(CFLAGS_MACHINE) $< 
.cpp.o:
	$(CC) -c -I$(LIBSRC_DIR) $(CFLAGS_GLOBAL) $(CFLAGS_MACHINE) $< 

#
libgoal: $(LIB_DIR) $(LIBB)
#
#################### to make library
#


.PHONY: link
#$(LIBB): 
.PRECIOUS: $(LIBB)
$(LIBB): $(src:.cc=.o) $(src:.cpp=.o)
	 $(AR) $(ARFLAGS) $@ $(addsuffix .o, $(basename $(src)))
	-@$(RANLIB) $@

$(LIB_DIR) :
	-@mkdir $(C_DIR)/lib
	@mkdir $(LIB_DIR)

link:
	 $(AR) $(ARFLAGS) $(LIBB) $(addsuffix .o, $(basename $(src)))
	-@$(RANLIB) $(LIBB)
#
########## Make clean for each directory and its subdirectories
#
clean:
	@for f in '$(shell find -name "*.o")'; do \
	(rm -f $$f; \
	$(AR) -d $(LIBB) $$f); done

allclean: clean
	@for f in '$(shell find -name "*.*~")'; do \
	(rm -f $$f); done
#

######### Make depend 
#
depend:
	@touch $(DEPEND_FILE)
	makedepend -I$(LIBSRC_DIR) $(INCDIR_OLAS) -f$(DEPEND_FILE) $(DEPEND_FLAG) $(src)

clean-depend:
	@rm -f $(DEPEND_FILE)
	@rm -f $(DEPEND_FILE).bak

# DO NOT DELETE THIS LINE -- make depend depends on it
