#
	export MACHINE=LINUXg
#
src = cfs.cc
#
C_DIR = ..
#
############ General options
#
include $(C_DIR)/main/option.make
#
############ Define specific flags
#
include $(LIBSRC_DIR)/makefilemach/make.mach.$(MACHINE)
#
############ Define library
#
mylibrary=  driver general algsys-own domain pde forms \
            elements algsys-own matrix utils datainout unverg \
            datfile datainout
#
syslib = -lm 
#
########### Rules
#
.SUFFIXES: .f .c .cc .o
#
.cc.o:
	$(CC) -o $@ $(CFLAGS1) $(CFLAGS2) $(CFLAGS3) $<
 
.c.o:
	$(CC) -o $@ $(CFLAGS1) $(CFLAGS2) $(CFLAGS3) $<
 
.f.o:
	$(F77) -o $@ $(FFLAGS1) $(FFLAGS2) $<
 
$(program):: $(src:.cc=.o) $(fsrc:.f=.o) makelibs
	$(CC) -o $@ $(src:.cc=.o) $(fsrc:.f=.o) -L$(LIB_DIR) $(syslib) $(LINKFLAGS2) $(LINKFLAGS3) $(addprefix -l, $(mylibrary)) $(LIB_NAG)
#
########### Go to make library
#
include $(LIBSRC_DIR)/Makefile
#
######################  Clean work directory
clean-driver:
	@cd $(LIBSRC_DIR)/Driver && \
	$(MAKE) clean
clean-datainout:
	@cd $(LIBSRC_DIR)/DataInOut&& \
	$(MAKE) clean
clean-dat:
	@cd $(LIBSRC_DIR)/DataInOut/DatFile && \
	$(MAKE) clean
clean-unverg:
	@cd $(LIBSRC_DIR)/DataInOut/Unverg && \
	$(MAKE) clean
clean-algebraicsys:
	@cd $(LIBSRC_DIR)/AlgebraicSystem/LinAlg && \
	$(MAKE) clean
clean-domain:
	@cd $(LIBSRC_DIR)/Domain&&\
	$(MAKE) clean
clean-general:
	@cd $(LIBSRC_DIR)/General&&\
	$(MAKE) clean
clean-pde:
	@cd $(LIBSRC_DIR)/PDE&&\
	$(MAKE) clean
clean-forms:
	@cd $(LIBSRC_DIR)/Forms&&\
	$(MAKE) clean
clean-elements:
	@cd $(LIBSRC_DIR)/Elements&&\
	$(MAKE) clean
clean-utils:
	@cd $(LIBSRC_DIR)/Utils&&\
	$(MAKE) clean
#clean-matrix
#	@cd $(LIBSRC_DIR)/Matrix&&\\
#	$(MAKE) clean

######################   Real clean

clean: clean-library clean-prog
	@echo clean all
.PHONY: clean clean-library clean-prog

clean-library: clean-driver clean-datainout clean-dat clean-unverg clean-algebraicsys clean-domain clean-general clean-pde clean-forms clean-elements clean-utils 
 
clean-prog:
	-rm *.o

######################   Clean library
 
libclean: 
	@echo clean library
	@cd $(LIB_DIR) && \
	rm *.a
.PHONY: libclean
 

######################  Make documentation  
doc: doxygen
	@echo To see documentation run: netscape  $(C_DIR)/doc/html/index.html 
doxygen:
	doxygen DoxygenConfig/doxy-config
FlagDoc= -nd -p -u -v -b -f 
doc++:
	@-mkdir $(C_DIR)/doc
	@-mkdir $(C_DIR)/doc/html
	@-rm $(C_DIR)/doc/html/*.*
	doc++ $(FlagDoc) -d $(C_DIR)/doc/online docfile
.PHONY: doc doxygen doc++
#
##################### Make dependencies
#
depend-driver:
	@cd $(LIBSRC_DIR)/Driver&& \
	$(MAKE) depend
depend-datainout:
	@cd $(LIBSRC_DIR)/DataInOut&& \
	$(MAKE) depend
depend-dat:
	@cd $(LIBSRC_DIR)/DataInOut/DatFile && \
	$(MAKE) depend
depend-unverg:
	@cd $(LIBSRC_DIR)/DataInOut/Unverg && \
	$(MAKE) depend
depend-algebraicsys:
	@cd $(LIBSRC_DIR)/AlgebraicSystem/LinAlg && \
	$(MAKE) depend
depend-domain:
	@cd $(LIBSRC_DIR)/Domain&&\
	$(MAKE) depend
depend-general:
	@cd $(LIBSRC_DIR)/General&&\
	$(MAKE) depend
depend-pde:
	@cd $(LIBSRC_DIR)/PDE&&\
	$(MAKE) depend
depend-forms:
	@cd $(LIBSRC_DIR)/Forms&&\
	$(MAKE) depend
depend-elements:
	@cd $(LIBSRC_DIR)/Elements&&\
	$(MAKE) depend 
depend-utils:
	@cd $(LIBSRC_DIR)/Utils&&\
	$(MAKE) depend
depend-matrix:
	@cd $(LIBSRC_DIR)/Matrix&&\\
	$(MAKE) depend

depend: depend-cleanfile depend-driver depend-datainout depend-dat depend-unverg depend-algebraicsys depend-domain depend-general depend-pde depend-forms depend-elements depend-utils
	rm -f $(Dependencies).bak
	@echo done
	@echo see file: $(Dependencies) 

depend-cleanfile:
	@rm -rf $(Dependencies)
 
# DO NOT DELETE THIS LINE -- make depend depends on it
