#
	export MACHINE=LINUXg
#
src = cfs.cc
#
C_DIR = ..
LIBDIR = $(C_DIR)/libsrc
SUBLIBDIR = General DataInOut DataInOut/DatFile \
          DataInOut/Unverg AlgebraicSystem/LinAlg \
          Domain  Driver Elements Forms \
          PDE Utils Matrix
PATHSUBLIBDIR = $(addprefix $(LIBDIR)/, $(SUBLIBDIR))
#
############ General options
#
include $(C_DIR)/main/Makefile.option
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
#include $(LIBSRC_DIR)/Makefile
#
makelibs:
	@for n in $(PATHSUBLIBDIR); do \
         (if cd $$n; \
                     then $(MAKE) -f Makefile; \
                     else echo "WARNING:" $$n " doesn't exist"; \
          fi); \
         done
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
clean-matrix:
	@cd $(LIBSRC_DIR)/Matrix&&\\
	$(MAKE) clean

######################   Real clean

clean: clean-library clean-prog
	@echo clean all
.PHONY: clean clean-library clean-prog

clean-library:
	@for n in $(PATHSUBLIBDIR); do \
         (if cd $$n; \
                     then $(MAKE) clean; \
                     else echo "WARNING:" $$n " doesn't exist"; \
          fi); \
         done 
clean-prog:
	-rm *.o
#
######################   Clean library
# 
libclean: 
	@echo clean library
	@cd $(LIB_DIR) && \
	rm *.a
.PHONY: libclean
#
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
depend: depend-cleanfile
	@for n in $(PATHSUBLIBDIR); do \
         (if cd $$n; \
                     then $(MAKE) $@; \
                     else echo "WARNING:" $$n " doesn't exist"; \
          fi); \
         done
	@rm -f $(Dependencies).bak
#
depend-cleanfile:
	@rm -rf $(Dependencies)
 
# DO NOT DELETE THIS LINE -- make depend depends on it
