#
#
src = cfs.cc
#
############ Directories
C_DIR= ..
#
LIBDIR = $(C_DIR)/libsrc
#
SUBLIBDIR = Driver PDE AlgebraicSystem/LinAlg \
            Forms Elements Domain/GridCFS Domain \
            DataInOut/DefineFiles \
            DataInOut/Unverg DataInOut/DatFile DataInOut\
            Utils Matrix General
PATHSUBLIBDIR = $(addprefix $(LIBDIR)/, $(SUBLIBDIR))
#
LIB_GRIDLIB = $(C_DIR)/libgrid
LIB_OWN = $(C_DIR)/lib
#
#
############
ifeq "$(ARCH)" "LINUX"
syslib :=
else
syslib := -lGridBase -lIO  -lUtil -lMesh -lC
endif
############ General options
#
include $(C_DIR)/main/Makefile.option
#
############ Define specific flags
#
include $(LIBSRC_DIR)/makefilemach/make.mach.$(MACHINE)
#
########### Rules
#
.SUFFIXES: .f .c .cc .o
#
.cc.o:
	$(CC) $(CFLAG_EXE) -o $@ $(CFLAGS1) $(CFLAGS2) $(CFLAGS3) $<
 
.c.o:
	$(CC) -o $@ $(CFLAGS1) $(CFLAGS2) $(CFLAGS3) $<
 
.f.o:
	$(F77) -o $@ $(FFLAGS1) $(FFLAGS2) $<
 
$(program):: $(src:.cc=.o) $(fsrc:.f=.o) makelibs 
	$(CC) $(CFLAG_EXE) -o $@ $(src:.cc=.o) $(fsrc:.f=.o) -L$(LIB_DIR) -L$(LIB_GRIDLIB)  $(LINKFLAGS2) $(LINKFLAGS3) $(addprefix -l, $(mylibrary)) $(syslib) $(LIB_NAG) 
#
########## Set varible MACHINE
setvar:
ifeq "$(ARCH)" "SGI"
	@echo setenv MACHINE SGI
else
	export MACHINE=LINUX
endif
########### Go to make library
#
makelibs:
	@for n in $(PATHSUBLIBDIR); do \
         (if cd $$n; \
                     then $(MAKE) -f Makefile; \
                     else echo "WARNING:" $$n " doesn't exist"; \
          fi); \
         done
#
#
############### Print info about variables that should be set before work
#
info: 
ifeq "$(ARCH)" "LINUX"
	@echo Before run of cfs set this variable on LINUX:
	@echo export LD_LIBRARY_PATH=$(LIB_OWN)/$(MACHINE)
else
	@echo Before run of cfs set this variable on SGI:
	@echo setenv LD_LIBRARYN32_PATH $(LIB_GRIDLIB):$(LIB_OWN)/$(MACHINE)
endif
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
clean-algsys:
	@cd $(LIBSRC_DIR)/AlgebraicSystem/LinAlg && \
	$(MAKE) clean
clean-domain:
	@cd $(LIBSRC_DIR)/Domain&&\
	$(MAKE) clean
	@cd $(LIBSRC_DIR)/Domain/GridCFS&& \
	$(MAKE) clean
	@cd $(LIBSRC_DIR)/Domain/Gridlib&& \
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
	rm *.so
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
