#
#
src = cfs.cc
#
C_DIR=..
########## set OS variable
#
OS=$(shell uname -s)
#
############ Directories
#
ifeq "$(OS)" "Linux"
SUBLIBDIR := Driver Estimator PDE \
             Forms Elements Elements/2D Elements/3D \
             Domain/GridCFS Domain/NetGen Domain \
             DataInOut/DefineFiles \
             DataInOut/Unverg DataInOut/GMV \
             DataInOut/DatFile DataInOut/AnsysFile DataInOut\
             Utils Matrix General
else
SUBLIBDIR := Driver Estimator PDE AlgebraicSystem/LinAlg \
            Forms Elements Elements/2D Elements/3D \
            Domain/GridCFS Domain/Gridlib Domain \
            DataInOut/DefineFiles \
            DataInOut/Unverg DataInOut/GMV \
            DataInOut/DatFile DataInOut\AnsysFile DataInOut\
            Utils Matrix General
endif
#
############ define my libraries
#
ifeq "$(OS)" "Linux"
mylibrary :=  driver estimator pde forms  \
            elements3d elements2d elements \
            gridcfs tonetgen\
            domain unverg gmv \
            definefiles datfile ansysfile datainout \
            utils matrix general \
            multigrid-cla general-cla superlu-cla \
            nginterface vis  \
            stlgeom mesh csg gprim opti \
            la gen cfsinterface tix8.1.8.3 tk8.3 tcl8.3\
            GL GLU X11 Xi m dl pthread stdc++
else
mylibrary :=  driver estimator pde algsys-own forms  \
            elements3d elements2d elements \
            gridcfs \
            domain unverg gmv \
            definefiles datfile ansysfile datainout \
            utils matrix \
            multigrid-cla general-cla superlu-cla \
            general
endif
#
########### define system libraries
#
ifeq "$(OS)" "Linux"
syslib :=  
else
syslib := -lGridBase -lIO  -lUtil -lMesh -lC
endif
#
############ General options
#
include $(C_DIR)/main/Makefile.option
#
############ Define specific flags
#
include $(LIBSRC_DIR)/makefilemach/make.mach.$(MACHINE)
#
##########
#
PATHSUBLIBDIR = $(addprefix $(LIBSRC_DIR)/, $(SUBLIBDIR))
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
	$(CC) $(CFLAG_EXE) -o $@ $(src:.cc=.o) $(fsrc:.f=.o) -L$(LIB_DIR) -L$(LIB_GRIDLIB)  -L$(LIB_NETGEN) -L$(LIB_CLA) $(LINKFLAGS2) $(LINKFLAGS3) $(addprefix -l, $(mylibrary)) $(syslib) $(LIB_NAG) 
#
########## Set varible MACHINE
setvar:
ifeq "$(OS)" "Linux"
	export MACHINE=LINUX
else
	@echo setenv MACHINE SGI
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
ifeq "$(OS)" "Linux"
	@echo Before run of cfs set this variable on LINUX:
	@echo export LD_LIBRARY_PATH=$(LIB_OWN)/$(MACHINE):$(NetGen)/lib
else
	@echo Before run of cfs set this variable on SGI:
	@echo setenv LD_LIBRARYN32_PATH $(LIB_GRIDLIB):$(LIB_OWN)/$(MACHINE)
endif
#
############## Make link to the source of another programm
#
link:
ifeq "$(OS)" "Linux"
	rm -rf $(IncNetGen) &&  ln -s $(NetGen)/libsrc/include $(IncNetGen) && \
        rm -rf $(LibNetGen) && ln -s $(NetGen)/lib $(LibNetGen) &&\
        rm -rf $(LibCLA) && ln -s $(CLA)/libsrc $(LibCLA) 
else
	rm -rf $(LibCLA) && ln -s $(CLA)/libsrc $(LibCLA) && \
        rm -rf $(IncGridlib) && ln -s $(Gridlib) $(IncGridlib)
endif
#        
############## Make links in dir include
#
make-include:
	@echo Creating links to header-files
	@(rm -rf $(Include); mkdir $(Include))
	@for n in $(PATHSUBLIBDIR); do \
        (if cd $$n; \
                  then $(MAKE) make-include ; \
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
clean-gmv:
	@cd $(LIBSRC_DIR)/DataInOut/GMV && \
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
	@cd $(LIBSRC_DIR)/Elements/2D&&\
        $(MAKE) clean
	@cd $(LIBSRC_DIR)/Elements/3D&&\
        $(MAKE) clean
clean-estimator:
	@cd $(LIBSRC_DIR)/Estimator && \
        $(MAKE) clean
clean-utils:
	@cd $(LIBSRC_DIR)/Utils&&\
	$(MAKE) clean
clean-matrix:
	@cd $(LIBSRC_DIR)/Matrix&&\
	$(MAKE) clean
clean-general-piles:
	@cd $(LIBSRC_DIR)/AlgebraicSystem/CLA/general && \
	$(MAKE) clean
clean-multigrid-piles:
	@cd $(LIBSRC_DIR)/AlgebraicSystem/CLA/multigrid && \
	$(MAKE) clean
clean-superlu-piles:
	@cd $(LIBSRC_DIR)/AlgebraicSystem/CLA/superlu && \
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
depend: #depend-cleanfile
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
