############# Makefile for project CFS++ ###############################
#
src = cfs.cc
#
C_DIR=..
#
########## set OS variable
#
OS=$(shell uname -s)
#

############ General options
#
include $(C_DIR)/main/Makefile.option
#

SUBLIBDIR := Driver PDE CoupledPDE Estimator\
             Forms Elements Elements/1D Elements/2D Elements/3D \
             Domain/GridCFS Domain \
             DataInOut/DefineFiles \
             DataInOut/Unverg DataInOut/GMV \
             DataInOut/DatFile DataInOut/AnsysFile DataInOut\
	     AlgebraicSystem/LinAlg\
             Utils Matrix General MpCCIcpl  

#
############ define my libraries
#

ifeq ($(MpCCI),yes) 
  mylibrary = driver domain pde mpccicpl coupledpde\
            estimator forms  \
            elements3d elements2d elements1d elements \
            gridcfs gmv \
            definefiles unverg ansysfile datainout \
            algsys-own\
            utils matrix general \
            multigrid-las 
else
  mylibrary = driver domain pde coupledpde\
            estimator forms  \
            elements3d elements2d elements1d elements \
            gridcfs \
            gmv \
            definefiles unverg ansysfile datainout \
            algsys-own\
            utils matrix general \
            multigrid-las
endif

#
ifeq ($(AdaptivLib),yes)
mylibrary = $(mylibrary) adaptgrid
endif


############ Directories for Adaptivity
#
ifeq ($AdaptivLib),yes)
SUBLIBDIR := $(SUBLIBDIR) Domain/AdaptGrid
endif

#
############ Define specific flags
#
include $(LIBSRC_DIR)/makefilemach/make.mach.$(MACHINE)
#

include $(C_DIR)/main/Makefile.option
ifeq ($(MpCCI),yes)
  CC=/opt/MpCCI/mpcci-2.0.0-linux-gcc32/bin/cciCC
  LINKFLAGS2 = -L/opt/MpCCI/mpcci-2.0.0-linux-gcc32/lib -L/opt/MPICH/mpich-1.2.5/lib
  sysliball = $(syslib) -lcci -lmpich -lg2c
else
  sysliball = $(syslib)
endif

##########
#
PATHSUBLIBDIR = $(addprefix $(LIBSRC_DIR)/, $(SUBLIBDIR))
#

########### Rules
#

.SUFFIXES: .f .c .cc .o

.cc.o:
	$(CC) -o $@ $(CFLAGS1) $(CFLAGS2) $(CFLAGS3) $(INCL) $<

#ifeq "$(MACHINE)" "SUN"
#$(program):: $(src:.cc=.o) $(mylibrary)
#	$(CC) $(CFLAG_EXE) -o $@ $(src:.cc=.o) $(fsrc:.f=.o) -L$(LIB_DIR) -L$(LIB_GRIDLIB)  -L$(LIB_NETGEN) -L$(LIB_CLA) $(LINKFLAGS2) $(LINKFLAGS3) $(INCL) $(addprefix -l, $(mylibrary)) $(syslib)
#else

ifeq ($(PETSc)),yes)
$(program):: $(src:.cc=.o) $(fsrc:.f=.o) $(mylibrary)
	$(CXX_CLINKER) -o $@ $(src:.cc=.o) $(fsrc:.f=.o) $(PETSC_LIB) \
	-L$(LIB_DIR) $(LINKFLAGS2) $(LINKFLAGS3) \
	-L$(LIB_DIR) -L$(LIB_GRIDLIB)  -L$(LIB_NETGEN) -L$(LIB_CLA) \
	$(INCL) $(addprefix -l, $(mylibrary)) $(sysliball) -Wl,-rpath,$(LIB_OWN)/$(MACHINE)
else
$(program):: $(src:.cc=.o) $(mylibrary)
	$(CC) $(CFLAG_EXE) -o $@ $(src:.cc=.o) $(fsrc:.f=.o) -L$(LIB_DIR) -L$(LIB_GRIDLIB)  -L$(LIB_NETGEN) -L$(LIB_CLA) $(LINKFLAGS2) $(LINKFLAGS3) $(INCL) $(addprefix -l, $(mylibrary)) $(sysliball) -Wl,-rpath,$(LIB_OWN)/$(MACHINE)
endif
#

########### Go to make library
#
include $(LIBSRC_DIR)/Makefile


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
	rm -rf $(LibCLA) && ln -s $(LAS)/libsrc $(LibCLA)
else
	rm -rf $(LibCLA) && ln -s $(LAS)/libsrc $(LibCLA)
endif
#
#ifeq "$(OS)" "Linux"
#	rm -rf $(IncNetGen) &&  ln -s $(NetGen)/libsrc/include $(IncNetGen) && \
#        rm -rf $(LibNetGen) && ln -s $(NetGen)/lib $(LibNetGen) &&\
#        rm -rf $(LibCLA) && ln -s $(CLA)/libsrc $(LibCLA) 
#else
#	rm -rf $(LibCLA) && ln -s $(CLA)/libsrc $(LibCLA) && \
#        rm -rf $(IncGridlib) && ln -s $(Gridlib) $(IncGridlib)
#endif
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
############## Make library adaptgrid
#
#adaptgrid:
#	@cd $(AdaptGrid) && $(MAKE) -f Makefile
clean-adaptgrid:
	@cd $(AdaptGrid) && $(MAKE) clean
#
######################  Clean work directory
clean-driver:
	@cd $(LIBSRC_DIR)/Driver && \
	$(MAKE) clean
clean-datainout:
	@cd $(LIBSRC_DIR)/DataInOut&& \
	$(MAKE) clean
clean-ansys:
	@cd $(LIBSRC_DIR)/DataInOut/AnsysFile && \
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
	@cd $(LIBSRC_DIR)/Domain/AdaptGrid&& \
        $(MAKE) clean
clean-general:
	@cd $(LIBSRC_DIR)/General&&\
	$(MAKE) clean
clean-pde:
	@cd $(LIBSRC_DIR)/PDE&&\
	$(MAKE) clean
clean-coupledpde:
	@cd $(LIBSRC_DIR)/CoupledPDE&&\
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
clean-mpccicpl:
	@cd $(LIBSRC_DIR)/MpCCIcpl&&\
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
######################   All clean --- for distributive

cleanall: cleanall-library cleanall-files
	@echo clean all
.PHONY: cleanall cleanall-library cleanall-files

cleanall-library:
	@for n in $(PATHSUBLIBDIR); do \
         (if cd $$n; \
                     then $(MAKE) allclean; \
                     else echo "WARNING:" $$n " doesn't exist"; \
          fi); \
         done 
cleanall-files:
	-rm *.o
#
######################   Clean library
# 
libclean: 
	@echo clean library
	@cd $(LIB_DIR) && \
	rm *
.PHONY: libclean
#
######################  Make documentation  
doc: doxygen

	@echo To see documentation run: konqueror  $(C_DIR)/doc/html/index.html 
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
##################### Make distributive
#
LIBSRC=$(addprefix libsrc/,$(SUBLIBDIR))
CVS_DIR=$(addsuffix /CVS,$(LIBSRC))
echo:
	echo $(CVS_DIR)
#
distrib: clean 
	cd $(C_DIR) && \
	tar -czvf cfs++.tgz $(LIBSRC) main/Makefile main/Makefile.option main/cfs.cc main/control 
#	tar --delete $(C_DIR) cfs++.tar 
#
##################### Make check
#
ControlDir = control
ControlExtOut = unverg
CheckFiles = channel 
#
check:
	@for n in $(CheckFiles); do\
	(cfs $$n ; \
	if diff $$n.$(ControlExtOut) $$n.test ;\
		then echo " $$n: Ok :) ";\
		else echo " $$n: Error or control files are not updated.";\
	fi);\
	done; 
#
##################### Make dependencies
#
depend: #depend-cleanfile
	@for n in $(PATHSUBLIBDIR); do \
         (if cd $$n; \
                     then $(MAKE) $@; \
                     else echo "WARNING:" $$n " doesn't exist"; \
          fi); 
	done
	@rm -f $(Dependencies).bak
#
depend-cleanfile:
	@rm -rf $(Dependencies)
 
# DO NOT DELETE THIS LINE -- make depend depends on it
