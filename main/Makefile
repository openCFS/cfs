# ############ Makefile for project CFS++ ###############################
#
src = cfs.cc
# src = bubble.cc
#
C_DIR= ..
#
########## set OS variable
#
OS=$(shell uname -s)
#
############ General options
#
include $(C_DIR)/main/Makefile.option
#
SUBLIBDIR := Driver PDE CoupledPDE \
             Forms Elements Elements/1D Elements/2D Elements/3D \
             Domain/GridCFS Domain \
             DataInOut/DefineFiles \
             DataInOut/ParamHandling \
             DataInOut/Unverg DataInOut/GMV \
             DataInOut/DatFile DataInOut/AnsysFile DataInOut\
             DataInOut/Database\
             Utils Matrix General ODESolve ODEDescr
#
ifeq ($(Adaptivity),yes)
SUBLIBDIR := $(SUBLIBDIR) Domain/AdaptGrid Estimator
mylibrary := $(mylibrary) estimator toadaptgrid
endif
#
ifeq ($(MpCCI),yes)
SUBLIBDIR := $(SUBLIBDIR) MpCCIcpl
mylibrary := $(mylibrary) mpccicpl
endif
#
############ define my libraries
#
mylibrary = driver domain pde coupledpde forms  \
            elements3d elements2d elements1d elements \
            gridcfs definefiles gmv unverg\
	    ansysfile datainout paramhandling \
            driver \
            utils matrix general  odesolve odedescr
#
ifeq "$(LASTYPE)" "OLAS"

else
mylibrary += multigrid-las
endif
#
#
ifeq ($(Adaptivity),yes)
mylibrary := $(mylibrary) adaptgrid 
endif
#
ifeq ($(MpCCI),yes)
SUBLIBDIR := $(SUBLIBDIR) MpCCIcpl
mylibrary := $(mylibrary) mpccicpl
endif
#
############ Define specific flags
#
include $(LIBSRC_DIR)/makefilemach/make.mach.$(MACHINE)
#
#
############ General options
#
include $(C_DIR)/main/Makefile.option
#

#####DataBase
ifeq ($(DATABASE),yes)
  syslibDB    =  -lmysqlclient
  LINKFLAGSDB =  -L/usr/lib/mysql 
  SUBLIBDIR  :=  $(SUBLIBDIR) DataInOut/Database
  mylibrary  :=  $(mylibrary) database
endif


#####GSI format
ifeq ($(GSI),yes)
  SUBLIBDIR  :=  $(SUBLIBDIR) DataInOut/GSI
  mylibrary  :=  $(mylibrary) gsi
endif

############ MpCCI
#
ifeq ($(MpCCI),yes)
  CC=/opt/MpCCI/mpcci-3.0.0-sdk-linux-glibc23-gcc33/bin/cciCC
  LINKFLAGS2 = -L/opt/MpCCI/mpcci-3.0.0-sdk-linux-glibc23-gcc33/lib -L/opt/mpich-1.2.5.2/lib $(LINKFLAGSDB)
  sysliball = $(syslib) -lcci -lmpich $(XERCESLIB) $(syslibDB) -lg2c
else
  sysliball = $(syslib) $(XERCESLIB) $(OLASLIB) $(syslibDB)
  LINKFLAGS2 = $(LINKFLAGSDB)
endif


#
##########
#
PATHSUBLIBDIR = $(addprefix $(LIBSRC_DIR)/, $(SUBLIBDIR))
#

########### Rules
#

.SUFFIXES: .f .c .cc .o

.cc.o:
	$(CC) -o $@ $(CFLAGS1) $(CFLAGS2) $(CFLAGS3) $(OLAS) $(TRACEFLAGS) $(INCL) $(INCLAdaptGrid) $(INCLMpCCI) $(INCLDATABASE) $<

#
#
#
ifeq ($(PETSc)),yes)
$(program):: $(src:.cc=.o) $(fsrc:.f=.o) $(mylibrary)
	$(CXX_CLINKER) -o $@ $(src:.cc=.o) $(fsrc:.f=.o) $(PETSC_LIB) \
	-L$(LIB_DIR) $(LINKFLAGS2) $(LINKFLAGS3) \
	-L$(CVSEXTERNAL)/$(MACHINE)/lib -L$(LIB_GRIDLIB) -L$(LIB_NETGEN) \
	-L$(LIB_CLA) $(INCL) $(addprefix -l, $(mylibrary)) $(sysliball) \
	-Wl,-rpath,$(LIB_OWN)/$(MACHINE) 
else
$(program):: $(src:.cc=.o) $(mylibrary)
	$(CC) $(CFLAG_EXE) -o $@ $(src:.cc=.o) $(fsrc:.f=.o) -L$(LIB_DIR) \
	-L$(LIB_GRIDLIB) -L$(LIB_NETGEN) -L$(LIB_CLA) \
	-L$(CVSEXTERNAL)/$(MACHINE)/lib \
	$(INCL) $(LINKFLAGS2) $(LINKFLAGS3) $(INCL) \
	$(addprefix -l, $(mylibrary)) $(sysliball) 
#-Wl,-rpath,$(LIB_OWN)/$(MACHINE)
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
clean-adaptgrid:
	@cd $(AdaptGrid) && $(MAKE) clean
clean-driver:
	@cd $(LIBSRC_DIR)/Driver && \
	$(MAKE) clean
clean-datainout:
	@cd $(LIBSRC_DIR)/DataInOut&& \
	$(MAKE) clean
	@cd $(LIBSRC_DIR)/DataInOut/DefineFiles && \
	$(MAKE) clean
	@cd $(LIBSRC_DIR)/DataInOut/ParamHandling && \
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
clean-database:
	@cd $(LIBSRC_DIR)/DataInOut/Database && \
	$(MAKE) clean
clean-gsi:
	@cd $(LIBSRC_DIR)/DataInOut/GSI && \
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
clean-newpde:
	$(MAKE) clean-pde
	$(MAKE) clean-driver
	$(MAKE) clean-domain
clean-odesolve:
	@cd $(LIBSRC_DIR)/ODESolve && \
	$(MAKE) clean
clean-odedescr:
	@cd $(LIBSRC_DIR)/ODEDescr && \
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
depend: 
	@for n in $(PATHSUBLIBDIR);do \
        (if cd $$n; \
                     then $(MAKE) $@;\
                     else echo "WARNING:" $$n " doesn't exist";\
        fi); done
	makedepend $(INCL) -fMakefile.depend cfs.cc

clean-depend:
	@for n in $(PATHSUBLIBDIR);do \
        (if cd $$n; \
                     then $(MAKE) $@;\
                     else echo "WARNING:" $$n " doesn't exist";\
        fi); done

# DO NOT DELETE THIS LINE -- make depend depends on it
