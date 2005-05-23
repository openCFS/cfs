###############################################################################
#                                                                             #
#                         Makefile for CFS++ project                          #
#                                                                             #
###############################################################################


# -----------------------------------------------------------------------------
# MAIN and EXECUTABLE
#
# Where to find the main routine and how to name the executable.
# -----------------------------------------------------------------------------
program = cfs
src = cfs.cc

# -----------------------------------------------------------------------------
# MACHINE-INDEPENDENT OPTIONS
#
# We include a file with machine-independent options.
# -----------------------------------------------------------------------------

include Makefile.option

# -----------------------------------------------------------------------------
# SUBDIRECTORIES
#
# We generate a list of all module subdirectories. The EXTRA_MODULES variable
# is composed in Makefile.option.
# -----------------------------------------------------------------------------

C_DIR = ..

SUBLIBDIR = Driver PDE CoupledPDE \
            Forms Elements Elements/1D Elements/2D Elements/3D \
            Domain/GridCFS Domain \
            DataInOut/DefineFiles \
            DataInOut/ParamHandling \
            DataInOut/Unverg DataInOut/GMV \
            DataInOut/AnsysFile DataInOut\
            DataInOut/Database\
            DataInOut/CommandLine\
            Utils Matrix General ODESolve ODEDescr \
	    ParamIdent\
            $(EXTRA_MODULES)

PATHSUBLIBDIR = $(addprefix $(LIBSRC_DIR)/, $(SUBLIBDIR))
TAGSPATH = $(addsuffix /*.??,$(PATHSUBLIBDIR))

# -----------------------------------------------------------------------------
# MYLIBRARIES
#
# We generate a list of all libraries that will be generate from the CFS++
# source code. The EXTRA_CFSLIBS variable is composed in Makefile.option.
# -----------------------------------------------------------------------------

mylibrary = driver domain pde coupledpde forms  \
            elements3d elements2d elements1d elements \
            gridcfs definefiles gmv unverg \
	    ansysfile datainout paramhandling commandline \
            driver \
            utils matrix general odesolve odedescr \
	    paramIdent \
            $(EXTRA_CFSLIBS)

# -----------------------------------------------------------------------------
# MACHINE-DEPENDENT OPTIONS
#
# We include a file with machine-dependent options. Actually MACHINE does not
# only refer to a certain architecture, but can also refer to a different
# compiler or a debug version.
# -----------------------------------------------------------------------------

include $(LIBSRC_DIR)/makefilemach/make.mach.$(MACHINE)

# -----------------------------------------------------------------------------
# GENERAL COMPILATION RULES
#
# We define a set of general compilation rules
# -----------------------------------------------------------------------------

.SUFFIXES: .f .c .cc .o

.cc.o:
	$(CC) -c -o $@ -I$(LIBSRC_DIR) $(CFLAGS_GLOBAL) $(CFLAGS_MACHINE) $<

# -----------------------------------------------------------------------------
# EXECUTABLE
#
# This is the default target and the one that will generate the CFS++
# executable. If you wonder what the :: means, see section "4.13 Double-Colon
# Rules" in the GNU make documentation.
# -----------------------------------------------------------------------------
$(program):: $(src:.cc=.o) $(mylibrary)
	$(LINK) $(LFLAGS_EXE) -o $@ $(src:.cc=.o) -L$(LIB_DIR) \
	$(LFLAGS_MACHINE) $(addprefix -l, $(mylibrary)) \
	$(LFLAGS_GLOBAL) $(syslib) 

# -----------------------------------------------------------------------------
# MODULE TARGETS
#
# We define a set of targets for the individual modules/subdirectories.
# -----------------------------------------------------------------------------
driver:
	@cd $(LIBSRC_DIR)/Driver; $(MAKE) -f Makefile
paramIdent:
	@cd $(LIBSRC_DIR)/ParamIdent; $(MAKE) -f Makefile 
domain:
	@cd $(LIBSRC_DIR)/Domain; $(MAKE) -f Makefile
gridcfs:
	@cd $(LIBSRC_DIR)/Domain/GridCFS; $(MAKE) -f Makefile
estimator:
	@cd $(LIBSRC_DIR)/Estimator; $(MAKE) -f Makefile
toadaptgrid:
	@cd $(LIBSRC_DIR)/Domain/AdaptGrid; $(MAKE) -f Makefile
adaptgrid:
	@cd $(AdaptGrid); $(MAKE) -f Makefile
forms:
	@cd $(LIBSRC_DIR)/Forms; $(MAKE) -f Makefile
pde:
	@cd $(LIBSRC_DIR)/PDE; $(MAKE) -f Makefile
coupledpde:
	@cd $(LIBSRC_DIR)/CoupledPDE; $(MAKE) -f Makefile
elements:
	@cd $(LIBSRC_DIR)/Elements; $(MAKE) -f Makefile
elements1d:
	@cd $(LIBSRC_DIR)/Elements/1D; $(MAKE) -f Makefile
elements2d:
	@cd $(LIBSRC_DIR)/Elements/2D; $(MAKE) -f Makefile
elements3d:
	@cd $(LIBSRC_DIR)/Elements/3D; $(MAKE) -f Makefile
definefiles:
	@cd $(LIBSRC_DIR)/DataInOut/DefineFiles; $(MAKE) -f Makefile
datainout:
	@cd $(LIBSRC_DIR)/DataInOut; $(MAKE) -f Makefile
ansysfile:
	@cd $(LIBSRC_DIR)/DataInOut/AnsysFile; $(MAKE) -f Makefile
ansys-out:
	@cd $(LIBSRC_DIR)/DataInOut/AnsysOut; $(MAKE) -f Makefile
unverg:
	@cd $(LIBSRC_DIR)/DataInOut/Unverg; $(MAKE) -f Makefile
gmv:
	@cd $(LIBSRC_DIR)/DataInOut/GMV; $(MAKE) -f Makefile
database:
	@cd $(LIBSRC_DIR)/DataInOut/Database; $(MAKE) -f Makefile
gsi:
	@cd $(LIBSRC_DIR)/DataInOut/GSI; $(MAKE) -f Makefile
utils:
	@cd $(LIBSRC_DIR)/Utils; $(MAKE) -f Makefile
matrix:
	@cd $(LIBSRC_DIR)/Matrix; $(MAKE) -f Makefile
general:
	@cd $(LIBSRC_DIR)/General; $(MAKE) -f Makefile
mpccicpl:
	@cd $(LIBSRC_DIR)/MpCCIcpl; $(MAKE) -f Makefile
paramhandling:
	@cd $(LIBSRC_DIR)/DataInOut/ParamHandling; $(MAKE) -f Makefile
commandline:
	@cd $(LIBSRC_DIR)/DataInOut/CommandLine; $(MAKE) -f Makefile
odesolve:
	@cd $(LIBSRC_DIR)/ODESolve; $(MAKE) -f Makefile
odedescr:
	@cd $(LIBSRC_DIR)/ODEDescr; $(MAKE) -f Makefile

# -----------------------------------------------------------------------------
# CLEANING TARGETS I
#
# We define a set of targets for cleaning individual modules/subdirectories or
# combinations thereof.
# -----------------------------------------------------------------------------

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
	@cd $(LIBSRC_DIR)/DataInOut/CommandLine && \
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
clean-odesolve:
	@cd $(LIBSRC_DIR)/ODESolve && \
	$(MAKE) clean
clean-odedescr:
	@cd $(LIBSRC_DIR)/ODEDescr && \
	$(MAKE) clean
clean-paramIdent:
	@cd $(LIBSRC_DIR)/ParamIdent && \
	$(MAKE) clean

# -----------------------------------------------------------------------------
# CLEANING TARGETS II
#
# We define a set of targets for cleaning most everything.
# -----------------------------------------------------------------------------

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

libclean: 
	@echo clean library
	@cd $(LIB_DIR) && \
	rm *
.PHONY: libclean

# -----------------------------------------------------------------------------
# DOCUMENTATION
#
# Targets related to the doxygen documentation for CFS++
# -----------------------------------------------------------------------------
.PHONY: doc clean-doc
doc:
	cd $(C_DIR)/Doc; ./makedoc.sh
clean-doc:
	@cd $(C_DIR)/Doc; rm -rf ./html

# -----------------------------------------------------------------------------
# DEPENDENCIES
#
# Automatic generation of dependencies resulting from #include directives
# with the makedepend tool
# -----------------------------------------------------------------------------
.Phony: depend clean-depend

depend:
	@for n in $(PATHSUBLIBDIR);do \
        (if cd $$n; \
                     then $(MAKE) $@;\
                     else echo "WARNING:" $$n " doesn't exist";\
        fi); done
	touch Makefile.depend
	@makedepend -I$(LIBSRC_DIR) $(INCDIR_OLAS)  -f$(DEPEND_FILE) $(DEPEND_FLAG) cfs.cc

clean-depend:
	@rm -f ./$(DEPEND_FILE)
	@rm -f ./$(DEPEND_FILE).bak
	@for n in $(PATHSUBLIBDIR);do \
        (if cd $$n; \
                     then $(MAKE) $@;\
                     else echo "WARNING:" $$n " doesn't exist";\
        fi); done

tags:
	etags --declarations --members --ignore-indentation  -o ../libsrc/TAGS $(TAGSPATH)

# DO NOT DELETE THIS LINE -- make depend depends on it
