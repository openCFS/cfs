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
# MODULES / SUBDIRECTORIES
#
# We generate a list of all module subdirectories. The EXTRA_MODULES variable
# is composed in Makefile.option.
# -----------------------------------------------------------------------------

C_DIR = ..

mymodules = Driver PDE CoupledPDE \
            Forms Elements Elements/1D Elements/2D Elements/3D \
            Domain/GridCFS Domain/GridAdaption Domain/AdaptGrid Domain \
            Domain/GridStruct \
            DataInOut/DefineFiles \
            DataInOut/ParamHandling \
            DataInOut/CommandLine \
            DataInOut/Unverg \
            DataInOut/cfsRST \
            DataInOut/GMV \
            DataInOut/AnsysFile \
            DataInOut/Scripting \
            DataInOut \
            DataInOut/Database\
            Utils Utils/mathParser Matrix General ODESolve ODEDescr \
	    ParamIdent Materials \
            $(EXTRA_MODULES)

PATHSUBLIBDIR = $(addprefix $(LIBSRC_DIR)/, $(mymodules))
TAGSPATH = $(addsuffix /*.??,$(PATHSUBLIBDIR))


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
.PHONY: $(program)
$(program):: $(src:.cc=.o) $(LIB_DIR)/lib$(CFS_LIB).a
	$(LINK) $(LFLAGS_EXE) -o $@ $(src:.cc=.o) -L$(LIB_DIR) \
	$(LFLAGS_MACHINE) -l$(CFS_LIB) \
	$(LFLAGS_GLOBAL) $(syslib) 


#
# This target triggers compilation of the CFS library
#
$(LIB_DIR)/lib$(CFS_LIB).a: $(mymodules)

# -----------------------------------------------------------------------------
# MODULE TARGETS
#
# We define a set of targets for the individual modules/subdirectories.
# -----------------------------------------------------------------------------
$(mymodules):
	@cd $(LIBSRC_DIR)/$@; $(MAKE) -f Makefile


# -----------------------------------------------------------------------------
# CLEANING TARGETS I
#
# We define a set of targets for cleaning individual modules/subdirectories or
# combinations thereof.
# -----------------------------------------------------------------------------

clean-targets = $(addprefix clean-, $(mymodules))

$(clean-targets):
	cd $(LIBSRC_DIR)/$(subst clean-,,$@) && $(MAKE) clean


# -----------------------------------------------------------------------------
# CLEANING TARGETS II
#
# We define a set of targets for cleaning most everything.
# -----------------------------------------------------------------------------

clean: clean-library clean-prog

.PHONY: clean clean-library clean-prog

clean-library:
	@for n in $(PATHSUBLIBDIR); do \
         (if cd $$n; \
                     then $(MAKE) clean; \
                     else echo "WARNING:" $$n " doesn't exist"; \
          fi); \
         done 
clean-prog:
	@-rm -f *.o

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
	cd $(C_DIR)/Doc/developer; ./makedoc.sh
clean-doc:
	@cd $(C_DIR)/Doc/developer; rm -rf ./html

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
