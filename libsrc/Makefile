#         Makefile for libsrc

what =

#amylibrary = $(addprefix a, $(mylibrary))
#makelibs:	$(what) $(amylibrary)
#$(amylibrary):
#cd $(LIBSRC_DIR)/$(patsubst a%,%,$@) && \
#$(MAKE) -f Makefile$(appl) $(what)

makelibs: $(what) adriver adomain aforms aelements ageneral apde autils aalgsys-own adatainout amatrix aunverg adatfile
	@echo done
adriver:
	@cd $(LIBSRC_DIR)/Driver; $(MAKE) -f Makefile $(appl) $(what)
adatainout:
	@cd $(LIBSRC_DIR)/DataInOut; $(MAKE) -f Makefile $(appl) $(what)
aalgsys-own:
	@cd $(LIBSRC_DIR)/AlgebraicSystem/LinAlg; $(MAKE) -f Makefile $(appl) $(what)
adomain:
	@cd $(LIBSRC_DIR)/Domain; $(MAKE) -f Makefile $(appl) $(what)
ageneral:
	@cd $(LIBSRC_DIR)/General; $(MAKE) -f Makefile $(appl) $(what)
apde:
	@cd $(LIBSRC_DIR)/PDE; $(MAKE) -f Makefile $(appl) $(what)
aforms:
	@cd $(LIBSRC_DIR)/Forms; $(MAKE) -f Makefile $(appl) $(what)
aelements:
	@cd $(LIBSRC_DIR)/Elements; $(MAKE) -f Makefile $(appl) $(what)
autils:
	@cd $(LIBSRC_DIR)/Utils; $(MAKE) -f Makefile $(appl) $(what)
adatfile:
	@cd $(LIBSRC_DIR)/DataInOut/DatFile; $(MAKE) -f Makefile $(appl) $(what)  
amatrix:
	@cd $(LIBSRC_DIR)/Matrix; $(MAKE) -f Makefile $(appl) $(what)  
aunverg:
	@cd $(LIBSRC_DIR)/DataInOut/Unverg; $(MAKE) -f Makefile $(appl) $(what)  
