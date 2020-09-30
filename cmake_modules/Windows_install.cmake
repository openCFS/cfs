#-----------------------------------------------------------------------------
# Perform installation step of CFS
#-----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# A) Platform-dependent tasks 
# ----------------------------------------------------------------------------

IF(WIN32)
  # ----------------------------------------------------------------------------
  # 1) Copy shared libraries, exe and batch files to the binary directory 
  # ----------------------------------------------------------------------------
  INSTALL(DIRECTORY ${CFS_BINARY_DIR}/bin DESTINATION . 
	  #     PATTERN "*.exe" EXCLUDE
	  #	PATTERN "*.lib" EXCLUDE
	  # PATTERN "*.res" EXCLUDE
	  # PATTERN "*.exp" EXCLUDE
	  # PATTERN "*.ilk" EXCLUDE
	  # PATTERN "py*.bat" EXCLUDE
	  # PATTERN "*.py" EXCLUDE
	  # PATTERN "*.manifest" EXCLUDE
	  # PATTERN "*.sh" EXCLUDE
	  # PATTERN "*.txt" EXCLUDE
	  # PATTERN "*.xml" EXCLUDE
	  # PATTERN "adf2hdf.bat" EXCLUDE
	  # PATTERN "cgnsupdate.bat" EXCLUDE
	  # PATTERN "hdf2adf.bat" EXCLUDE
	  # PATTERN "xmdump*" EXCLUDE
	  )
          
  # ----------------------------------------------------------------------------
  # 2) Copy python libraries 
  # ----------------------------------------------------------------------------
  INSTALL(DIRECTORY ${CFS_BINARY_DIR}/lib64 DESTINATION . 
	  PATTERN "*-gd-*" EXCLUDE
	  )

  # ----------------------------------------------------------------------------
  # 3) Copy content of testsuite examples
  # ----------------------------------------------------------------------------
  INSTALL(DIRECTORY ${CFS_BINARY_DIR}/testsuite/TESTSUIT/Singlefield DESTINATION ./examples 
	  PATTERN "CMakeFiles" EXCLUDE
	  PATTERN "*.cmake" EXCLUDE
	  PATTERN "Makefile" EXCLUDE
	  )
  INSTALL(DIRECTORY ${CFS_BINARY_DIR}/testsuite/TESTSUIT/Coupledfield DESTINATION ./examples 
	  PATTERN "CMakeFiles" EXCLUDE
	  PATTERN "*.cmake" EXCLUDE
	  PATTERN "Makefile" EXCLUDE
	  )


ELSE(WIN32)

  # ----------------------------------------------------------------------------
  # 1) Copy shared libraries from $LIB_SUFFIX directory 
  # ----------------------------------------------------------------------------
  INSTALL(DIRECTORY ${CFS_BINARY_DIR}/${LIB_SUFFIX} DESTINATION . 
          PATTERN "*.a" EXCLUDE
          PATTERN "*.prl" EXCLUDE
          PATTERN "*.la" EXCLUDE
          PATTERN "*.sh" EXCLUDE
          PATTERN "libxslt-plugins" EXCLUDE
          PATTERN "pkgconfig" EXCLUDE
          PATTERN "libpng" EXCLUDE
          PATTERN "python2.7/test" EXCLUDE)
  
  # ----------------------------------------------------------------------------
  # 2) Copy qt-plugins directory to ${BIN_DIR} 
  # ----------------------------------------------------------------------------
  INSTALL(DIRECTORY ${CFS_BINARY_DIR}/${BIN_DIR}/plugins DESTINATION ${BIN_DIR})
ENDIF(WIN32)



# ----------------------------------------------------------------------------
# B) Platform-independent tasks 
# ----------------------------------------------------------------------------

# Run post-install script
#INSTALL(SCRIPT "cmake_modules/post-install.cmake")
