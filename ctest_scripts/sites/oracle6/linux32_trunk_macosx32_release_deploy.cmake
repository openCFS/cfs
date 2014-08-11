#-----------------------------------------------------------------------------
# Set binary directory
#-----------------------------------------------------------------------------
SET(CTEST_BINARY_DIRECTORY "$ENV{HOME}/Documents/dev/CFS_BUILD_NIGHTLY")
SET(ENV{PATH} "/opt/pckg/macosx_10.6_dev/bin:$ENV{PATH}")

message("  Searching for cfs and cfstool executables...")
SET(ARCH "I386")
SET(CFS_BIN_DIR "${CTEST_BINARY_DIRECTORY}/bin/MACOSX_10.6_${ARCH}")
SET(CFS_BINARY "${CFS_BIN_DIR}/cfsbin")
SET(CFSTOOL_BINARY "${CFS_BIN_DIR}/cfstoolbin")

MESSAGE("CFS_BINARY ${CFS_BINARY}")
MESSAGE("CFSTOOL_BINARY ${CFSTOOL_BINARY}")

IF(EXISTS "${CFS_BINARY}" AND
   EXISTS "${CFSTOOL_BINARY}")

  message("  Cleaning up binary directory...")

  # Remove unnecessary subdirs before packing archive.
  FILE(REMOVE_RECURSE
    "${CTEST_BINARY_DIRECTORY}/cfsdeps"
    "${CTEST_BINARY_DIRECTORY}/include"
    "${CTEST_BINARY_DIRECTORY}/tmp"
    "${CTEST_BINARY_DIRECTORY}/testsuite"
    "${CTEST_BINARY_DIRECTORY}/source"
    "${CTEST_BINARY_DIRECTORY}/CMakeFiles"
    "${CTEST_BINARY_DIRECTORY}/Testing"
  )

  file(GLOB_RECURSE STATIC_ARCHIVES "${CTEST_BINARY_DIRECTORY}/lib/*.a")
  FILE(REMOVE ${STATIC_ARCHIVES})

  message("  Packing binaries...")

  FILE(REMOVE_RECURSE "$ENV{HOME}/Documents/dev/deploy")
  FILE(MAKE_DIRECTORY "$ENV{HOME}/Documents/dev/deploy")

  EXECUTE_PROCESS(
    COMMAND mv CFS_BUILD_NIGHTLY deploy
    WORKING_DIRECTORY "$ENV{HOME}/Documents/dev"
    RESULT_VARIABLE RETVAL
  )
  
  EXECUTE_PROCESS(
    COMMAND zip -yr CFS_BUILD_NIGHTLY.zip CFS_BUILD_NIGHTLY
    WORKING_DIRECTORY "$ENV{HOME}/Documents/dev/deploy"
    RESULT_VARIABLE RETVAL
  )

  FILE(REMOVE_RECURSE "$ENV{HOME}/Documents/dev/deploy/CFS_BUILD_NIGHTLY")

  EXECUTE_PROCESS(
    COMMAND date +%F
    OUTPUT_VARIABLE DATE
    RESULT_VARIABLE RETVAL
  )
  STRING(REPLACE "\n"  " " DATE ${DATE})
  STRING(STRIP ${DATE} DATE)

  CONFIGURE_FILE("/vagrant/readme_macosx.md"
                 "$ENV{HOME}/Documents/dev/deploy/readme_macosx.md"
                 @ONLY)
  EXECUTE_PROCESS(
    COMMAND pandoc --highlight-style tango -s -S -H "/vagrant/readme_macosx.css" "readme_macosx.md" -o README.html
    WORKING_DIRECTORY "$ENV{HOME}/Documents/dev/deploy"
    RESULT_VARIABLE RETVAL
  )
  FILE(REMOVE "$ENV{HOME}/Documents/dev/deploy/readme_macosx.md")

  FILE(MAKE_DIRECTORY "$ENV{HOME}/Documents/dev/deploy/example")
  EXECUTE_PROCESS(
    COMMAND ${CMAKE_COMMAND} -E copy_directory "$ENV{HOME}/Documents/dev/NIGHTLY/CFS_TESTSUITE_FESPACE/TESTSUIT/Singlefield/Electrostatics/Cube3d" example
    WORKING_DIRECTORY "$ENV{HOME}/Documents/dev/deploy"
    RESULT_VARIABLE RETVAL
  )
  FILE(REMOVE_RECURSE "$ENV{HOME}/Documents/dev/deploy/example/.svn")
  FILE(REMOVE "$ENV{HOME}/Documents/dev/deploy/example/CMakeLists.txt")
  FILE(READ "$ENV{HOME}/Documents/dev/deploy/example/Cube3d.xml" CUBE3D_XML)
  FILE(REMOVE "$ENV{HOME}/Documents/dev/deploy/example/Cube3d.xml")
  STRING(REPLACE "<hdf5/>" "<gmsh/>" EXAMPLE_XML ${CUBE3D_XML})
  STRING(REPLACE "</pdeList>" "</pdeList>
    <linearSystems>
      <system>
        <solverList>
          <pardiso>
            <logging>yes</logging>
            <stats>yes</stats>
          </pardiso>
        </solverList>
      </system>
    </linearSystems>" EXAMPLE_XML ${EXAMPLE_XML})
  FILE(WRITE "$ENV{HOME}/Documents/dev/deploy/example/example.xml" ${EXAMPLE_XML})

  EXECUTE_PROCESS(
    COMMAND genisoimage -D -V "CFS_NIGHTLY_${ARCH}_${DATE}" -no-pad -r -apple -o "${HOSTNAME}_${TEST_NAME}.dmg" deploy
    WORKING_DIRECTORY "$ENV{HOME}/Documents/dev"
    RESULT_VARIABLE RETVAL
  )

  EXECUTE_PROCESS(
    COMMAND dmg dmg 
      "${HOSTNAME}_${TEST_NAME}.dmg"
      "${NIGHTLY_ARCHIVES_DIR}/${HOSTNAME}_${TEST_NAME}.dmg"
    WORKING_DIRECTORY "$ENV{HOME}/Documents/dev"
    RESULT_VARIABLE RETVAL
  )

  FILE(REMOVE "$ENV{HOME}/Documents/dev/${HOSTNAME}_${TEST_NAME}.dmg")
  FILE(REMOVE_RECURSE "$ENV{HOME}/Documents/dev/deploy")

ENDIF()
