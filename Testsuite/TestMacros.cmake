#-------------------------------------------------------------------------------
# Emit the CPU model name as a CDash measurement so it appears in the test
# result table (next to Processors). Linux: lscpu. Windows: in-process registry
# CPU_MODEL e.g. "Intel(R) Xeon(R) Platinum 8581C CPU @ 2.10GHz", "AMD EPYC 7F32 8-Core Processor"
#-------------------------------------------------------------------------------
macro(emit_cpu_dart)
  set(CPU_MODEL "")
  if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
    execute_process(
      COMMAND sh -c "lscpu | sed -n 's/^Model name:[[:space:]]*//p'"
      OUTPUT_VARIABLE CPU_MODEL OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  elseif(CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")
    # the alternative is > 1 sec: powershell -NoProfile -Command "(Get-CimInstance Win32_Processor | Select-Object -First 1).Name"
    cmake_host_system_information(RESULT CPU_MODEL
      QUERY WINDOWS_REGISTRY "HKLM/HARDWARE/DESCRIPTION/System/CentralProcessor/0"
      VALUE "ProcessorNameString")
  endif()
  if(CPU_MODEL)
    execute_process(COMMAND ${CMAKE_COMMAND} -E echo
      "<DartMeasurement name=\"CPU\" type=\"text/string\">${CPU_MODEL}</DartMeasurement>")
  endif()
endmacro()

#-------------------------------------------------------------------------------
# Define a macro to generate the test name from the directory name given by
# SRCDIR.
#-------------------------------------------------------------------------------
MACRO(GENERATE_TEST_NAME_AND_FILE SRCDIR)
#  MESSAGE("SRCDIR ${SRCDIR}")
  STRING(REGEX MATCH "(TESTSUIT|PYTHON)/.*$" CURRENT_TEST_SUBDIR "${SRCDIR}")
  # MESSAGE("CURRENT_TEST_SUBDIR ${CURRENT_TEST_SUBDIR}")
  STRING(REGEX REPLACE "TESTSUIT/" "" CURRENT_TEST_SUBSUBDIR "${CURRENT_TEST_SUBDIR}")
  # MESSAGE("CURRENT_TEST_SUBSUBDIR ${CURRENT_TEST_SUBSUBDIR}")
  STRING(REGEX REPLACE "/" "_" TEST_NAME "${CURRENT_TEST_SUBSUBDIR}")
  # MESSAGE("TEST_NAME ${TEST_NAME}")
  STRING(REGEX MATCH "[^/]+$" TEST_FILE_BASENAME "${SRCDIR}")
ENDMACRO(GENERATE_TEST_NAME_AND_FILE)

#-------------------------------------------------------------------------------
# Determine a list of files which are necessary to run a test. Some input files
# need to be copied to TESTSUITE_BIN_DIR afterwards by COPY_TEST_FILES.
#-------------------------------------------------------------------------------
MACRO(DETERMINE_TEST_FILE_LIST)
  FILE(GLOB FILES "${CURRENT_TEST_DIR}/*")

  SET(TEST_FILE_LIST "")

  # Disregard temporary and other files.
  FOREACH(file ${FILES})
    # pattern to ignore
    IF(${file} MATCHES "[.]h5ref$" OR
       ${file} MATCHES ".+/[.].+$" OR
       ${file} MATCHES "~" OR
       ${file} MATCHES "#" OR
       ${file} MATCHES "[.]las$" OR
       (${file} MATCHES "[.]mtx$" AND NOT
        ${file} MATCHES "[.]ref.mtx$") OR
       ${file} MATCHES "[.]vec$" OR
       ${file} MATCHES "[.]pdf$" OR
       ${file} MATCHES "CMakeLists.txt$" OR
       ${file} MATCHES "results_" OR
       ${file} MATCHES "history" OR
       ${file} MATCHES "ansys_" OR
       ${file} MATCHES "gid_" OR
       ${file} MATCHES "gmsh_" OR
       ${file} MATCHES "[.]geo$" OR
       ${file} MATCHES "[.]png$")
      MESSAGE("ignoring ${file}")
    ELSE()
      SET(TEST_FILE_LIST "${file};${TEST_FILE_LIST}")
    ENDIF()
  ENDFOREACH()

  # Append files contained in global variable ADD_TEST_FILES
  FOREACH(F ${ADD_TEST_FILES})
    SET(TEST_FILE_LIST "${CURRENT_TEST_DIR}/${F};${TEST_FILE_LIST}") 
  ENDFOREACH()
  
  # Determine name of .h5ref
  FILE(GLOB H5REF_FILE "${CURRENT_TEST_DIR}/${TEST_FILE_BASENAME}.h5ref")
ENDMACRO(DETERMINE_TEST_FILE_LIST)


#-------------------------------------------------------------------------------
# Copy files to TESTSUITE_BIN_DIR. 
#-------------------------------------------------------------------------------
MACRO(COPY_TEST_FILES FILE_LIST SKIP_H5REF)
  FILE(MAKE_DIRECTORY "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}")
  FOREACH(F ${FILE_LIST})
    FILE(COPY "${F}" DESTINATION "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}")
  ENDFOREACH()

  FILE(REMOVE "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/${TEST_FILE_BASENAME}.h5ref")
  
  #IF(NOT (${SKIP_H5REF} STREQUAL \"ON\"))
  IF(NOT SKIP_H5REF)
	  IF(UNIX)
		EXECUTE_PROCESS(
		  COMMAND "${CMAKE_COMMAND}" -E create_symlink "${H5REF_FILE}" "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/${TEST_FILE_BASENAME}.h5ref"
		  WORKING_DIRECTORY "."
		  RESULT_VARIABLE RETVAL)
		 # MESSAGE(FATAL_ERROR "On Unix")
	  ELSE()
		EXECUTE_PROCESS(
		  COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${H5REF_FILE}" "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}"
		  WORKING_DIRECTORY "."
		  RESULT_VARIABLE RETVAL)
		 # MESSAGE(FATAL_ERROR "Not on Unix")
	  ENDIF()

	  IF(NOT RETVAL EQUAL 0)
		MESSAGE(FATAL_ERROR "Could not create symlink to .h5ref file.")
	  ENDIF()
  ELSE()
	SET(RETVAL 0)
  ENDIF()
ENDMACRO(COPY_TEST_FILES)




#-------------------------------------------------------------------------------
# Remove remnants of previous runs.
#-------------------------------------------------------------------------------
MACRO(CLEANUP_TEST_DIR)
  SET(FILE_LIST "results_hdf5;${TEST_FILE_BASENAME}.info.xml;${TEST_FILE_BASENAME}.las")

  FOREACH(F ${FILE_LIST})
    # MESSAGE("F ${F}")
    FILE(REMOVE_RECURSE ${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/${F})
  ENDFOREACH()
ENDMACRO(CLEANUP_TEST_DIR)

# Removed to reduce complexity: Fabian, 23.1.2016
# MACRO(CHECK_NUMBER_OF_NODES_AND_ELEMS H5REF)


#-------------------------------------------------------------------------------
# Run CFS++.
#-------------------------------------------------------------------------------
MACRO(RUN_TEST_SIMULATION)
  EXECUTE_PROCESS(
    COMMAND "${CFS_BINARY}" ${CFS_ARGS} --noColor "${TEST_FILE_BASENAME}"
    WORKING_DIRECTORY "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}"
    ERROR_VARIABLE CFS_ERROR
    ECHO_ERROR_VARIABLE          # also stream cfs stderr (warnings!) to the test log, not only capture it for the crash message
    RESULT_VARIABLE CFS_RETVAL)
  IF(NOT CFS_RETVAL EQUAL 0)
    MESSAGE("ERROR: RV=${CFS_RETVAL} : ${CFS_ERROR}")
    message("COMMAND = ${CFS_BINARY} ${CFS_ARGS} --noColor ${TEST_FILE_BASENAME}")
    message("WORKING_DIRECTORY = ${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}")
    if(PROCEED_AFTER_SIMULATION_CRASH)
      message(WARNING "CFS simulation for test case '${TEST_NAME}' failed. Continuing anyways ... (PROCEED_AFTER_SIMULATION_CRASH=${PROCEED_AFTER_SIMULATION_CRASH})")
    else()
      message(FATAL_ERROR "CFS simulation for test case '${TEST_NAME}' failed.")
    endif()
  ENDIF()
ENDMACRO(RUN_TEST_SIMULATION)

MACRO(RUN_MPI_TEST_SIMULATION)
#  EXECUTE_PROCESS(COMMAND nproc OUTPUT_VARIABLE NPROC)
  EXECUTE_PROCESS(
    COMMAND mpiexec -n 4 "${CFS_BINARY}" ${CFS_ARGS} --noColor "${TEST_FILE_BASENAME}"
    WORKING_DIRECTORY "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}"
    ERROR_VARIABLE CFS_ERROR
    RESULT_VARIABLE CFS_RETVAL)
  IF(NOT CFS_RETVAL EQUAL 0)
    MESSAGE("ERROR: RV=${CFS_RETVAL} : ${CFS_ERROR}")
    MESSAGE(FATAL_ERROR "CFS++ simulation for test case '${TEST_NAME}' failed.")
  ENDIF()
ENDMACRO(RUN_MPI_TEST_SIMULATION)

#-------------------------------------------------------------------------------
# Run a preCICE-coupled test: launch the mock partner (MOCK_BINARY) and cfs
# concurrently in the test working directory. EXECUTE_PROCESS runs several
# COMMANDs as a parallel pipeline and waits for all of them; RESULTS_VARIABLE
# yields one exit code per command, so a crash of *either* fails the test.
# The CTest TIMEOUT property guards against a hung m2n socket.
#-------------------------------------------------------------------------------
MACRO(RUN_COUPLED_SIMULATION)
  SET(_COUPLED_WD "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}")
  # preCICE leaves an exchange directory behind; start clean
  FILE(REMOVE_RECURSE "${_COUPLED_WD}/precice-run")
  EXECUTE_PROCESS(
    COMMAND "${MOCK_BINARY}" "${PRECICE_CONFIG}" "${MOCK_SPEC}"
    COMMAND "${CFS_BINARY}" ${CFS_ARGS} --noColor "${TEST_FILE_BASENAME}"
    WORKING_DIRECTORY "${_COUPLED_WD}"
    RESULTS_VARIABLE COUPLED_RETVALS
    ERROR_VARIABLE COUPLED_ERROR)
  FOREACH(_rv ${COUPLED_RETVALS})
    IF(NOT _rv EQUAL 0)
      MESSAGE("ERROR: exit codes=${COUPLED_RETVALS} : ${COUPLED_ERROR}")
      MESSAGE("COMMAND = ${MOCK_BINARY} ${PRECICE_CONFIG} ${MOCK_SPEC} | ${CFS_BINARY} ${CFS_ARGS} --noColor ${TEST_FILE_BASENAME}")
      MESSAGE("WORKING_DIRECTORY = ${_COUPLED_WD}")
      IF(PROCEED_AFTER_SIMULATION_CRASH)
        MESSAGE(WARNING "coupled run for test case '${TEST_NAME}' failed. Continuing anyways ...")
      ELSE()
        MESSAGE(FATAL_ERROR "coupled run for test case '${TEST_NAME}' failed.")
      ENDIF()
    ENDIF()
  ENDFOREACH()
ENDMACRO(RUN_COUPLED_SIMULATION)

#-------------------------------------------------------------------------------
# Like RUN_COUPLED_SIMULATION, but the coupling partner is a SECOND openCFS
# instance (cfs+cfs) instead of the mock_fluid binary. Used for openCFS<->openCFS
# preCICE tests (e.g. the two-domain characteristic acoustic coupling). The partner
# participant's problem name is PARTNER_BASENAME (<PARTNER_BASENAME>.xml); the main
# participant (compared against the .h5ref) is TEST_FILE_BASENAME. Both run in the
# same working directory and share the preCICE exchange directory.
#-------------------------------------------------------------------------------
MACRO(RUN_COUPLED_CFS_SIMULATION)
  SET(_COUPLED_WD "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}")
  # preCICE leaves an exchange directory behind; start clean
  FILE(REMOVE_RECURSE "${_COUPLED_WD}/precice-run")
  # IMPORTANT: execute_process chains multiple COMMANDs as a PIPELINE
  # (partner stdout -> main stdin). Both cfs participants are verbose, and the main
  # cfs never reads stdin, so the partner would block once the ~64 KB pipe buffer
  # fills -> deadlock (the run "hangs" under ctest even though each works manually).
  # Redirect the partner's output to a log file via `sh` so the pipe carries nothing;
  # the two processes still start concurrently (required for the preCICE handshake).
  # `exec` makes the partner's exit status propagate as the sh exit status.
  # (preCICE socket-coupled tests target Unix.)
  SET(_PARTNER_CFS_ARGS "${CFS_ARGS}")
  STRING(REPLACE ";" " " _PARTNER_CFS_ARGS "${_PARTNER_CFS_ARGS}")
  EXECUTE_PROCESS(
    COMMAND sh -c "exec \"${CFS_BINARY}\" ${_PARTNER_CFS_ARGS} --noColor \"${PARTNER_BASENAME}\" > \"${PARTNER_BASENAME}.log\" 2>&1"
    COMMAND "${CFS_BINARY}" ${CFS_ARGS} --noColor "${TEST_FILE_BASENAME}"
    WORKING_DIRECTORY "${_COUPLED_WD}"
    RESULTS_VARIABLE COUPLED_RETVALS
    ERROR_VARIABLE COUPLED_ERROR)
  FOREACH(_rv ${COUPLED_RETVALS})
    IF(NOT _rv EQUAL 0)
      MESSAGE("ERROR: exit codes=${COUPLED_RETVALS} : ${COUPLED_ERROR}")
      MESSAGE("COMMAND = ${CFS_BINARY} ${CFS_ARGS} --noColor ${PARTNER_BASENAME} (-> ${PARTNER_BASENAME}.log)  +  ${CFS_BINARY} ${CFS_ARGS} --noColor ${TEST_FILE_BASENAME}")
      MESSAGE("WORKING_DIRECTORY = ${_COUPLED_WD} (partner output in ${PARTNER_BASENAME}.log)")
      IF(PROCEED_AFTER_SIMULATION_CRASH)
        MESSAGE(WARNING "coupled cfs+cfs run for test case '${TEST_NAME}' failed. Continuing anyways ...")
      ELSE()
        MESSAGE(FATAL_ERROR "coupled cfs+cfs run for test case '${TEST_NAME}' failed.")
      ENDIF()
    ENDIF()
  ENDFOREACH()
ENDMACRO(RUN_COUPLED_CFS_SIMULATION)

#-------------------------------------------------------------------------------
# Like RUN_COUPLED_SIMULATION, but the coupling partner is a real OpenFOAM solver
# using the cfsdeps OpenFOAM + openfoam-adapter (requires OPENFOAM_ENV_SCRIPT,
# the openfoam-env.sh of the cfs build dir). The partner runs from the SAME
# working directory as cfs via `-case ${FOAM_CASE_DIR}`, so the (relative)
# preCICE m2n exchange directory is identical for both participants. blockMesh
# generates the fluid mesh at test time (no polyMesh files in the repo). All
# partner output goes to log files (see the pipe-deadlock note in
# RUN_COUPLED_CFS_SIMULATION); `exec` propagates the solver's exit status.
#-------------------------------------------------------------------------------
MACRO(RUN_COUPLED_FOAM_SIMULATION)
  SET(_COUPLED_WD "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}")
  # preCICE leaves an exchange directory behind; start clean
  FILE(REMOVE_RECURSE "${_COUPLED_WD}/precice-run")
  EXECUTE_PROCESS(
    COMMAND bash -c "source '${OPENFOAM_ENV_SCRIPT}' > log.openfoam-env 2>&1 && blockMesh -case '${FOAM_CASE_DIR}' > log.blockMesh 2>&1 && exec ${FOAM_SOLVER} -case '${FOAM_CASE_DIR}' > log.${FOAM_SOLVER} 2>&1"
    COMMAND "${CFS_BINARY}" ${CFS_ARGS} --noColor "${TEST_FILE_BASENAME}"
    WORKING_DIRECTORY "${_COUPLED_WD}"
    RESULTS_VARIABLE COUPLED_RETVALS
    ERROR_VARIABLE COUPLED_ERROR)
  FOREACH(_rv ${COUPLED_RETVALS})
    IF(NOT _rv EQUAL 0)
      MESSAGE("ERROR: exit codes=${COUPLED_RETVALS} : ${COUPLED_ERROR}")
      MESSAGE("COMMAND = ${FOAM_SOLVER} -case ${FOAM_CASE_DIR} (logs: log.openfoam-env, log.blockMesh, log.${FOAM_SOLVER})  +  ${CFS_BINARY} ${CFS_ARGS} --noColor ${TEST_FILE_BASENAME}")
      MESSAGE("WORKING_DIRECTORY = ${_COUPLED_WD}")
      IF(PROCEED_AFTER_SIMULATION_CRASH)
        MESSAGE(WARNING "coupled cfs+OpenFOAM run for test case '${TEST_NAME}' failed. Continuing anyways ...")
      ELSE()
        MESSAGE(FATAL_ERROR "coupled cfs+OpenFOAM run for test case '${TEST_NAME}' failed.")
      ENDIF()
    ENDIF()
  ENDFOREACH()
ENDMACRO(RUN_COUPLED_FOAM_SIMULATION)

#-------------------------------------------------------------------------------
# Run cfstool in specified mode.
#-------------------------------------------------------------------------------
MACRO(DIFF_TEST_RESULTS_CFSTOOL EPSILON)
  
  set(h5_res "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/results_hdf5/${TEST_FILE_BASENAME}")
  IF(NOT EXISTS "${h5_res}.h5")
    if(EXISTS "${h5_res}.cfs") # must be renamed because cfstool exepts a file called *.h5
      message("found '${h5_res}.cfs' renaming to '${h5_res}.h5'")
      file(RENAME "${h5_res}.cfs" "${h5_res}.h5")
    else()
      MESSAGE("No current h5 file found! Test will fail.")
    endif(EXISTS "${h5_res}.cfs")
  ENDIF()
  message("TESTSUITE_CFSTOOL_MODE=${TESTSUITE_CFSTOOL_MODE}")
  message("Used CFSTOOL_MODE=${CFSTOOL_MODE}")
  # I strongly dissagree with the paternalism to limit the number of cells for the testsuite.
  # I'm old enough to judge myself! Therefore I set checkLimits to false. Fabian, 25.10.2017.
  EXECUTE_PROCESS(
    COMMAND "${CFSTOOL_BINARY}" -m "${CFSTOOL_MODE}" --checkLimits false --eps "${EPSILON}" "${H5REF_FILE}" "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/results_hdf5/${TEST_FILE_BASENAME}.h5"
    WORKING_DIRECTORY "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}"
    ERROR_VARIABLE CFSTOOL_ERROR
    RESULT_VARIABLE CFSTOOL_RETVAL
  )

  # Since EXECUTE_PROCESS does not respect the /B(atch) flag to exit in some 
  # versions of Windows CMD.exe we explicitly search for an EXITCODE text in
  # stderr of cfstool.bat.
  IF(WIN32 AND CFSTOOL_ERROR)
    STRING(REPLACE "\n" ";" CFSTOOL_ERROR ${CFSTOOL_ERROR})

    FOREACH(LINE IN ITEMS ${CFSTOOL_ERROR})
      MESSAGE("LINE ${LINE}")
      IF(LINE MATCHES "EXITCODE=")
        SET(CFSTOOL_RETVAL 1)
      ENDIF()
    ENDFOREACH()

    MESSAGE("CFSTOOL_BINARY ${CFSTOOL_BINARY}")
    MESSAGE("CFSTOOL_RETVAL ${CFSTOOL_RETVAL}")
  ENDIF()

  IF(NOT CFSTOOL_RETVAL EQUAL 0)
    MESSAGE("executed command: ${CFSTOOL_BINARY} -m ${CFSTOOL_MODE} --checkLimits false --eps ${EPSILON} ${H5REF_FILE} ${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/results_hdf5/${TEST_FILE_BASENAME}.h5")
    MESSAGE(WARNING "${CFSTOOL_ERROR}")
    MESSAGE(FATAL_ERROR "Results for test case '${TEST_NAME}' could not be validated.")
  ENDIF()
ENDMACRO(DIFF_TEST_RESULTS_CFSTOOL)

#-------------------------------------------------------------------------------
# Compare result .h5 against .h5ref using h5diff (from hdf5-tools).
#
# Tolerance variables (set via -D in the calling test's CMakeLists.txt):
#
#   EPSILON   relative tolerance, mapped to h5diff -p.
#             A pair (a, b) is equal if  |a - b| / |b|  <=  EPSILON.
#             Same semantic as cfstool's --eps in relL2diff mode, so the
#             convention carries over from existing CFS tests.
#             Use for FP arithmetic output (postprocessing, IFFTs, etc.).
#
#   ABSTOL    absolute tolerance, mapped to h5diff -d.
#             A pair (a, b) is equal if  |a - b|  <=  ABSTOL.
#             Useful as a floor next to EPSILON when reference values can
#             be zero (relative comparison breaks down at 0).
#
#   When both EPSILON and ABSTOL are set, h5diff treats a pair as equal if
#   *either* criterion holds — the standard "absolute OR relative" rule.
#   When neither is set, comparison is bit-exact (correct for deterministic
#   geometry-only outputs like the extrudeMesh tests).
#
#   H5DIFF_OBJ   optional Object within the HDF5 file to restrict comparison
#                 to, e.g. "/Mesh". Empty = compare the whole file.
#-------------------------------------------------------------------------------
MACRO(DIFF_TEST_RESULTS_H5DIFF)
  set(h5_res "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/results_hdf5/${TEST_FILE_BASENAME}.h5")
  if(NOT EXISTS "${h5_res}")
    message(FATAL_ERROR "Python test produced no '${h5_res}' — nothing to diff.")
  endif()
  if(NOT H5DIFF_EXECUTABLE)
    message(FATAL_ERROR "h5diff not found. Install hdf5-tools and reconfigure.")
  endif()

  set(H5DIFF_ARGS "")

  # /FileInfo holds the write timestamp, creator, and version of the writing
  # tool — never useful to compare and the Date dataset will always differ
  # between the .h5ref (written once) and a fresh test result (written today).
  list(APPEND H5DIFF_ARGS "--exclude-path" "/FileInfo")

  # Additional user-supplied exclusions (semicolon-separated list).
  foreach(p ${H5DIFF_EXCLUDE})
    list(APPEND H5DIFF_ARGS "--exclude-path" "${p}")
  endforeach()

  if(EPSILON)
    list(APPEND H5DIFF_ARGS "-p" "${EPSILON}")
  endif()
  if(ABSTOL)
    list(APPEND H5DIFF_ARGS "-d" "${ABSTOL}")
  endif()
  list(APPEND H5DIFF_ARGS "${H5REF_FILE}" "${h5_res}")
  if(H5DIFF_OBJ)
    string(REGEX REPLACE "^\"|\"$" "" H5DIFF_OBJ "${H5DIFF_OBJ}")
    list(APPEND H5DIFF_ARGS "${H5DIFF_OBJ}")
  endif()

  message(STATUS "Running: ${H5DIFF_EXECUTABLE} ${H5DIFF_ARGS} in ${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}.")
  execute_process(
    COMMAND "${H5DIFF_EXECUTABLE}" ${H5DIFF_ARGS}
    WORKING_DIRECTORY "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}"
    OUTPUT_VARIABLE H5DIFF_OUTPUT
    ERROR_VARIABLE  H5DIFF_ERROR
    RESULT_VARIABLE H5DIFF_RETVAL
  )

  if(NOT H5DIFF_RETVAL EQUAL 0)
    message("executed: ${H5DIFF_EXECUTABLE} ${H5DIFF_ARGS}")
    message("${H5DIFF_OUTPUT}")
    if(H5DIFF_ERROR)
      message("${H5DIFF_ERROR}")
    endif()
    message(FATAL_ERROR
      "h5diff: results for '${TEST_NAME}' differ from reference "
      "(exit ${H5DIFF_RETVAL}).")
  else()
    message(STATUS "h5diff: results for '${TEST_NAME}' match reference.")
  endif()
ENDMACRO(DIFF_TEST_RESULTS_H5DIFF)

#-------------------------------------------------------------------------------
# Run compare_info_xml.py which searches for special stuff to compare in info.xml
# SKIP_NOISE: optional
# LAST: optional see Optimization/MultiSequence/eigen_stiffness for an example. 
#-------------------------------------------------------------------------------
macro(DIFF_TEST_RESULTS_INFO_XML EPSILON SKIP_NOISE LAST)
  message(STATUS "compare selected data from info.xml via ${COMPARE_INFO_XML} using ${Python_EXECUTABLE}") 
  # have no "" around LAST in the command, otherwise command fails when no last is given
  message(STATUS "WORKING_DIRECTORY ${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}")  
  execute_process(
    COMMAND ${Python_EXECUTABLE} "${COMPARE_INFO_XML}" "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/${TEST_FILE_BASENAME}.ref.info.xml" "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/${TEST_FILE_BASENAME}.info.xml" --eps ${EPSILON} --skip_noise ${SKIP_NOISE} ${LAST} 
    WORKING_DIRECTORY "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}"
    ERROR_VARIABLE COMPARE_INFO_XML_ERROR
    RESULT_VARIABLE COMPARE_INFO_XML_RETVAL
  )
  
  if(NOT COMPARE_INFO_XML_RETVAL EQUAL 0)
    message(STATUS "CMD ${Python_EXECUTABLE} ${COMPARE_INFO_XML} ${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/${TEST_FILE_BASENAME}.ref.info.xml ${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/${TEST_FILE_BASENAME}.info.xml --eps ${EPSILON} --skip_noise ${NOISE_OPT} ${LAST} ")
    message(WARNING "COMPARE_INFO_XML_RETVAL: ${COMPARE_INFO_XML_RETVAL}")
    message(WARNING "COMPARE_INFO_XML_ERROR: ${COMPARE_INFO_XML_ERROR}")
    message(FATAL_ERROR "Results for test case '${TEST_NAME}' could not be validated.")
  endif()
endmacro(DIFF_TEST_RESULTS_INFO_XML)

#-------------------------------------------------------------------------------
# Run compare_hist_file.py compares the *.hist (text) file of the test output 
# to the equally named *.histref file in the test folder 
#-------------------------------------------------------------------------------
macro(DIFF_TEST_RESULTS_HIST_FILE EPSILON)
  message(STATUS "compare output .hist file to the provided .hist.ref file via ${COMPARE_HIST_FILE}\n") 
  message(STATUS "WORKING_DIRECTORY ${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}")
  execute_process(
    COMMAND ${Python_EXECUTABLE} "${COMPARE_HIST_FILE}" "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/${TEST_FILE_BASENAME}.hist.ref" "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/${TEST_FILE_BASENAME}.hist" --eps ${EPSILON}
    WORKING_DIRECTORY "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}"
    ERROR_VARIABLE COMPARE_HIST_FILE_ERROR
    RESULT_VARIABLE COMPARE_HIST_FILE_RETVAL
  )
  
  if(NOT COMPARE_HIST_FILE_RETVAL EQUAL 0)
    message(STATUS "CMD ${Python_EXECUTABLE} ${COMPARE_HIST_FILE} ${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/${TEST_FILE_BASENAME}.hist.ref ${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/${TEST_FILE_BASENAME}.hist --eps ${EPSILON}")
    message(WARNING "COMPARE_HIST_FILE_RETVAL: ${COMPARE_HIST_FILE_RETVAL}")
    message(WARNING "COMPARE_HIST_FILE_ERROR: ${COMPARE_HIST_FILE_ERROR}")
    message(FATAL_ERROR "compare_hist_file.py failed for ${TEST_NAME}.")
  endif()
endmacro(DIFF_TEST_RESULTS_HIST_FILE)

#-------------------------------------------------------------------------------
# Run python tests
# For mativz: run matviz.py and call compare_image.py afterwards
# This is called via PythonTest.cmake(.in) Only variables defined there
# or the the test-case call are known here
#-------------------------------------------------------------------------------
MACRO(RUN_AND_TEST_PYTHON EPSILON NO_COMPARE)
  MESSAGE("perform python tests\n") 

  IF (MATVIZ_ARGS)
    RUN_AND_TEST_MATVIZ()
  ENDIF()
  
  ######################## basecell###############
  IF (BASECELL_ARGS AND NOT CFS_HOMOGENIZE_ARGS)
    RUN_AND_TEST_BASECELL()
  ENDIF()
  IF (CFS_HOMOGENIZE_ARGS)
    RUN_AND_TEST_BASECELL_HOMOGENIZATION()
  ENDIF()
  
  ###################### interpretation and evaluation with CFS++ ###############
  if (INTERPRETATION_ARGS)
    if(THREADS)
      message("set CFS, OMP and MKL threads to ${THREADS}")
      set(ENV{CFS_NUM_THREADS} "${THREADS}")
      set(ENV{MKL_NUM_THREADS} "${THREADS}")
      set(ENV{OMP_NUM_THREADS} "${THREADS}")
    endif()
    RUN_TEST_SIMULATION()
    if(TEST_INFO_XML) # set in the testcase via -DTEST_INFO_XML:STRING="ON"
      DIFF_TEST_RESULTS_INFO_XML("${EPSILON}" "${SKIP_NOISE}" "${LAST}")
    elseif(TEST_HIST_FILE)
      DIFF_TEST_RESULTS_HIST_FILE("${EPSILON}")
    else()
      # Run cfstool to diff against the .h5ref.
      DIFF_TEST_RESULTS_CFSTOOL("${EPSILON}")
    endif()
  endif()
  
  if(PYTHON_ARGS)
    COPY_TEST_FILES("${TEST_FILE_LIST}" "OFF")
    separate_arguments(NO_QUOTES UNIX_COMMAND "${PYTHON_ARGS}")
    separate_arguments(ARGS_LIST UNIX_COMMAND "${NO_QUOTES}")
    EXECUTE_PROCESS(
      # create the hdf5_results folder
      COMMAND ${CMAKE_COMMAND} -E make_directory "${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/results_hdf5"
      # important to have to quotation marks around the variable
      COMMAND ${Python_EXECUTABLE} ${ARGS_LIST}
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      ERROR_VARIABLE PYTHON_ERROR
      RESULT_VARIABLE PYTHON_RETVAL
    )
    if(NOT PYTHON_RETVAL EQUAL 0)
      message("CMD ${Python_EXECUTABLE} ${ARGS_LIST}")
      message(WARNING "PYTHON_RETVAL: ${PYTHON_RETVAL}")
      message(WARNING "PYTHON_ERROR: ${PYTHON_ERROR}")
      message(FATAL_ERROR "Python test failed for ${TEST_NAME}.")
    endif()
    if(TEST_H5DIFF) # set in the testcase via -DTEST_H5DIFF:BOOL=ON
      DIFF_TEST_RESULTS_H5DIFF()
    else()
      message(STATUS "skip h5diff on request by -DTEST_H5DIFF=OFF")
    endif()
  endif()
ENDMACRO(RUN_AND_TEST_PYTHON)

MACRO(RUN_AND_TEST_MATVIZ)
   # remove quotes
  separate_arguments(MATVIZ_NO_QUOTES UNIX_COMMAND "${MATVIZ_ARGS}")
  # create a list
  separate_arguments(MATVIZ_ARGS_LIST UNIX_COMMAND "${MATVIZ_NO_QUOTES}")

  # perform the matviz test
  EXECUTE_PROCESS(
    # important to have to quotation marks around the variable
    COMMAND ${Python_EXECUTABLE} "${MATVIZ_PY}" ${MATVIZ_ARGS_LIST}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    ERROR_VARIABLE MATVIZ_ERROR
    RESULT_VARIABLE MATVIZ_RETVAL
  )
  
  IF(NOT MATVIZ_RETVAL EQUAL 0)
    MESSAGE("MATVIZ_ERROR = ${MATVIZ_ERROR}")
    MESSAGE("MATVIZ_RETVAL = ${MATVIZ_RETVAL}")
    message("CMD: ${Python_EXECUTABLE} ${MATVIZ_PY} ${MATVIZ_ARGS_LIST}")
    MESSAGE(FATAL_ERROR "Matviz tests failed for ${TEST_NAME}.")
  ENDIF()

  # compare the matviz result with the reference image
  IF(NO_COMPARE)
    MESSAGE("skip comparison on request by NO_COMPARE")
  ELSE()  
    # if not set the argument is not given for compare_images
    IF(EPSILON)
      SET(EPS_TEST "--eps" "${EPSILON}")
    ENDIF()
  
    EXECUTE_PROCESS(
      COMMAND ${Python_EXECUTABLE} "${COMPARE_IMAGES}" ${EPS_TEST} ${TEST}.ref.png ${TEST}.png
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      ERROR_VARIABLE COMPARE_IMAGES_ERROR
      RESULT_VARIABLE COMPARE_IMAGES_RETVAL
    )
  
    IF(NOT COMPARE_IMAGES_RETVAL EQUAL 0)
      message("CMD: ${Python_EXECUTABLE} ${COMPARE_IMAGES} ${EPS_TEST} ${TEST}.ref.png ${TEST}.png")
      MESSAGE("COMPARE_IMAGES_ERROR = ${COMPARE_IMAGES_ERROR}")
      MESSAGE("COMPARE_IMAGES_RETVAL = ${COMPARE_IMAGES_RETVAL}")
      MESSAGE(FATAL_ERROR "compare_images.py failed for ${TEST_NAME}.")
    ENDIF()
  ENDIF()
ENDMACRO(RUN_AND_TEST_MATVIZ)

MACRO(RUN_AND_TEST_BASECELL)
  # remove quotes
  separate_arguments(BASECELL_NO_QUOTES UNIX_COMMAND "${BASECELL_ARGS}")
  # create a list
  separate_arguments(BASECELL_ARGS_LIST UNIX_COMMAND "${BASECELL_NO_QUOTES}")
  
  # perform basecell tests
  EXECUTE_PROCESS(
    # important to have to quotation marks around the variable
    COMMAND  ${Python_EXECUTABLE} "${BASECELL_PY}" ${BASECELL_ARGS_LIST}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    ERROR_VARIABLE BASECELL_ERROR
    RESULT_VARIABLE BASECELL_RETVAL
  )
  
  IF(NOT BASECELL_RETVAL EQUAL 0)
    MESSAGE("CMD = ${Python_EXECUTABLE} ${BASECELL_PY} ${BASECELL_ARGS_LIST}")
    MESSAGE("BASECELL_ERROR = ${BASECELL_ERROR}")
    MESSAGE("BASECELL_RETVAL = ${BASECELL_RETVAL}")
    MESSAGE(FATAL_ERROR "Basecell tests failed for ${TEST_NAME}.")
  ENDIF()
  
  IF(EPSILON)
    SET(EPS_TEST "--eps" "${EPSILON}")
  ENDIF()
  
  DIFF_TEST_RESULTS_INFO_XML("${EPSILON}" "${SKIP_NOISE}" "${LAST}")  
ENDMACRO(RUN_AND_TEST_BASECELL)

MACRO(RUN_AND_TEST_BASECELL_HOMOGENIZATION)
  # remove quotes
  separate_arguments(BASECELL_NO_QUOTES UNIX_COMMAND "${BASECELL_ARGS}")
  # create a list
  separate_arguments(BASECELL_ARGS_LIST UNIX_COMMAND "${BASECELL_NO_QUOTES}")
  
  # perform basecell tests
  EXECUTE_PROCESS(
    # important to have to quotation marks around the variable
    COMMAND ${Python_EXECUTABLE} "${BASECELL_PY}" ${BASECELL_ARGS_LIST}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    ERROR_VARIABLE BASECELL_ERROR
    RESULT_VARIABLE BASECELL_RETVAL
  )  
    
  IF(NOT BASECELL_RETVAL EQUAL 0)
    message("CMD = ${Python_EXECUTABLE} ${BASECELL_PY} ${BASECELL_ARGS_LIST}")
    MESSAGE("BASECELL_ERROR = ${BASECELL_ERROR}")
    MESSAGE("BASECELL_RETVAL = ${BASECELL_RETVAL}")
    MESSAGE(FATAL_ERROR "Basecell tests failed for ${TEST_NAME}.")
  ENDIF()
  
  # remove quotes
  separate_arguments(CFS_NO_QUOTES UNIX_COMMAND "${CFS_HOMOGENIZE_ARGS}")
  # create a list
  separate_arguments(CFS_ARGS_LIST UNIX_COMMAND "${CFS_NO_QUOTES}")
  
  EXECUTE_PROCESS(
    COMMAND "${CFS_BINARY}" --noColor ${CFS_ARGS_LIST}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    ERROR_VARIABLE CFS_ERROR
    RESULT_VARIABLE CFS_RETVAL)
    
  IF(NOT CFS_RETVAL EQUAL 0)
    MESSAGE("CFS_ERROR = ${CFS_ERROR}")
    MESSAGE("CFS_RETVAL = ${CFS_RETVAL}")
    MESSAGE(FATAL_ERROR "Basecell tests failed for ${TEST_NAME}.")
  ENDIF()  
  
  IF(EPSILON)
    SET(EPS_TEST "--eps" "${EPSILON}")
  ENDIF()
  
  SET(TEST_FILE_BASENAME "homogenize")
  EXECUTE_PROCESS(
    COMMAND ${Python_EXECUTABLE} "${COMPARE_INFO_XML}" "${TEST}_homogenize.ref.info.xml" "${TEST}_homogenize.info.xml" --eps ${EPSILON} ${LAST}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    ERROR_VARIABLE COMPARE_INFO_XML_ERROR
    RESULT_VARIABLE COMPARE_INFO_XML_RETVAL
  )
  
  IF(NOT COMPARE_INFO_XML_RETVAL EQUAL 0)
    message("CMD = ${Python_EXECUTABLE} ${COMPARE_INFO_XML} ${TEST}_homogenize.ref.info.xml ${TEST}_homogenize.info.xml --eps ${EPSILON} ${LAST}")
    MESSAGE(WARNING "COMPARE_INFO_XML_RETVAL: ${COMPARE_INFO_XML_RETVAL}")
    MESSAGE(WARNING "COMPARE_INFO_XML_ERROR: ${COMPARE_INFO_XML_ERROR}")
    MESSAGE(FATAL_ERROR "Homogenization results for test case '${TEST_NAME}' differ from reference.")
  ENDIF()
  
ENDMACRO(RUN_AND_TEST_BASECELL_HOMOGENIZATION)


MACRO(DUMP)
  get_cmake_property(_variableNames VARIABLES)
  foreach (_variableName ${_variableNames})
    message("${_variableName}=${${_variableName}}")
  endforeach()
ENDMACRO()

