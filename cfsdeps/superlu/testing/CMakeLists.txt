#######################################################################
#  CMakeLists.txt file for the test programs for the linear equation
#  routines in SuperLU. 
#######################################################################

macro(LOCAL_TARGET execname)
  if(NOT CMAKE_CONFIGURATION_TYPES)
    set_target_properties(${execname} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  else(NOT CMAKE_CONFIGURATION_TYPES)
    foreach(CFG_TYPE ${CMAKE_CONFIGURATION_TYPES})
      string(TOUPPER "${CFG_TYPE}" CFG_TYPE_UPPER)
      set_target_properties(${execname}	PROPERTIES RUNTIME_OUTPUT_DIRECTORY_${CFG_TYPE_UPPER} ${CMAKE_CURRENT_BINARY_DIR}/${CFG_TYPE})
    endforeach(CFG_TYPE ${CMAKE_CONFIGURATION_TYPES})
  endif(NOT CMAKE_CONFIGURATION_TYPES)
endmacro(LOCAL_TARGET execname)

macro(SUPERLU_TEST testname)
  add_custom_command(
    OUTPUT ${testname}.out
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/${testname}.cmake
    DEPENDS superlu_${testname}
    )
  add_custom_target(${testname}_test ALL DEPENDS ${testname}.out)
endmacro(SUPERLU_TEST testname)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../SRC)

add_subdirectory(MATGEN) 

SET(TARGET_LL tmglib superlu ${BLAS_LIBRARIES} ${M_LIBRARY})
IF(WIN32)
  IF(MSVC OR CMAKE_C_COMPILER_ID STREQUAL "Intel")
    LIST(APPEND TARGET_LL xgetopt)
  ENDIF()
ENDIF()


# SLINTST -- Single precision real test routines
set(stest_SRCS sp_ienv.c sdrive.c sp_sconvert.c sgst01.c sgst02.c sgst04.c sgst07.c)
add_executable(superlu_stest ${stest_SRCS})
target_link_libraries(superlu_stest ${TARGET_LL})
LOCAL_TARGET(superlu_stest)
set(SUPERLU_STEST_EXEC ${CMAKE_CURRENT_BINARY_DIR}/superlu_stest)
configure_file(stest.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/stest.cmake @ONLY)
SUPERLU_TEST(stest)

# DLINTST -- Double precision real test routines
set(dtest_SRCS sp_ienv.c ddrive.c sp_dconvert.c dgst01.c dgst02.c dgst04.c dgst07.c)
add_executable(superlu_dtest ${dtest_SRCS})
target_link_libraries(superlu_dtest ${TARGET_LL})
LOCAL_TARGET(superlu_dtest)
set(SUPERLU_DTEST_EXEC ${CMAKE_CURRENT_BINARY_DIR}/superlu_dtest)
configure_file(dtest.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/dtest.cmake @ONLY)
SUPERLU_TEST(dtest)

# CLINTST -- Double precision complex test routines
set(ctest_SRCS sp_ienv.c cdrive.c sp_cconvert.c cgst01.c cgst02.c cgst04.c cgst07.c)
add_executable(superlu_ctest ${ctest_SRCS})
target_link_libraries(superlu_ctest ${TARGET_LL})
LOCAL_TARGET(superlu_ctest)
set(SUPERLU_CTEST_EXEC ${CMAKE_CURRENT_BINARY_DIR}/superlu_ctest)
configure_file(ctest.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/ctest.cmake @ONLY)
SUPERLU_TEST(ctest)

# ZLINTST -- Double precision complex test routines
set(ztest_SRCS sp_ienv.c zdrive.c sp_zconvert.c zgst01.c zgst02.c zgst04.c zgst07.c)
add_executable(superlu_ztest ${ztest_SRCS})
target_link_libraries(superlu_ztest ${TARGET_LL})
LOCAL_TARGET(superlu_ztest)
set(SUPERLU_ZTEST_EXEC ${CMAKE_CURRENT_BINARY_DIR}/superlu_ztest)
configure_file(ztest.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/ztest.cmake @ONLY)
SUPERLU_TEST(ztest)


# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
