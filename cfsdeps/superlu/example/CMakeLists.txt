#######################################################################
#  CMakeLists.txt file for creating the example programs for the 
#  linear equation routines in SuperLU.  The files are grouped as follows:
#
#       SLINEXM -- Single precision real example routines
#       DLINEXM -- Double precision real example routines
#       CLINEXM -- Double precision complex example routines
#       ZLINEXM -- Double precision complex example routines
#
#######################################################################

include_directories(../SRC)

set(SLINEXM	slinsol.c)
set(SLINEXM1	slinsol1.c)
set(SLINXEXM 	slinsolx.c)
set(SLINXEXM1 	slinsolx1.c)
set(SLINXEXM2	slinsolx2.c)
set(SITSOL     	sitersol.c sfgmr.c)
# set(SITSOL1    	sitersol1.c sfgmr.c)

set(DLINEXM	dlinsol.c)
set(DLINEXM1	dlinsol1.c)
set(DLINXEXM 	dlinsolx.c)
set(DLINXEXM1 	dlinsolx1.c)
set(DLINXEXM2 	dlinsolx2.c)
set(SUPERLUEXM 	superlu.c sp_ienv.c)
set(DITSOL	ditersol.c dfgmr.c)
# set(DITSOL1	ditersol1.c dfgmr.c)

set(CLINEXM 	clinsol.c)
set(CLINEXM1 	clinsol1.c)
set(CLINXEXM 	clinsolx.c)
set(CLINXEXM1 	clinsolx1.c)
set(CLINXEXM2 	clinsolx2.c)
set(CITSOL	citersol.c cfgmr.c)
# set(CITSOL1	citersol1.c cfgmr.c)

set(ZLINEXM 	zlinsol.c)
set(ZLINEXM1 	zlinsol1.c)
set(ZLINXEXM 	zlinsolx.c)
set(ZLINXEXM1 	zlinsolx1.c)
set(ZLINXEXM2 	zlinsolx2.c)
set(ZITSOL	zitersol.c zfgmr.c)
# set(ZITSOL1	zitersol1.c zfgmr.c)

SET(TARGET_LL 
  superlu
  ${BLAS_LIBRARIES}
  ${M_LIBRARY}
  )

IF(WIN32)
  IF(MSVC OR CMAKE_C_COMPILER_ID STREQUAL "Intel")
    ADD_LIBRARY(xgetopt STATIC xgetopt.c)
    LIST(APPEND TARGET_LL xgetopt)
  ENDIF()
ENDIF()


add_executable(slinsol ${SLINEXM})
target_link_libraries(slinsol ${TARGET_LL})

add_executable(slinsol1 ${SLINEXM1})
target_link_libraries(slinsol1 ${TARGET_LL}) 

add_executable(slinsolx ${SLINXEXM})
target_link_libraries(slinsolx ${TARGET_LL}) 

add_executable(slinsolx1 ${SLINXEXM1})
target_link_libraries(slinsolx1 ${TARGET_LL}) 

add_executable(slinsolx2 ${SLINXEXM2})
target_link_libraries(slinsolx2 ${TARGET_LL}) 

add_executable(sitersol ${SITSOL})
target_link_libraries(sitersol ${TARGET_LL}) 

#add_executable(sitersol1 ${SITSOL1})
#target_link_libraries(sitersol1 ${TARGET_LL}) 

add_executable(dlinsol ${DLINEXM})
target_link_libraries(dlinsol ${TARGET_LL}) 

add_executable(dlinsol1 ${DLINEXM1})
target_link_libraries(dlinsol1 ${TARGET_LL}) 

add_executable(dlinsolx ${DLINXEXM})
target_link_libraries(dlinsolx ${TARGET_LL}) 

add_executable(dlinsolx1 ${DLINXEXM1})
target_link_libraries(dlinsolx1 ${TARGET_LL}) 

add_executable(dlinsolx2 ${DLINXEXM2})
target_link_libraries(dlinsolx2 ${TARGET_LL}) 

add_executable(superlu_example ${SUPERLUEXM})
target_link_libraries(superlu_example ${TARGET_LL}) 

add_executable(ditersol ${DITSOL})
target_link_libraries(ditersol ${TARGET_LL}) 

#add_executable(ditersol1 ${DITSOL1})
#target_link_libraries(ditersol1 ${TARGET_LL}) 

add_executable(clinsol ${CLINEXM})
target_link_libraries(clinsol ${TARGET_LL}) 

add_executable(clinsol1 ${CLINEXM1})
target_link_libraries(clinsol1 ${TARGET_LL}) 

add_executable(clinsolx ${CLINXEXM})
target_link_libraries(clinsolx ${TARGET_LL}) 

add_executable(clinsolx1 ${CLINXEXM1})
target_link_libraries(clinsolx1 ${TARGET_LL}) 

add_executable(clinsolx2 ${CLINXEXM2})
target_link_libraries(clinsolx2 ${TARGET_LL}) 

add_executable(citersol ${CITSOL})
target_link_libraries(citersol ${TARGET_LL}) 

#add_executable(citersol1 ${CITSOL1})
#target_link_libraries(citersol1 ${TARGET_LL}) 

add_executable(zlinsol ${ZLINEXM})
target_link_libraries(zlinsol ${TARGET_LL}) 

add_executable(zlinsol1 ${ZLINEXM1})
target_link_libraries(zlinsol1 ${TARGET_LL}) 

add_executable(zlinsolx ${ZLINXEXM})
target_link_libraries(zlinsolx ${TARGET_LL}) 

add_executable(zlinsolx1 ${ZLINXEXM1})
target_link_libraries(zlinsolx1 ${TARGET_LL}) 

add_executable(zlinsolx2 ${ZLINXEXM2})
target_link_libraries(zlinsolx2 ${TARGET_LL}) 

add_executable(zitersol ${ZITSOL})
target_link_libraries(zitersol ${TARGET_LL}) 

#add_executable(zitersol1 ${ZITSOL1})
#target_link_libraries(zitersol1 ${TARGET_LL}) 
