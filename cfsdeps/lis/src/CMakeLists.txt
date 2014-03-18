# src/
# Get source files
FILE(GLOB temp_Files esolver/*.h esolver/*.c)
SOURCE_GROUP( lis\\src\\esolver FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files}) 
IF (LIS_ENABLE_SAAMG OR LIS_ENABLE_FORTRAN)
	FILE(GLOB temp_Files fortran/*.h fortran/*.c fortran/*.f fortran/*.F)
	SOURCE_GROUP( lis\\src\\fortran_c FILES ${temp_Files} )
	LIST(APPEND lis_src_Files ${temp_Files}) 
ENDIF()
IF (LIS_ENABLE_SAAMG)
	IF (LIS_ENABLE_MPI)
		SET(temp_Files 
			fortran/lis_m_data_type_AMGCG.F90
			fortran/amg/lis_m_sort.F90
			fortran/amg/lis_m_queue.F90
			fortran/amg/lis_m_hash.F90
			fortran/amg/lis_m_count_time_mod.F90
			fortran/amg/lis_m_solver_SR2.F90
			fortran/amg/lis_m_solver_Gnum.F90
			fortran/amg/lis_m_aggregate.F90
			fortran/amg/lis_m_data_creation_AMGCG.F90
			fortran/amg/lis_m_solver_AMGCG.F90
			fortran/amg/lis_m_finit.F90
		)
	ELSE()
		SET(temp_Files 
			fortran/amg/lis_s_data_type_AMGCG.F90
			fortran/amg/lis_s_sort.F90
			fortran/amg/lis_s_queue.F90
			fortran/amg/lis_s_aggregate.F90
			fortran/amg/lis_s_data_creation_AMGCG.F90
			fortran/amg/lis_s_solver_AMGCG.F90
			fortran/amg/lis_s_finit.F90
		)
	ENDIF()	
	SOURCE_GROUP( lis\\src\\fortran\\amg FILES ${temp_Files} )
	LIST(APPEND lis_src_Files ${temp_Files}) 
ENDIF()
FILE(GLOB temp_Files matrix/*.h matrix/*.c)
SOURCE_GROUP( lis\\src\\matrix FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files}) 
FILE(GLOB temp_Files matvec/*.h matvec/*.c)
SOURCE_GROUP( lis\\src\\matvec FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files}) 
FILE(GLOB temp_Files precision/*.h precision/*.c)
SOURCE_GROUP( lis\\src\\precision FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files}) 
FILE(GLOB temp_Files precon/*.h precon/*.c)
SOURCE_GROUP( lis\\src\\precon FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files}) 
FILE(GLOB temp_Files solver/*.h solver/*.c)
SOURCE_GROUP( lis\\src\\solver FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files}) 
FILE(GLOB temp_Files system/*.h system/*.c)
SOURCE_GROUP( lis\\src\\system FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files}) 
FILE(GLOB temp_Files vector/*.h vector/*.c)
SOURCE_GROUP( lis\\src\\vector FILES ${temp_Files} )
LIST(APPEND lis_src_Files ${temp_Files}) 

FILE(GLOB lis_include_Files ../include/*.h)
SOURCE_GROUP( lis\\include FILES ${lis_include_Files} )


# Create the executable
INCLUDE_DIRECTORIES(
	.
	../include
)

ADD_LIBRARY( lis
        ${lis_src_Files}
)
#SET_TARGET_PROPERTIES(lis PROPERTIES LINKER_LANGUAGE C)

# Link the standard C math library
# find_library(M_LIB m)
#TARGET_LINK_LIBRARIES( lis ${M_LIB} )

# Link MPI
IF (LIS_ENABLE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
	TARGET_LINK_LIBRARIES( lis ${MPI_C_LIBRARIES} )
ENDIF ()

INSTALL(TARGETS lis
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}"
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}"
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}")

SET(LIS_INCLUDE_FILES
  ../include/lis.h
  ../include/lisf.h
  ${PROJECT_BINARY_DIR}/lis_config.h
  ../include/lis_precon.h
)

install(FILES ${LIS_INCLUDE_FILES} DESTINATION ${INSTALL_INC_DIR})

