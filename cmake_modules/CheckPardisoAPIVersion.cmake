IF(CFS_PARDISO STREQUAL MKL)
  SET(PARDISO_API_VER 3)

  IF(UNIX)
    IF(DEFINED MKL_UPDATE)
      MESSAGE(STATUS "Intel MKL version ${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}.${MKL_UPDATE} detected. Setting Pardiso API version to ${PARDISO_API_VER}.")
    ELSE()
      MESSAGE(STATUS "Intel MKL version ${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION} detected. Setting Pardiso API version to ${PARDISO_API_VER}.")
    ENDIF()
  ELSE()
    MESSAGE(STATUS "Using Intel MKL. Setting Pardiso API version to ${PARDISO_API_VER}.")
  ENDIF()
ELSE()
#-----------------------------------------------------------------------------
# Set source code for testing API version of PARDISO.
#-----------------------------------------------------------------------------
SET(PARDISO_API_CHECK_SOURCE 
"#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>

/* PARDISO prototypes according to API version. */
 
#if PARDISO_API_VER == 3
extern \"C\" {
  void pardisoinit_ (void *pt[64], int *mtype, int iparm[64]);

  void pardiso_ (void *pt[64], int *maxfct, int *mnum, int *mtype,
                 int *phase, int *n, double a[], int ia[], int ja[],
                 int perm[], int *nrhs, int iparm[64], int *msglvl,
                 double b[], double x[], int *error);
}
#endif

#if PARDISO_API_VER == 4
extern \"C\" {
  void pardisoinit_ (void *pt[64], int *mtype, int *solver,
                     int iparm[64], double dparm[64], int *error);

  void pardiso_ (void *pt[64], int *maxfct, int *mnum, int *mtype,
                 int *phase, int *n, double a[], int ia[], int ja[],
                 int perm[], int *nrhs, int iparm[64], int *msglvl,
                 double b[], double x[], int *error, double dparm[64] );
}
#endif

/* -------------------------------------------------------------------- */
/*      Example program to show the use of the \"PARDISO\" routine        */
/*      on symmetric linear systems                                     */
/* -------------------------------------------------------------------- */
/*      This program can be downloaded from the following site:         */
/*      http://www.pardiso-project.org                                  */
/*                                                                      */
/*  (C) Olaf Schenk, Department of Computer Science,                    */
/*      University of Basel, Switzerland.                               */
/*      Email: olaf.schenk@unibas.ch                                    */
/* -------------------------------------------------------------------- */

int main (int argc, char* argv[])
{
    /* Matrix data. */
    int    n = 8;
    int    ia[ 9] = { 0, 4, 7, 9, 11, 14, 16, 17, 18 };
    int    ja[18] = { 0,    2,       5, 6,
                         1, 2,    4,
                            2,             7,
                               3,       6,
                                  4, 5, 6,
                                     5,    7,
                                        6,
                                           7 };
    double  a[18] = { 7.0,      1.0,           2.0, 7.0,
                          -4.0, 8.0,           2.0,
                                1.0,                     5.0,
                                     7.0,           9.0,
                                          5.0, 1.0, 5.0,
                                               0.0,      5.0,
                                                   11.0,
                                                         5.0 };

    int      nnz = ia[n];
    int      mtype = -2;        /* Real symmetric matrix */

    /* RHS and solution vectors. */
    double   b[8], x[8];
    std::fill(b, b+8, 0);
    std::fill(x, x+8, 0);
    int      nrhs = 1;          /* Number of right hand sides. */

    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */ 
    void    *pt[64];
    for(int i=0; i<64; i++) {
      pt[i] = NULL;
    }

    // We interpret the parameters to pardisoinit_ as consecutive data.
    // This way we can make sure, that the solver parameter will be overwritten
    // by default parameters by pardisoinit_ in the case of a Pardiso 3.x
    // library being called through the Pardiso 4.0 interface. We can afterwards
    // check for this condition.
    
    int param_array[66];
    
    /* Pardiso control parameters. */
    int* solver = param_array;
    int* iparm = param_array+1;
    
    std::fill(param_array, param_array+66, -666);
    solver[0] = 0; // Use direct solver    
    iparm[0] = 0; // IPARM(1) sets iparm array to default values if 0

    double dparm[64];
    std::fill(dparm, dparm+64, -999);
    
    int      maxfct, mnum, phase, error, msglvl;

    /* Number of processors. */
    int      num_procs;

    /* Auxiliary variables. */
    int      i;

    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */

   
	/* -------------------------------------------------------------------- */
	/* ..  Setup Pardiso control parameters.                                */
	/* -------------------------------------------------------------------- */
	error = 0;
#if PARDISO_API_VER == 3
    pardisoinit_ (pt,  &mtype, iparm); 
//    std::string fn = \"/tmp/pardiso_log_3.txt\";
#endif
#if PARDISO_API_VER == 4
    pardisoinit_ (pt,  &mtype, solver, iparm, dparm, &error); 
//    std::string fn = \"/tmp/pardiso_log_4.txt\";
#endif
   

//    std::ofstream os(fn.c_str());
//    os << \"solver : \" << (*solver) << std::endl;
//    for(int i=0; i<66; i++) {
//       os << i << \": \" << iparm[i] << \"  \" << dparm[i] << std::endl;
//    }
//    os.close();
   
    if (error != 0) 
    {
        if (error == -10 )
           std::cout << \"No license file found\" << std::endl;
        if (error == -11 )
           std::cout << \"License is expired\" << std::endl;
        if (error == -12 )
           std::cout << \"Wrong username or hostname\" << std::endl;
	   
        std::cout << std::endl << \"At least we have a defined error value and therefore know that the API is working.\" << std::endl;	   
        return -error; 
    }
    else
        std::cout << \"PARDISO license check was successful...\" << std::endl;

    // Check if everything is where it is supposed to be
    if(*solver != 0  || // Stands for the solver, must always be 0 (direct)
       iparm[0] != 1  || // IPARM(1) must always be 1 after call to pardisoinit_
       iparm[63] != 0)   // If calling a 3.x lib with 4.x API iparm[64] == -666
    {
      std::cout << \"IPARM array is not the way it should be...\" << std::endl;
      return 66;
    }

    num_procs = 1;
    iparm[2]  = num_procs;


    maxfct = 1;		/* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */
    
    msglvl = 1;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */

	/* -------------------------------------------------------------------- */
	/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
	/*     notation.                                                        */
	/* -------------------------------------------------------------------- */
    for (i = 0; i < n+1; i++) {
        ia[i] += 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] += 1;
    }
 

	/* -------------------------------------------------------------------- */
	/* ..  Reordering and Symbolic Factorization.  This step also allocates */
	/*     all memory that is necessary for the factorization.              */
	/* -------------------------------------------------------------------- */
    phase = 11; 
    error = 0;
    
#if PARDISO_API_VER == 3
    pardiso_ (pt, &maxfct, &mnum, &mtype, &phase,
		       &n, a, ia, ja, &idum, &nrhs,
		       iparm, &msglvl, &ddum, &ddum, &error);
#endif
#if PARDISO_API_VER == 4
    pardiso_ (pt, &maxfct, &mnum, &mtype, &phase,
		       &n, a, ia, ja, &idum, &nrhs,
		       iparm, &msglvl, &ddum, &ddum, &error, dparm);
#endif
    
    if (error != 0) {
        std::cerr << std::endl << \"ERROR during symbolic factorization: \" << error << std::endl;
        switch(error) {
          case -10:
          case -11:
          case -12:
            std::cout << \"There exists a licensing problem which needs to be solved.\" << std::endl;
            return -error;
	    
          default:
            return phase;
        }
    }
    std::cout << std::endl << \"Reordering completed ... \" << std::endl;
    std::cout << std::endl << \"Number of nonzeros in factors  = \" << iparm[17] << std::endl;
    std::cout << std::endl << \"Number of factorization MFLOPS = \" << iparm[18] << std::endl;

	/* -------------------------------------------------------------------- */
	/* ..  Numerical factorization.                                         */
	/* -------------------------------------------------------------------- */    
    phase = 22;

    error = 0;
#if PARDISO_API_VER == 3
    pardiso_ (pt, &maxfct, &mnum, &mtype, &phase,
		       &n, a, ia, ja, &idum, &nrhs,
		       iparm, &msglvl, &ddum, &ddum, &error);
#endif
#if PARDISO_API_VER == 4
    pardiso_ (pt, &maxfct, &mnum, &mtype, &phase,
		       &n, a, ia, ja, &idum, &nrhs,
		       iparm, &msglvl, &ddum, &ddum, &error, dparm);
#endif
   
    if (error != 0) {
        std::cerr << std::endl << \"ERROR during numerical factorization: \" << error << std::endl;
        return phase; 
    }
    std::cout << std::endl << \"Factorization completed ...\" << std::endl;

	/* -------------------------------------------------------------------- */    
	/* ..  Back substitution and iterative refinement.                      */
	/* -------------------------------------------------------------------- */    
    phase = 33;

    iparm[7] = 1;       /* Max numbers of iterative refinement steps. */

    /* Set right hand side to one. */
    for (i = 0; i < n; i++) {
        b[i] = 1;
    }
   
    error = 0;
#if PARDISO_API_VER == 3
    pardiso_ (pt, &maxfct, &mnum, &mtype, &phase,
		       &n, a, ia, ja, &idum, &nrhs,
		       iparm, &msglvl, b, x, &error);
#endif
#if PARDISO_API_VER == 4
    pardiso_ (pt, &maxfct, &mnum, &mtype, &phase,
		       &n, a, ia, ja, &idum, &nrhs,
		       iparm, &msglvl, b, x, &error, dparm);
#endif
   
    if (error != 0) {
        std::cerr << std::endl << \"ERROR during solution: \" << error << std::endl;
        return phase;
    }

    std::cout << std::endl << \"Solve completed ... \" << std::endl;
    std::cout << std::endl << \"The solution of the system is: \";
    for (i = 0; i < n; i++) {
        std::cout << std::endl << \" x[\" << i << \"] = \" << x[i];
    }
    std::cout << std::endl;

    
	/* -------------------------------------------------------------------- */    
	/* ..  Convert matrix back to 0-based C-notation.                       */
	/* -------------------------------------------------------------------- */ 
    for (i = 0; i < n+1; i++) {
        ia[i] -= 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] -= 1;
    }

	/* -------------------------------------------------------------------- */    
	/* ..  Termination and release of memory.                               */
	/* -------------------------------------------------------------------- */    
    phase = -1;                 /* Release internal memory. */
    
#if PARDISO_API_VER == 3
    pardiso_ (pt, &maxfct, &mnum, &mtype, &phase,
		       &n, &ddum, ia, ja, &idum, &nrhs,
		       iparm, &msglvl, &ddum, &ddum, &error);
#endif
#if PARDISO_API_VER == 4
    pardiso_ (pt, &maxfct, &mnum, &mtype, &phase,
		       &n, &ddum, ia, ja, &idum, &nrhs,
		       iparm, &msglvl, &ddum, &ddum, &error, dparm);
#endif

    return 0;
}"
)

#-----------------------------------------------------------------------------
# Set required paths and libraries in order for the tests to find them.
#-----------------------------------------------------------------------------
# CMAKE_REQUIRED_FLAGS = string of compile command line flags
# CMAKE_REQUIRED_DEFINITIONS = list of macros to define (-DFOO=bar)
# CMAKE_REQUIRED_INCLUDES = list of include directories
# CMAKE_REQUIRED_LIBRARIES = list of libraries to link
# CFS_REQUIRED_LINK_DIRECTORIES = list of LINK_DIRECTORIES

IF(CFS_FORTRAN_LIBS)
  string(REPLACE ":"
         ";" CFS_FORTRAN_LIBS
         ${CFS_FORTRAN_LIBS})
ENDIF(CFS_FORTRAN_LIBS)

SET(CMAKE_REQUIRED_LIBRARIES 
  ${PARDISO_LIBRARY}
  ${LAPACK_LIBRARY}
  ${BLAS_LIBRARY}
  ${CFS_FORTRAN_LIBS}
  ${CFS_FORTRAN_DYNRT_LIBS})
  
#-----------------------------------------------------------------------------
# Check for the API version each time the user changes CFS_PARDISO.
#-----------------------------------------------------------------------------
IF(NOT PARDISO_API_VER_3_LAST_CFS_PARDISO STREQUAL CFS_PARDISO)
  #---------------------------------------------------------------------------
  # Reset internal cache variable in order for CFS_CHECK_CXX_SOURCE_RUNS to
  # reperform testing.
  #---------------------------------------------------------------------------
  SET(PARDISO_API_VER_3 "" CACHE INTERNAL "PARDISO_API_VER_3" FORCE)

  #---------------------------------------------------------------------------
  # Set a define to determine that we want to check for API version 3.
  #---------------------------------------------------------------------------
  SET(CMAKE_REQUIRED_DEFINITIONS "-DPARDISO_API_VER=3")

  #---------------------------------------------------------------------------
  # Perform test
  #---------------------------------------------------------------------------
  CFS_CHECK_CXX_SOURCE_RUNS("${PARDISO_API_CHECK_SOURCE}"
                            PARDISO_API_VER_3
                            "${IFORT_LIB_PATH}")


  #---------------------------------------------------------------------------
  # Only set last openCFS PARDISO variable if test was successful
  #---------------------------------------------------------------------------
  IF(PARDISO_API_VER_3_EXITCODE EQUAL 0)
    SET(PARDISO_API_VER_3_LAST_CFS_PARDISO ${CFS_PARDISO} CACHE INTERNAL "${CFS_PARDISO}" FORCE)
  ELSE()
	  IF(USE_BLAS_LAPACK STREQUAL "MKL" AND CMAKE_CXX_COMPILER_ID MATCHES "Intel")
       MESSAGE("Please try using the library MKL provided by the Intel compiler!")
       MESSAGE("Mixing runtime libraries provided by standalone MKL and Intel")
       MESSAGE("compiler may cause problems.")
    ENDIF()
  ENDIF()

  #---------------------------------------------------------------------------
  # Warn user about missing or corrupt license files.
  #---------------------------------------------------------------------------
  IF(PARDISO_API_VER_3_EXITCODE EQUAL 10)
    MESSAGE("Could not find PARDISO license.lic!\n"
            "Please get the file at http://www.pardiso-project.org"
            "and set the environment variable PARDISO_LIC_PATH to\n"
            "the directory containing pardiso.lic.")
  ENDIF(PARDISO_API_VER_3_EXITCODE EQUAL 10)
  IF(PARDISO_API_VER_3_EXITCODE EQUAL 11)
    MESSAGE("PARDISO license is expired.")
  ENDIF(PARDISO_API_VER_3_EXITCODE EQUAL 11)
  IF(PARDISO_API_VER_3_EXITCODE EQUAL 12)
    MESSAGE("Wrong username or hostname for PARDISO.")
  ENDIF(PARDISO_API_VER_3_EXITCODE EQUAL 12)

ENDIF(NOT PARDISO_API_VER_3_LAST_CFS_PARDISO STREQUAL CFS_PARDISO)


#-----------------------------------------------------------------------------
# Check for the API version each time the user changes CFS_PARDISO.
#-----------------------------------------------------------------------------
IF(NOT PARDISO_API_VER_4_LAST_CFS_PARDISO STREQUAL CFS_PARDISO)
  #---------------------------------------------------------------------------
  # Reset internal cache variable in order for CFS_CHECK_CXX_SOURCE_RUNS to
  # reperform testing.
  #---------------------------------------------------------------------------
  SET(PARDISO_API_VER_4 "" CACHE INTERNAL "PARDISO_API_VER_4" FORCE)

  #---------------------------------------------------------------------------
  # Set a define to determine that we want to check for API version 4.
  #---------------------------------------------------------------------------
  SET(CMAKE_REQUIRED_DEFINITIONS "-DPARDISO_API_VER=4")

  #---------------------------------------------------------------------------
  # Perform test
  #---------------------------------------------------------------------------
  CFS_CHECK_CXX_SOURCE_RUNS("${PARDISO_API_CHECK_SOURCE}"
                            PARDISO_API_VER_4
                            "${IFORT_LIB_PATH}")

  #---------------------------------------------------------------------------
  # Only set last openCFS PARDISO variable if test was successful
  #---------------------------------------------------------------------------
  IF(PARDISO_API_VER_4_EXITCODE EQUAL 0)
    SET(PARDISO_API_VER_4_LAST_CFS_PARDISO ${CFS_PARDISO}
        CACHE INTERNAL "${CFS_PARDISO}" FORCE)
  ENDIF(PARDISO_API_VER_4_EXITCODE EQUAL 0)
  
  #---------------------------------------------------------------------------
  # Warn user about missing or corrupt license files.
  #---------------------------------------------------------------------------
  IF(PARDISO_API_VER_4_EXITCODE EQUAL 10)
    MESSAGE("Could not find PARDISO license.lic!\n"
            "Please get the file at http://www.pardiso-project.org"
            "and set the environment variable PARDISO_LIC_PATH to\n"
            "the directory containing pardiso.lic.")
  ENDIF(PARDISO_API_VER_4_EXITCODE EQUAL 10)
  IF(PARDISO_API_VER_4_EXITCODE EQUAL 11)
    MESSAGE("PARDISO license is expired.")
  ENDIF(PARDISO_API_VER_4_EXITCODE EQUAL 11)
  IF(PARDISO_API_VER_4_EXITCODE EQUAL 12)
    MESSAGE("Wrong username or hostname for PARDISO.")
  ENDIF(PARDISO_API_VER_4_EXITCODE EQUAL 12)
  
ENDIF(NOT PARDISO_API_VER_4_LAST_CFS_PARDISO STREQUAL CFS_PARDISO)


#-----------------------------------------------------------------------------
# If both checks failed issue an error
#-----------------------------------------------------------------------------
IF(NOT PARDISO_API_VER_3 AND
   NOT PARDISO_API_VER_4)
   
  SET(PARDISO_API_VER_3_LAST_CFS_PARDISO "" CACHE INTERNAL "PARDISO_API_VER_3_LAST_CFS_PARDISO" FORCE)
  SET(PARDISO_API_VER_4_LAST_CFS_PARDISO "" CACHE INTERNAL "PARDISO_API_VER_4_LAST_CFS_PARDISO" FORCE)

  MESSAGE(FATAL_ERROR "Could not determine Pardiso API version!\n\n"
                      "If you are using the Intel compilers you\n"
                      "should make sure that you source iccvars.sh\n"
                      "and ifortvars.sh before running CMake.\n\n"
                      "You should investigate the problem in CMakeFiles/CMakeError.log.\n\n"
                      "If you know the PARDISO API version you may also\n"
                      "set it by hand in CMakeCache.txt(e.g. PARDISO_API_VER_3:INTERNAL=1)")

ELSE(NOT PARDISO_API_VER_3 AND
     NOT PARDISO_API_VER_4)

  #---------------------------------------------------------------------------
  # If both checks succeded also issue an error
  #---------------------------------------------------------------------------
  IF(PARDISO_API_VER_3 AND
     PARDISO_API_VER_4)
    SET(PARDISO_API_VER_3_LAST_CFS_PARDISO "" CACHE INTERNAL
        "PARDISO_API_VER_3_LAST_CFS_PARDISO" FORCE)
    SET(PARDISO_API_VER_4_LAST_CFS_PARDISO "" CACHE INTERNAL
        "PARDISO_API_VER_4_LAST_CFS_PARDISO" FORCE)

    MESSAGE(FATAL_ERROR "Checks suggest that both Pardiso API versions are present!\n"
                        "Please try to re-run CMake in order to get a correct result.\n\n"
                        "The case you have discovered is nonsense but can happen due to the \"excellent\" (<- note the irony)\n"
                        "software engineering practices used to implement the PARDISO API.\n"
                        "There is no function to determine the version, and the number and ordering\n"
                        "of parameters changed from version 3 to 4. THANK YOU VERY MUCH PARDISO DEVELOPERS!\n\n"
                        "If you know the PARDISO API version you may also\n"
                        "set it by hand in CMakeCache.txt(e.g. PARDISO_API_VER_3:INTERNAL=1)")
  ELSE(PARDISO_API_VER_3 AND
       PARDISO_API_VER_4)
       
    #-------------------------------------------------------------------------
    # Otherwise set the corresponding API version. And make sure that the
    # the other API version does not get rechecked in every CMake run.
    #-------------------------------------------------------------------------
    IF(PARDISO_API_VER_3)
      SET(PARDISO_API_VER 3)
      SET(PARDISO_API_VER_4_LAST_CFS_PARDISO "${CFS_PARDISO}" CACHE INTERNAL "PARDISO_API_VER_4_LAST_CFS_PARDISO" FORCE)      
    ENDIF(PARDISO_API_VER_3)

    IF(PARDISO_API_VER_4)
      SET(PARDISO_API_VER 4)
      SET(PARDISO_API_VER_3_LAST_CFS_PARDISO "${CFS_PARDISO}" CACHE INTERNAL "PARDISO_API_VER_3_LAST_CFS_PARDISO" FORCE)      
    ENDIF(PARDISO_API_VER_4)
       
  ENDIF(PARDISO_API_VER_3 AND
       PARDISO_API_VER_4)
       
ENDIF(NOT PARDISO_API_VER_3 AND
      NOT PARDISO_API_VER_4)

ENDIF() # the not MKL version - practivally dead!
