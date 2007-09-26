// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <def_use_lapack.hh>

#include "matvec/matvec.hh"
#include "precond/precond.hh"
#include "precond/idprecond.hh"
#include "precond/generateprecond.hh"

// include source code for templated preconditioners
#include "precond/jacprecond.cc"
#include "precond/ssorprecond.cc"
#include "precond/ilu0precond.cc"
#include "precond/ilutpprecond.cc"
#include "precond/ilukprecond.cc"
//#include "precond/mgprecond.cc"
#include "precond/ILDLPrecond/ildlprecond.cc"
#include "precond/ic0precond.cc"

#if defined( __INTEL_COMPILER)
#ifdef USE_LAPACK
#include "external/lapack/lapack.hh"
#include "external/lapack/lapackgbmatrix.cc"
#endif
#endif

#include "matvec/diag_matrix.hh"

namespace OLAS {


  // *********************************************************
  //   Central macro for generation of preconditioner object
  //   Its use also takes care of the instantiation in the
  //   case of the templated preconditioners
  // *********************************************************
#define PRECOND_OBJ(matEntry,matStore,numDof,precondObjType)\
if ((entryType==matEntry) && (storageType==matStore) && (dofNum==numDof)) {\
retVal = new precondObjType( mat, myParams, myReport );\
(*cla) << " GenerateStdPrecondObject: Generated "\
       << MACRO2STRING(precondObjType)\
       << " preconditioner" << std::endl; }


  // **************************
  //   Jacobi Preconditioners
  // **************************

  // Compressed Row Storage Format
  typedef JacPrecond<CRS_Matrix<Double>,Double> RCRSJacDof1;
  typedef JacPrecond<CRS_Matrix<Complex>,Complex> CCRSJacDof1;

  // Symmetric Compressed Row Storage Format
  typedef JacPrecond<SCRS_Matrix<Double>,Double> RSCRSJacDof1;
  typedef JacPrecond<SCRS_Matrix<Complex>,Complex> CSCRSJacDof1;

  // Diagonal matrix
  typedef JacPrecond<Diag_Matrix<Double>,Double> RDIAGJacDof1;
  typedef JacPrecond<Diag_Matrix<Complex>,Complex> CDIAGJacDof1;


  // **********************
  //   Generation routine
  // **********************
  BaseStdPrecond* GenerateStdPrecondObject( const StdMatrix &mat, 
                                            PrecondType ptype,
					    OLAS_Params *myParams,
					    OLAS_Report *myReport ) {


    BaseStdPrecond *retVal = NULL;

    // Obtain matrix type information
    MatrixEntryType entryType = mat.GetEntryType();
    MatrixStorageType storageType = mat.GetStorageType();
    Integer dofNum = mat.GetBlockSize();

    // Branch depending on desired preconditioner
    switch( ptype ) {

      // ==============================
      //   Not really preconditioners
      // ==============================

      // Note: If no preconditioner was demanded we will currently also
      //       generate an identity preconditioner to remain consistent.
    case NOPRECOND:
    case ID:
      retVal = new IdPrecond;
      (*cla) << " GenerateStdPrecondObject: Generated Identity preconditioner"
	     << std::endl;
      break;

      // ========================= 
      //   Jacobi Preconditioner
      // ========================= 
    case JACOBI:

      // real crs jacobi
      PRECOND_OBJ(DOUBLE,SPARSE_NONSYM,1,RCRSJacDof1);
      //complex crs jacobi
      PRECOND_OBJ(COMPLEX,SPARSE_NONSYM,1,CCRSJacDof1);

      //real scrs jacobi
      PRECOND_OBJ(DOUBLE,SPARSE_SYM,1,RSCRSJacDof1);
      //complex scrs jacobi
      PRECOND_OBJ(COMPLEX,DIAG,1,CSCRSJacDof1);

      //real diag jacobi
      PRECOND_OBJ(DOUBLE,DIAG,1,RDIAGJacDof1);
      //complex diag jacobi
      PRECOND_OBJ(COMPLEX,SPARSE_SYM,1,CDIAGJacDof1);

      break;

      // ======================= 
      //   SSOR Preconditioner
      // ======================= 
    case SSOR:

      // real crs SSOR
      PRECOND_OBJ(DOUBLE,SPARSE_NONSYM,1,SSORPrecond<Double>);
      //complex crs SSOR
      PRECOND_OBJ(COMPLEX,SPARSE_NONSYM,1,SSORPrecond<Complex>);

      break;

      // ======================= 
      //   ILU0 Preconditioner
      // ======================= 
    case ILU0:

      // real scalar ILU(0) (currently only in CRS)
      PRECOND_OBJ( DOUBLE, SPARSE_NONSYM, 1, ILU0Precond<Double> );

      // complex scalar ILU(0) (currently only in CRS)
      PRECOND_OBJ( COMPLEX, SPARSE_NONSYM, 1, ILU0Precond<Complex> );

      // Test, if preconditioner object could be generated. If not, then
      // this is most likely due to a request for an unsupported version
      // of ILUPrecond.
      if ( retVal == NULL ) { 
	(*error) << " Request for unsupported version of ILU0Precond";
	Error( __FILE__, __LINE__ );
      }
      break;

      // ============================= 
      //   ILU(tau,p) Preconditioner
      // ============================= 
    case ILUTP:

      // real scalar ILU(tau,p) (works with CRS_Matrix only)
      PRECOND_OBJ( DOUBLE, SPARSE_NONSYM, 1, ILUTP_Precond<Double> );

      // complex scalar ILU(tau,p) (works with CRS_Matrix only)
      PRECOND_OBJ( COMPLEX, SPARSE_NONSYM, 1, ILUTP_Precond<Complex> );

      // Test, if preconditioner object could be generated. If not, then
      // this is most likely due to a request for an unsupported version
      // of ILUTP_Precond.
      if ( retVal == NULL ) { 
	(*error) << " Generation of ILUTP_Precond object failed! "
		 << "Most likely this is due to request for an unsupported "
		 << "version, i.e. non CRS_Matrix or non-scalar entries!";
	Error( __FILE__, __LINE__ );
      }
      break;

      // ========================= 
      //   ILU(k) Preconditioner
      // ========================= 
    case ILUK:

      // real scalar ILU(k) (works with CRS_Matrix only)
      PRECOND_OBJ( DOUBLE, SPARSE_NONSYM, 1, ILUK_Precond<Double> );

      // complex scalar ILU(k) (works with CRS_Matrix only)
      PRECOND_OBJ( COMPLEX, SPARSE_NONSYM, 1, ILUK_Precond<Complex> );

      // Test, if preconditioner object could be generated. If not, then
      // this is most likely due to a request for an unsupported version
      // of ILUK_Precond.
      if ( retVal == NULL ) { 
	(*error) << " Generation of ILUK_Precond object failed! "
		 << "Most likely this is due to request for an unsupported "
		 << "version, i.e. non CRS_Matrix or non-scalar entries!";
	Error( __FILE__, __LINE__ );
      }
      break;

      // ========================
      //   ILDL Preconditioners
      // ========================

      // The following are all sub-variants of the same preconditioner
      // class. We therefore generate for all variants an object of type
      // ILDLPrecond. The latter generates a fitting factoriser internally
      // to reflect the specific variant
    case ILDL0:
    case ILDLK:
    case ILDLTP:
    case ILDLCN:

      // Set corresponding subtype
      myParams->SetValue( "ILDLPRECOND_subType", ptype );

      // Generate preconditioner
      PRECOND_OBJ( DOUBLE,  SPARSE_SYM, 1, ILDLPrecond<Double> );
      PRECOND_OBJ( COMPLEX, SPARSE_SYM, 1, ILDLPrecond<Complex> );

      // Test, if preconditioner object could be generated. If not, then
      // this is most likely due to a request for an unsupported version.
      if ( retVal == NULL ) { 
	(*error) << "Generation of ILDLPrecond object (variant = "
                 << Enum2String(ptype) << ") failed! "
		 << "Most likely this is due to request for an unsupported "
		 << "version, i.e. non SCRS_Matrix or non-scalar entries!";
	Error( __FILE__, __LINE__ );
      }
      break;

      // ======================= 
      //   IC0 Preconditioner
      // ======================= 
    case IC0:

      // real scalar IC(0) (currently only in CRS)
      PRECOND_OBJ( DOUBLE, SPARSE_SYM, 1, IC0Precond<Double> );

      // complex scalar IC(0) (currently only in CRS)
      PRECOND_OBJ( COMPLEX, SPARSE_SYM, 1, IC0Precond<Complex> );

      // Test, if preconditioner object could be generated. If not, then
      // this is most likely due to a request for an unsupported version
      // of ICPrecond.
      if ( retVal == NULL ) { 
	(*error) << " Request for unsupported version of IC0Precond";
	Error( __FILE__, __LINE__ );
      }
      break;

      // ============================
      //   Multigrid Preconditioner
      // ============================
    case MG:

      // multigrid preconditioners can be used with CRS_Matrix only.
      // PRECOND_OBJ( DOUBLE, SPARSE_NONSYM, 1, MGPrecond<Double> );

      // Test, if preconditioner object could be generated. If not, then
      // this is most likely due to a request for an unsupported version
      // of MGPrecond.
      if ( retVal == NULL ) { 
          (*error) << " Request for unsupported version of MGPrecond";
          Error( __FILE__, __LINE__ );
      }
      break;

      // ======================
      //   Something's broken
      // ======================

    default:
      Error( "GeneratePrecondObject failed: Preconditioner type unknown",
             __FILE__, __LINE__ );
    }

    // test, if preconditioner object could be generated
    if ( retVal == NULL ) {
      (*error) << "Something's broke. No preconditioner was generated!";
      Error( __FILE__, __LINE__ );
    }

    // Finished
    return retVal;
  }
}
