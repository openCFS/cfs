// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <def_use_lapack.hh>

#include "OLAS/precond/idprecond.hh"
#include "OLAS/precond/generateprecond.hh"

// include source code for templated preconditioners
#include "OLAS/precond/jacprecond.hh"
#include "OLAS/precond/ssorprecond.hh"
#include "OLAS/precond/ilu0precond.hh"
#include "OLAS/precond/ilutpprecond.hh"
#include "OLAS/precond/ilukprecond.hh"
//#include "precond/mgmakeprecond.cc"
#include "OLAS/precond/ILDLPrecond/ildlprecond.hh"
#include "OLAS/precond/ic0precond.hh"

#include "MatVec/diag_matrix.hh"

namespace CoupledField {


  // *********************************************************
  //   Central macro for generation of preconditioner object
  //   Its use also takes care of the instantiation in the
  //   case of the templated preconditioners
  // *********************************************************
#define PRECOND_OBJ(matEntry,matStore,precondObjType)\
if ((entryType==matEntry) && (storageType==matStore) ) {\
retVal = new precondObjType( mat, solverNode, olasInfo );\
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
                                            ParamNode* solverNode,
                                            InfoNode *olasInfo ) {


    BaseStdPrecond *retVal = NULL;

    // Obtain matrix type information
    BaseMatrix::EntryType entryType = mat.GetEntryType();
    BaseMatrix::StorageType storageType = mat.GetStorageType();
    std::string precondStr = "";

    ParamNode* solverXML = NULL;
    solverXML = solverNode->Get("solver", false);

    EnumMap::iterator it, end;
    it = BasePrecond::precondType.map.begin();
    end = BasePrecond::precondType.map.end();
    
    for( ; it != end; it++ ) {
      if( solverXML->Has("precond", it->second) ) {
        if(precondStr != "")
          EXCEPTION("Two preconditioners have been specified: " << precondStr << " and " << (it->second))
        
        precondStr = it->second;
      }
    }
  
    if(precondStr == "") {
      precondStr = "noPrecond";
    }

    BasePrecond::PrecondType ptype = BasePrecond::NOPRECOND;
    ptype = BasePrecond::precondType.Parse(precondStr);
    
    // Branch depending on desired preconditioner
    switch( ptype ) {

      // ==============================
      //   Not really preconditioners
      // ==============================

      // Note: If no preconditioner was demanded we will currently also
      //       generate an identity preconditioner to remain consistent.
    case BasePrecond::NOPRECOND:
    case BasePrecond::ID:
      retVal = new IdPrecondStd;
      (*cla) << " GenerateStdPrecondObject: Generated Identity preconditioner"
	     << std::endl;
      break;

      // =========================
      //   Jacobi Preconditioner
      // =========================
    case BasePrecond::JACOBI:

      // real crs jacobi
      PRECOND_OBJ(BaseMatrix::DOUBLE,BaseMatrix::SPARSE_NONSYM,RCRSJacDof1);
      //complex crs jacobi
      PRECOND_OBJ(BaseMatrix::COMPLEX,BaseMatrix::SPARSE_NONSYM,CCRSJacDof1);

      //real scrs jacobi
      PRECOND_OBJ(BaseMatrix::DOUBLE,BaseMatrix::SPARSE_SYM,RSCRSJacDof1);
      //complex scrs jacobi
      PRECOND_OBJ(BaseMatrix::COMPLEX,BaseMatrix::DIAG,CSCRSJacDof1);

      //real diag jacobi
      PRECOND_OBJ(BaseMatrix::DOUBLE,BaseMatrix::DIAG,RDIAGJacDof1);
      //complex diag jacobi
      PRECOND_OBJ(BaseMatrix::COMPLEX,BaseMatrix::SPARSE_SYM,CDIAGJacDof1);

      break;

      // =======================
      //   SSOR Preconditioner
      // =======================
    case BasePrecond::SSOR:

      // real crs SSOR
      PRECOND_OBJ(BaseMatrix::DOUBLE,BaseMatrix::SPARSE_NONSYM,SSORPrecond<Double>);
      //complex crs SSOR
      PRECOND_OBJ(BaseMatrix::COMPLEX,BaseMatrix::SPARSE_NONSYM,SSORPrecond<Complex>);

      break;

      // =======================
      //   ILU0 Preconditioner
      // =======================
    case BasePrecond::ILU0:

      // real scalar ILU(0) (currently only in CRS)
      PRECOND_OBJ( BaseMatrix::DOUBLE, BaseMatrix::SPARSE_NONSYM,  ILU0Precond<Double> );

      // complex scalar ILU(0) (currently only in CRS)
      PRECOND_OBJ( BaseMatrix::COMPLEX, BaseMatrix::SPARSE_NONSYM, ILU0Precond<Complex> );

      // Test, if preconditioner object could be generated. If not, then
      // this is most likely due to a request for an unsupported version
      // of ILUPrecond.
      if ( retVal == NULL ) {
        EXCEPTION( " Request for unsupported version of ILU0Precond" );
      }
      break;

      // =============================
      //   ILU(tau,p) Preconditioner
      // =============================
    case BasePrecond::ILUTP:

      // real scalar ILU(tau,p) (works with CRS_Matrix only)
      PRECOND_OBJ( BaseMatrix::DOUBLE, BaseMatrix::SPARSE_NONSYM, ILUTP_Precond<Double> );

      // complex scalar ILU(tau,p) (works with CRS_Matrix only)
      PRECOND_OBJ( BaseMatrix::COMPLEX, BaseMatrix::SPARSE_NONSYM, ILUTP_Precond<Complex> );

      // Test, if preconditioner object could be generated. If not, then
      // this is most likely due to a request for an unsupported version
      // of ILUTP_Precond.
      if ( retVal == NULL ) {
        EXCEPTION( " Generation of ILUTP_Precond object failed! "
            << "Most likely this is due to request for an unsupported "
            << "version, i.e. non CRS_Matrix or non-scalar entries!" );
      }
      break;

      // =========================
      //   ILU(k) Preconditioner
      // =========================
    case BasePrecond::ILUK:

      // real scalar ILU(k) (works with CRS_Matrix only)
      PRECOND_OBJ( BaseMatrix::DOUBLE, BaseMatrix::SPARSE_NONSYM, ILUK_Precond<Double> );

      // complex scalar ILU(k) (works with CRS_Matrix only)
      PRECOND_OBJ( BaseMatrix::COMPLEX, BaseMatrix::SPARSE_NONSYM, ILUK_Precond<Complex> );

      // Test, if preconditioner object could be generated. If not, then
      // this is most likely due to a request for an unsupported version
      // of ILUK_Precond.
      if ( retVal == NULL ) {
        EXCEPTION( " Generation of ILUK_Precond object failed! "
            << "Most likely this is due to request for an unsupported "
            << "version, i.e. non CRS_Matrix or non-scalar entries!" );
      }
      break;

      // ========================
      //   ILDL Preconditioners
      // ========================

      // The following are all sub-variants of the same preconditioner
      // class. We therefore generate for all variants an object of type
      // ILDLPrecond. The latter generates a fitting factoriser internally
      // to reflect the specific variant
    case BasePrecond::ILDL0:
    case BasePrecond::ILDLK:
    case BasePrecond::ILDLTP:
    case BasePrecond::ILDLCN:

      // Generate preconditioner
      PRECOND_OBJ( BaseMatrix::DOUBLE,  BaseMatrix::SPARSE_SYM, ILDLPrecond<Double> );
      PRECOND_OBJ( BaseMatrix::COMPLEX, BaseMatrix::SPARSE_SYM, ILDLPrecond<Complex> );

      // Test, if preconditioner object could be generated. If not, then
      // this is most likely due to a request for an unsupported version.
      if ( retVal == NULL ) {
        std::string tmp;
        tmp = BasePrecond::precondType.ToString(ptype);
        

        EXCEPTION( "Generation of ILDLPrecond object (variant = "
                   << tmp << ") failed! "
                   << "Most likely this is due to request for an unsupported "
                   << "version, i.e. non SCRS_Matrix or non-scalar entries!" );
      }
      break;

      // =======================
      //   IC0 Preconditioner
      // =======================
    case BasePrecond::IC0:

      // real scalar IC(0) (currently only in CRS)
      PRECOND_OBJ( BaseMatrix::DOUBLE, BaseMatrix::SPARSE_SYM, IC0Precond<Double> );

      // complex scalar IC(0) (currently only in CRS)
      PRECOND_OBJ( BaseMatrix::COMPLEX, BaseMatrix::SPARSE_SYM, IC0Precond<Complex> );

      // Test, if preconditioner object could be generated. If not, then
      // this is most likely due to a request for an unsupported version
      // of ICPrecond.
      if ( retVal == NULL ) {
        EXCEPTION( " Request for unsupported version of IC0Precond" );
      }
      break;

      // ============================
      //   Multigrid Preconditioner
      // ============================
    case BasePrecond::MG:

      // multigrid preconditioners can be used with CRS_Matrix only.
      // PRECOND_OBJ( DOUBLE, SPARSE_NONSYM, 1, MGPrecond<Double> );

      // Test, if preconditioner object could be generated. If not, then
      // this is most likely due to a request for an unsupported version
      // of MGPrecond.
      if ( retVal == NULL ) {
        EXCEPTION( " Request for unsupported version of MGPrecond" );
      }
      break;

      // ======================
      //   Something's broken
      // ======================

    default:
      EXCEPTION( "GeneratePrecondObject failed: Preconditioner type unknown" );
    }

    // test, if preconditioner object could be generated
    if ( retVal == NULL ) {
      EXCEPTION( "Something's broke. No preconditioner was generated!" );
    }

    // Finished
    return retVal;
  }

  BaseSBMPrecond* GenerateSBMPrecondObject( const SBM_Matrix &mat,
                                            ParamNode *solverNode,
                                            InfoNode *olasInfo ) {

    BaseSBMPrecond *retVal = NULL;

    std::string precondStr = "";

    ParamNode* solverXML = NULL;
    solverXML = solverNode->Get("solver", false);
  
    EnumMap::iterator it, end;
    it = BasePrecond::precondType.map.begin();
    end = BasePrecond::precondType.map.end();
    
    for( ; it != end; it++ ) {
      if( solverXML->Has("precond", it->second) ) {
        if(precondStr != "")
          EXCEPTION("Two preconditioners have been specified: " << precondStr << " and " << (it->second))
        
        precondStr = it->second;
      }
    }
  
    if(precondStr == "") {
      precondStr = "noPrecond";
    }

    BasePrecond::PrecondType ptype = BasePrecond::NOPRECOND;
    ptype = BasePrecond::precondType.Parse(precondStr);

    // Branch depending on desired preconditioner
    switch( ptype ) {

    // ==============================
    //   Not really preconditioners
    // ==============================

    // Note: If no preconditioner was demanded we will currently also
    //       generate an identity preconditioner to remain consistent.
    case BasePrecond::NOPRECOND:
    case BasePrecond::ID:
      retVal = new IdPrecondSBM;
      (*cla) << " GenerateStdPrecondObject: Generated Identity preconditioner"
      << std::endl;
      break;
    default:
      EXCEPTION( "GeneratePrecondObject failed: Preconditioner type unknown" );
    }

    // test, if preconditioner object could be generated
    if ( retVal == NULL ) {
      EXCEPTION( "Something's broke. No preconditioner was generated!" );
    }

    // Finished
    return retVal;

  }
}
