// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <def_use_lapack.hh>
#include <def_use_pardiso.hh>
#include <def_use_cholmod.hh>

#include "OLAS/precond/idprecond.hh"
#include "OLAS/precond/generateprecond.hh"
#include "MatVec/basematrix.hh"


#ifdef USE_PARDISO
#include "OLAS/external/pardiso/pardisosolver.hh"
#endif

#ifdef USE_CHOLMOD
#include "OLAS/external/cholmod/CholMod.hh"
#endif
#include "OLAS/solver/ldlsolver.hh"

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
#include "MatVec/vbr_matrix.hh"
#include "MatVec/sbmmatrix.hh"

namespace CoupledField {


  // *********************************************************
  //   Central macro for generation of preconditioner object
  //   Its use also takes care of the instantiation in the
  //   case of the templated preconditioners
  // *********************************************************
#define PRECOND_OBJ(matEntry,matStore,precondObjType)\
if ((entryType==matEntry) && (storageType==matStore) ) {\
retVal = new precondObjType( mat, solverXML, olasInfo );\
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
  
  // Block Jacobi preconditioner for VBR matrix
  typedef BlockJacPrecond<VBR_Matrix<Double>,Double> RVBRBlockJac;
  typedef BlockJacPrecond<VBR_Matrix<Complex>,Complex> CVBRBlockJac;


  // **********************
  //   Generation routine
  // **********************
  BasePrecond* GenerateStdPrecondObject( const StdMatrix &mat,
                                            PtrParamNode systemNode,
                                            PtrParamNode olasInfo ) {


    BasePrecond *retVal = NULL;
    
    // Obtain matrix type information
    BaseMatrix::EntryType entryType = mat.GetEntryType();
    BaseMatrix::StorageType storageType = mat.GetStorageType();
    std::string precondStr = "";

    PtrParamNode solverXML;
    solverXML = systemNode->Get("solver", ParamNode::INSERT );

    EnumMap::iterator it, end;
    it = BasePrecond::precondType.map.begin();
    end = BasePrecond::precondType.map.end();
    
    for( ; it != end; it++ ) {
      if( solverXML->HasByVal("precond", it->second) ) {
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
      
      // ===============================
      //   Block Jacobi Preconditioner
      // ===============================
    case BasePrecond::BLOCK_JACOBI:
      
      // real VBR Block Jacobi
      PRECOND_OBJ(BaseMatrix::DOUBLE, BaseMatrix::VAR_BLOCK_ROW, RVBRBlockJac);
      
      // complex VBR Block Jacobi
      PRECOND_OBJ(BaseMatrix::COMPLEX, BaseMatrix::VAR_BLOCK_ROW, CVBRBlockJac);

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
      
      
      // ============================
      //   Pardiso Preconditioner
      // ============================
    case BasePrecond::PARDISO:
#ifdef USE_PARDISO
      // Check suitability of matrix
      if ( mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX ) {
        EXCEPTION( "PardisoSolver only works with (S)CRS_Matrix class!" );
      }
      else {
        const StdMatrix &stdmat = dynamic_cast<const StdMatrix &>(mat);
        if ( stdmat.GetStorageType() != BaseMatrix::SPARSE_NONSYM &&
            stdmat.GetStorageType() != BaseMatrix::SPARSE_SYM  ) {
          EXCEPTION( "PardisoSolver only works with (S)CRS_Matrix class!" );
        }
      }

      if ( entryType == BaseMatrix::DOUBLE ) {
        retVal = new PardisoSolver<Double>( solverXML, olasInfo );
        ASSERTMEM( retVal, sizeof(PardisoSolver<Double>) );
        (*cla) << " GeneratePrecond: Generated real Pardiso precond"
            << std::endl;
      }
      if ( entryType == BaseMatrix::COMPLEX ) {
        retVal = new PardisoSolver<Complex>( solverXML, olasInfo );
        ASSERTMEM( retVal, sizeof(PardisoSolver<Complex>) );
        (*cla) << " GeneratePrecond: Generated complex Pardiso precond"
            << std::endl;
      }
#else

      EXCEPTION( "Compile with USE_PARDISO to enable interface to Pardiso "
          "library" );
#endif
      break;

   
      
      // ============================
      //   Pardiso Preconditioner
      // ============================
    case BasePrecond::CHOLMOD:
#ifdef USE_CHOLMOD
  {
    if(mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX || 
        dynamic_cast<const StdMatrix &>(mat).GetStorageType() != BaseMatrix::SPARSE_SYM){
      EXCEPTION("CholMod only works with SCRS_Matrix class!");
    }
    if(entryType == BaseMatrix::DOUBLE){
      retVal = new CholMod<Double>(solverXML, olasInfo, entryType);
      (*cla) << " GenerateSolver: Generated real CholMod solver" << std::endl;
    }else if(entryType == BaseMatrix::COMPLEX){
      retVal = new CholMod<Complex>(solverXML, olasInfo, entryType);
      (*cla) << " GenerateSolver: Generated complex CholMod solver" << std::endl;
    }
  }
#else
    EXCEPTION("Compile with USE_CHOLMOD to enable interface to CholMod");
#endif
    break;
    
    case BasePrecond::LDL_SOLVER:
      if(mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX 
          || dynamic_cast<const StdMatrix &>(mat).GetStorageType() != BaseMatrix::SPARSE_SYM){
         EXCEPTION("DirectLDL solver only works with SCRS_Matrix class!");
       }
       if(entryType == BaseMatrix::DOUBLE){
         retVal = new LDLSolver<Double>(solverXML, olasInfo );
         (*cla) << " GenerateSolver: Generated real CholMod solver" << std::endl;
       }else if(entryType == BaseMatrix::COMPLEX){
         retVal = new LDLSolver<Complex>(solverXML, olasInfo);
         (*cla) << " GenerateSolver: Generated complex CholMod solver" << std::endl;
       }
       break;
       
       // ======================
       //   Something's broken
       // ======================
    case BasePrecond::LDL_SOLVER2:
    case BasePrecond::DIRECT:
    case BasePrecond::RICHARDSON:
    case BasePrecond::CG:
    case BasePrecond::LANCZOS:
    case BasePrecond::QMR:
    case BasePrecond::GMRES:
    case BasePrecond::MINRES:
    case BasePrecond::SYMMLQ:
    case BasePrecond::LAPACK_LU: 
    case BasePrecond::LAPACK_LL:
    case BasePrecond::ILUPACK:
    case BasePrecond::LU_SOLVER: 
    case BasePrecond::DIAGSOLVER:
      EXCEPTION("Missing cases for solvers as preconditioners");
      break;
   
    default:
      EXCEPTION( "GeneratePrecondObject failed: Preconditioner type unknown" );
      break;
    }

    // test, if preconditioner object could be generated
    if ( retVal == NULL ) {
      EXCEPTION( "Something's broke. No preconditioner was generated!" );
    }

    // Finished
    return retVal;
  }

  BaseSBMPrecond* GenerateSBMPrecondObject( const SBM_Matrix &mat,
                                            PtrParamNode solverNode,
                                            PtrParamNode olasInfo ) {

    BaseSBMPrecond *retVal = NULL;

    std::string precondStr = "";

    PtrParamNode solverXML;
    solverXML = solverNode->Get("solver", ParamNode::INSERT );
  
    EnumMap::iterator it, end;
    it = BasePrecond::precondType.map.begin();
    end = BasePrecond::precondType.map.end();
    
    for( ; it != end; it++ ) {
      if( solverXML->HasByVal("precond", it->second) ) {
        if(precondStr != "")
          EXCEPTION("Two preconditioners have been specified: " << precondStr << " and " << (it->second))
        
        precondStr = it->second;
      }
    }
  
    if(precondStr == "") {
      precondStr = "noPrecond";
    }
    
    UInt numBlocks = mat.GetNumRows();

    BasePrecond::PrecondType ptype = BasePrecond::NOPRECOND;
    ptype = BasePrecond::precondType.Parse(precondStr);

    // Branch depending on desired preconditioner
    BasePrecond * p = NULL;
    PtrParamNode pardisoNode(new ParamNode(ParamNode::INSERT));
    PtrParamNode blockJacobiNode(new ParamNode(ParamNode::INSERT));
    PtrParamNode infoNode;
    switch( ptype ) {

    // ==============================
    //   Not really preconditioners
    // ==============================

    // Note: If no preconditioner was demanded we will currently also
    //       generate an identity preconditioner to remain consistent.
    case BasePrecond::NOPRECOND:
    case BasePrecond::ID:
//      retVal = new IdPrecondSBM;
//      (*cla) << " GenerateStdPrecondObject: Generated Identity preconditioner"
//      << std::endl;
//      break;
    case BasePrecond::JACOBI:
      retVal = new BaseSBMPrecond(numBlocks, olasInfo);
       infoNode = retVal->GetInfoNode();
     
      // ========================
      // MAGNETIC CASE
      // ========================
      // Block #0: Pardiso solver
      
      pardisoNode->Get("solver")->Get("precond")->SetValue("pardiso");
      p = GenerateStdPrecondObject(mat(0,0), pardisoNode, infoNode );
      retVal->SetPrecond(0, p);
      
      // Block #1: Block-Jacobi-preconditioner
      blockJacobiNode->Get("solver")->Get("precond")->SetValue("BlockJacobi");
      //p = GenerateStdPrecondObject(mat(1,1), pardisoNode, olasInfo );
      p = GenerateStdPrecondObject(mat(1,1), blockJacobiNode, infoNode );
      retVal->SetPrecond(1, p);
      
      // if we have static condensation, we have no third row
      if( mat.GetNumRows() == 3 ) {
        // Block #1: Block-Jacobi preconditioner
        //p = GenerateStdPrecondObject(mat(2,2), pardisoNode, olasInfo );
        p = GenerateStdPrecondObject(mat(2,2), blockJacobiNode, infoNode );
        retVal->SetPrecond(2, p);
      }
      
      
      // ========================
      // Standard Case
      // ========================
//      // hard-coded: set for each block the specified standard preconditioner
//      for(UInt i = 0; i < numBlocks; ++i ) {
//        BasePrecond * p = GenerateStdPrecondObject(mat(i,i), solverNode, olasInfo );
//        retVal->SetPrecond(i, p);
//      }
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

// ========================
