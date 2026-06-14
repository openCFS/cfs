#include <def_use_pardiso.hh>
#include <def_use_suitesparse.hh>

#include "OLAS/precond/IdPrecondStd.hh"
#include "OLAS/precond/generateprecond.hh"
#include "MatVec/BaseMatrix.hh"
#include "MatVec/CRS_Matrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "MatVec/Diag_Matrix.hh"
#include "MatVec/VBR_Matrix.hh"
#include "MatVec/SBM_Matrix.hh"
#include "OLAS/algsys/SolStrategy.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#ifdef USE_PARDISO
#include "OLAS/external/pardiso/PardisoSolver.hh"
#endif

#ifdef USE_SUITESPARSE
#include "OLAS/external/cholmod/CholMod.hh"
#endif

#include "OLAS/solver/LDLSolver.hh"

// include source code for templated preconditioners
#include "OLAS/precond/JacPrecond.hh"
#include "OLAS/precond/SSORPrecond.hh"
#include "OLAS/precond/ILU0Precond.hh"
#include "OLAS/precond/ILUTP_Precond.hh"
#include "OLAS/precond/ILUK_Precond.hh"
#include "OLAS/precond/MGPrecond.hh"
#include "OLAS/precond/ILDLPrecond/ildlprecond.hh"
#include "OLAS/precond/IC0Precond.hh"
#include "OLAS/precond/SBMDiagPrecond.hh"
#include "OLAS/precond/SBMJacobiPrecond.hh"
#include "OLAS/precond/LDLPrecond.hh"

namespace CoupledField {

// define logging stream
DEFINE_LOG(genPrecond, "genPrecond")



  // *********************************************************
  //   Central macro for generation of preconditioner object
  //   Its use also takes care of the instantiation in the
  //   case of the templated preconditioners
  // *********************************************************
#define PRECOND_OBJ(matEntry,matStore,precondObjType)\
if ((entryType==matEntry) && (storageType==matStore) ) {\
retVal = new precondObjType( mat, precondNode, olasInfo );\
LOG_DBG(genPrecond) << " GenerateStdPrecondObject: Generated "\
       << MACRO2STRING(precondObjType)\
       << " preconditioner"; }


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
  
  // Block Jacobi preconditioner 
  typedef BlockJacPrecond<SCRS_Matrix<Double>,Double> RSCRSBlockJac;
  typedef BlockJacPrecond<SCRS_Matrix<Complex>,Complex> CSCRSBlockJac;
  typedef BlockJacPrecond<VBR_Matrix<Double>,Double> RVBRBlockJac;
  typedef BlockJacPrecond<VBR_Matrix<Complex>,Complex> CVBRBlockJac;


  // **********************
  //   Generation routine
  // **********************
  BasePrecond* GenerateStdPrecondObject( const StdMatrix &mat,
                                         const std::string& precondId,
                                         PtrParamNode precondList,
                                         PtrParamNode olasInfo ) {


    BasePrecond *retVal = NULL;
    
    
    // Try to get xml node based on precond id
    ParamNodeList pNodes =  precondList->GetChildren();
    PtrParamNode precondNode;
    for( UInt i = 0; i < pNodes.GetSize(); ++i ) {
      if( pNodes[i]->Get("id")->As<std::string>() == precondId ) {
        precondNode = pNodes[i]; 
      }
    }
    
    if(!precondNode) {
      EXCEPTION("Preconditioner with id '" << precondId 
                << "' was not found!");
    }
    
    // Obtain matrix type information
    BaseMatrix::EntryType entryType = mat.GetEntryType();
    BaseMatrix::StorageType storageType = mat.GetStorageType();
    std::string precondStr = "";

    EnumMap::iterator it, end;
    it = BasePrecond::precondType.map.begin();
    end = BasePrecond::precondType.map.end();
    
    for( ; it != end; it++ ) {
      if( precondNode->GetName() ==  it->second ) {
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
      retVal = new IdPrecondStd(precondNode, olasInfo);
      LOG_DBG(genPrecond) << " GenerateStdPrecondObject: Generated Identity preconditioner";
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
      PRECOND_OBJ(BaseMatrix::COMPLEX,BaseMatrix::SPARSE_SYM,CSCRSJacDof1);
      
      //real diag jacobi
      PRECOND_OBJ(BaseMatrix::DOUBLE,BaseMatrix::DIAG,RDIAGJacDof1);
      //complex diag jacobi
      PRECOND_OBJ(BaseMatrix::COMPLEX,BaseMatrix::DIAG,CDIAGJacDof1);
      

      break;
      
      // ===============================
      //   Block Jacobi Preconditioner
      // ===============================
    case BasePrecond::BLOCK_JACOBI:

      // real SCRS Block Jacobi
      PRECOND_OBJ(BaseMatrix::DOUBLE, BaseMatrix::SPARSE_SYM, RSCRSBlockJac);
      
      // complex SCRS Block Jacobi
      PRECOND_OBJ(BaseMatrix::COMPLEX, BaseMatrix::SPARSE_SYM, CSCRSBlockJac);

      
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
      PRECOND_OBJ( BaseMatrix::DOUBLE, BaseMatrix::SPARSE_NONSYM, MGPrecond<Double> );
      PRECOND_OBJ( BaseMatrix::COMPLEX, BaseMatrix::SPARSE_NONSYM, MGPrecond<Complex> );

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
    case BasePrecond::PARDISO_PRECOND:
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
        retVal = new PardisoSolver<Double>( precondNode, olasInfo );
        ASSERTMEM( retVal, sizeof(PardisoSolver<Double>) );
        LOG_DBG(genPrecond) << " GeneratePrecond: Generated real Pardiso precond";
      }
      
      if ( entryType == BaseMatrix::COMPLEX ) {
        retVal = new PardisoSolver<Complex>( precondNode, olasInfo );
        ASSERTMEM( retVal, sizeof(PardisoSolver<Complex>) );
        LOG_DBG(genPrecond) << " GeneratePrecond: Generated complex Pardiso precond";
      }
#else

      EXCEPTION( "Compile with USE_PARDISO to enable interface to Pardiso "
          "library" );
#endif
      break;
      
      // ============================
      //   LDL Preconditioner
      // ============================
    case BasePrecond::LDL:
#ifdef USE_PARDISO
      if ( mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX ) {
        EXCEPTION( "LDLPrecond only works with (S)CRS_Matrix class!" );
      }
      else {
        const StdMatrix &stdmat = dynamic_cast<const StdMatrix &>(mat);
        if ( stdmat.GetStorageType() != BaseMatrix::SPARSE_SYM ) {
          EXCEPTION( "LDLPrecond requires a sparseSym (SCRS) matrix!" );
        }
      }

      if ( entryType == BaseMatrix::DOUBLE ) {
        retVal = new LDLPrecond<Double>( precondNode, olasInfo );
        ASSERTMEM( retVal, sizeof(LDLPrecond<Double>) );
        LOG_DBG(genPrecond) << " GeneratePrecond: Generated real LDL precond";
      }
      if ( entryType == BaseMatrix::COMPLEX ) {
        retVal = new LDLPrecond<Complex>( precondNode, olasInfo );
        ASSERTMEM( retVal, sizeof(LDLPrecond<Complex>) );
        LOG_DBG(genPrecond) << " GeneratePrecond: Generated complex LDL precond";
      }
#else
      EXCEPTION( "Compile with USE_PARDISO to enable the LDL preconditioner" );
#endif
      break;

      // ============================
      //   Cholmod Preconditioner
      // ============================
    case BasePrecond::CHOLMOD:
#ifdef USE_SUITESPARSE
  {
    if(mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX || 
        dynamic_cast<const StdMatrix &>(mat).GetStorageType() != BaseMatrix::SPARSE_SYM){
      EXCEPTION("CholMod only works with SCRS_Matrix class!");
    }
    if(entryType == BaseMatrix::DOUBLE){
      retVal = new CholMod<Double>(precondNode, olasInfo, entryType);
      LOG_DBG(genPrecond) << " GenerateSolver: Generated real CholMod solver";
    }else if(entryType == BaseMatrix::COMPLEX){
      retVal = new CholMod<Complex>(precondNode, olasInfo, entryType);
      LOG_DBG(genPrecond) << " GenerateSolver: Generated complex CholMod solver";
    }
  }
#else
    EXCEPTION("Compile with USE_SUITESPARSE to enable interface to CholMod");
#endif
    break;
    
    case BasePrecond::LDL_SOLVER:
      if(mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX 
          || dynamic_cast<const StdMatrix &>(mat).GetStorageType() != BaseMatrix::SPARSE_SYM){
         EXCEPTION("DirectLDL solver only works with SCRS_Matrix class!");
       }
       if(entryType == BaseMatrix::DOUBLE){
         retVal = new LDLSolver<Double>(precondNode, olasInfo );
         LOG_DBG(genPrecond) << " GenerateSolver: Generated real CholMod solver";
       }else if(entryType == BaseMatrix::COMPLEX){
         retVal = new LDLSolver<Complex>(precondNode, olasInfo);
         LOG_DBG(genPrecond) << " GenerateSolver: Generated complex CholMod solver";
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


    // Call post-initialization method
    retVal->BasePrecond::PostInit();
    
    // Finished
    return retVal;
  }

  BaseSBMPrecond* GenerateSBMPrecondObject( const SBM_Matrix &mat,
                                            const std::string& precondId,
                                            PtrParamNode precondList,
                                            PtrParamNode olasInfo ) {

    BaseSBMPrecond *retVal = NULL;

    // Try to get xml node based on precond id
    ParamNodeList pNodes =  precondList->GetChildren();
    PtrParamNode precondNode;
    for( UInt i = 0; i < pNodes.GetSize(); ++i ) {
      if( pNodes[i]->Get("id")->As<std::string>() == precondId ) {
        precondNode = pNodes[i]; 
      }
    }
    if(!precondNode) {
      EXCEPTION("Preconditioner with id '" << precondId 
                << "' was not found!");
    }
    
  
    EnumMap::iterator it, end;
    it = BasePrecond::precondType.map.begin();
    end = BasePrecond::precondType.map.end();
    std::string precondStr = "";
    for( ; it != end; it++ ) {
      if( precondNode->GetName() ==  it->second ) {
        if(precondStr != "")
          EXCEPTION("Two preconditioners have been specified: " 
              << precondStr << " and " << (it->second))
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
    PtrParamNode infoNode, blockInfoNode;
    SBMDiagPrecond * dp = NULL;
    SBMJacobiPrecond * jp = NULL;
    
    switch( ptype ) {

    // ==============================
    //   Not really preconditioners
    // ==============================

    // Note: If no preconditioner was demanded we will currently also
    //       generate an identity preconditioner to remain consistent.
    case BasePrecond::NOPRECOND:
    case BasePrecond::ID:
      WARN("In this case we create a SBM-Diag-Precond which does nothing");
//      retVal = new IdPrecondSBM;
//      LOG_DBG(genPrecond) << " GenerateStdPrecondObject: Generated Identity preconditioner";
     break;
      
     // ============================
     //   SBM-Jacobi Preconditioner
     // ============================
    case BasePrecond::SBM_JACOBI:
    {
      jp = new SBMJacobiPrecond(numBlocks, olasInfo);
      infoNode = jp->GetInfoNode();

      // Loop over all blocks and assign to each diagonal block a
      // StdPrecond object with the same preconditioner, otherwise we
      // would have a SBMDiagPrecond
      std::string precondId;
      if(precondNode->Get("precond")->Has("id")){
        precondId = precondNode->Get("precond")->Get("id")->As<std::string>();
      }else{
        EXCEPTION("No preconditioner id was defined for the SBM_Jacobi preconditioner!!");
      }

      for( UInt i = 0; i < numBlocks; ++i ) {
        // just for the info-xml
        blockInfoNode = infoNode->Get("block",ParamNode::APPEND);
        blockInfoNode->Get("num")->SetValue(i+1);
        // generate preconditioner object
        BasePrecond * actPrecond = GenerateStdPrecondObject( mat(i,i), precondId, precondList, blockInfoNode );
        // pass precond object to sbm-preconditioner
        jp->SetPrecond(i, actPrecond );
      }

      retVal = jp;
    }
      break;


      // ============================
      //   SBM-Diag Preconditioner
      // ============================ 
    case BasePrecond::SBM_DIAG:
      dp = new SBMDiagPrecond(numBlocks, olasInfo);
      infoNode = dp->GetInfoNode();

      // Loop over all blocks and assign to each diagonal block a StdPrecond object 
      for( UInt i = 0; i < numBlocks; ++i ) {

        blockInfoNode = infoNode->Get("block",ParamNode::APPEND);
        blockInfoNode->Get("num")->SetValue(i+1);
        
        // get hold of preconditioner subnode for given block
        PtrParamNode actNode = 
            precondNode->GetByVal( "precond", "block", std::to_string(i+1));
        if( !actNode ) {
          EXCEPTION("No preconditioner was defined for block #" << i+1 );
        }
        std::string precondId = actNode->Get("id")->As<std::string>();

        // generate preconditioner object
        BasePrecond * actPrecond = 
            GenerateStdPrecondObject( mat(i,i), precondId, precondList,
                                      blockInfoNode );

        // pass precond object to sbm-preconditioner
        dp->SetPrecond(i, actPrecond );

      }
      retVal = dp;
      break;
       
    default:
      EXCEPTION( "GeneratePrecondObject failed: Preconditioner type unknown" );
      break;
    }

    // test, if preconditioner object could be generated
    if ( retVal == NULL ) {
      EXCEPTION( "Something's broke. No preconditioner was generated!" );
    }

    // Call post-initialization method
    retVal->PostInit();
    
    // Finished
    return retVal;

  }
  
  
  std::set<BaseMatrix::StorageType> 
   GetPrecondCompatMatrixFormats(BasePrecond::PrecondType pt) {
    
    std::set<BaseMatrix::StorageType> ret;
    switch(pt) {
      
      case BasePrecond::ID:
        ret.insert(BaseMatrix::SPARSE_SYM);
        ret.insert(BaseMatrix::SPARSE_NONSYM);
        ret.insert(BaseMatrix::DIAG);
        ret.insert(BaseMatrix::VAR_BLOCK_ROW);
        break;
        
      case BasePrecond::MG:
        ret.insert(BaseMatrix::SPARSE_NONSYM);
        break;
      
      case BasePrecond::JACOBI:
        ret.insert(BaseMatrix::SPARSE_SYM);
        ret.insert(BaseMatrix::SPARSE_NONSYM);
        ret.insert(BaseMatrix::DIAG);
        break;

      case BasePrecond::BLOCK_JACOBI :
        ret.insert(BaseMatrix::SPARSE_SYM);
        ret.insert(BaseMatrix::VAR_BLOCK_ROW);
        break;

      case BasePrecond::SSOR:
        ret.insert(BaseMatrix::SPARSE_NONSYM);
        break;

      case BasePrecond::ILU0:
      case BasePrecond::ILUTP:
      case BasePrecond::ILUK:
        ret.insert(BaseMatrix::SPARSE_NONSYM);
        break;

      case BasePrecond::ILDL0:
      case BasePrecond::ILDLK:
      case BasePrecond::ILDLTP:
      case BasePrecond::ILDLCN:
        ret.insert(BaseMatrix::SPARSE_SYM);
        break;
            
      case BasePrecond::IC0 :
        ret.insert(BaseMatrix::SPARSE_SYM);
        break;
        
      case BasePrecond::SBM_DIAG:
      case BasePrecond::SBM_JACOBI:
        // works with all SBM matrices
        break;
            
        // Solver Based preconditioners
      case BasePrecond::PARDISO_PRECOND :
        ret.insert(BaseMatrix::SPARSE_SYM);
        ret.insert(BaseMatrix::SPARSE_NONSYM);
        break;

      case BasePrecond::CHOLMOD:
        ret.insert(BaseMatrix::SPARSE_SYM);
        break;
      
      case BasePrecond::LDL:
        ret.insert(BaseMatrix::SPARSE_SYM);
        break;
      
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
      case BasePrecond::LU_SOLVER: 
      case BasePrecond::DIAGSOLVER:
        EXCEPTION("Missing case statement for solvers as preconditioners");
        break;
      default:
        EXCEPTION("Missing case statement");
        break;
    } 
    return ret;
  }

  bool IsPreconndSBMCapable(BasePrecond::PrecondType pt) {
    bool ret = false;
    switch(pt) {

      case BasePrecond::JACOBI:
        ret = true;
        break;

      default:
        ret = false;
        break;
    }
    return ret;
   }

}

// ========================
