#include <def_use_pardiso.hh>
#include <def_use_suitesparse.hh>
#include <def_use_superlu.hh>
#include <def_use_lis.hh>
#include <def_use_ginkgo.hh>
#include <def_use_petsc.hh>
#include <def_use_phist_cg.hh>

#ifdef USE_GINKGO
#include "OLAS/external/ginkgo/GinkgoSolver.hh"
#endif

#include "OLAS/algsys/SolStrategy.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "generatesolver.hh"
#include "General/Exception.hh"

#include "OLAS/external/lapack/Lapack_LU.hh"
#include "OLAS/external/lapack/Lapack_LL.hh"

#ifdef USE_PARDISO
#include "OLAS/external/pardiso/PardisoSolver.hh"
#include "OLAS/external/pardiso/PardisoSolverPrimitive.hh"
#endif

#ifdef USE_SUITESPARSE
#include "OLAS/external/cholmod/CholMod.hh"
#include "OLAS/external/umfpack/UMFPACKSolver.hh"
#endif

#ifdef USE_SUPERLU
#include "OLAS/external/superlu/SuperLUSolver.hh"
#endif

#ifdef USE_LIS
#include "OLAS/external/lis/LISSolver.hh"
#endif

#ifdef USE_PETSC
#include "OLAS/external/petsc/PETSCSolver.hh"
#endif

#ifdef USE_PHIST_CG
#include "OLAS/external/phist/PhistLinearSolver.hh"
#endif

// include source code for templated solvers
#include "BaseSolver.hh"
#include "RichardsonSolver.hh"
#include "CGSolver.hh"
#include "GMRESSolver.hh"
#include "MINRESSolver.hh"
#include "LUSolver.hh"
#include "LDLSolver.hh"
#include "DiagSolver.hh"
#include "ExternalSolver.hh"

namespace CoupledField {

// define logging stream
DEFINE_LOG(genSolver, "genSolver")

// ****************************
//   Generate a solver object
// ****************************
BaseSolver* GenerateSolverObject( const BaseMatrix &mat,
                                  shared_ptr<SolStrategy> strat,
                                  PtrParamNode xml,
                                  PtrParamNode  olasInfo ){

  BaseSolver *retSolver = NULL;
  BaseMatrix::EntryType eType = mat.GetEntryType();
  std::string solverStr = "";

  
  // Obtain current solver id from strategy object;
  std::string solverId = strat->GetSolverId();
  ParamNodeList sNodes =  xml->GetChildren();
  PtrParamNode solverNode;
  for( UInt i = 0; i < sNodes.GetSize(); ++i ) {
    if( sNodes[i]->Get("id")->As<std::string>() == solverId ) {
      solverNode = sNodes[i]; 
    }
  }

  if(!solverNode) {
    EXCEPTION("Solver with id '" << solverId << "' was not found!");
  }

  // Convert string to enum object
  EnumMap::iterator it, end;
  it = BaseSolver::solverType.map.begin();
  end = BaseSolver::solverType.map.end();

  for( ; it != end; it++ ) {
    if( solverNode->GetName() ==  it->second ) {
      if(solverStr != "")
        EXCEPTION("Two solvers have been specified: " << solverStr << " and " << (it->second))

        solverStr = it->second;
    }
  }
  
  if(solverStr == "") {
    EXCEPTION("Could not determine solver!")
  }

  BaseSolver::SolverType solver = BaseSolver::NOSOLVER;
  solver = BaseSolver::solverType.Parse(solverStr);
  
  // Branch depending on desired solver
  switch( solver ) {
  case BaseSolver::RICHARDSON:
    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new RichardsonSolver<Double>( solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(RichardsonSolver<Double>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated real Richardson solver";
    }
    else if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new RichardsonSolver<Complex>( solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(RichardsonSolver<Complex>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated complex Richardson solver";
    }
    break;
  case BaseSolver::DIAGSOLVER:
    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new DiagSolver<Double>(solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(DiagSolver<Double>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated real Diagonal solver";
    }
    else if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new DiagSolver<Complex>(solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(DiagSolver<Complex>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated complex Diagonal solver";
    }
    break;

  case BaseSolver::CG:
    if(eType == BaseMatrix::DOUBLE) {
      retSolver = new CGSolver<Double>(solverNode, olasInfo );
      LOG_DBG(genSolver) << " GenerateSolver: Generated real CG solver" ;
    }
    if(eType == BaseMatrix::COMPLEX) {
      retSolver = new CGSolver<Complex>(solverNode, olasInfo );
      LOG_DBG(genSolver) << " GenerateSolver: Generated complex CG solver";
    }
    break;

  case BaseSolver::GMRES:
    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new GMRESSolver<Double>(solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(GMRESSolver<Double>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated real GMRES solver";
    }
    else if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new GMRESSolver<Complex>(solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(GMRESSolver<Complex>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated complex GMRES solver";
    }
    break;

  case BaseSolver::MINRES:
    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new MINRESSolver<Double>( solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(MINRESSolver<Double>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated real MINRES solver";
    }
    else if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new MINRESSolver<Complex>( solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(MINRESSolver<Complex>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated complex MINRES solver";
    }
    break;

  case BaseSolver::LU_SOLVER:
    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new LUSolver<Double>( solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(LUSolver<Double>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated real "
      << BaseSolver::solverType.ToString(solver);
    }
    else if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new LUSolver<Complex>( solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(LUSolver<Complex>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated complex "
                         << BaseSolver::solverType.ToString(solver);
    }
    break;
  case BaseSolver::LDL_SOLVER:
    // Get block size of matrix entries (we need an StdMatrix for this,
    // but, if it's not, an LDL_Solver makes no sense anyhow)
  {
    if(eType == BaseMatrix::DOUBLE) {
      retSolver = new LDLSolver<Double>( solverNode, olasInfo );
      LOG_DBG(genSolver) << " GenerateSolver: Generated LDLSolver<Double>";
    }
    else if(eType == BaseMatrix::COMPLEX)  {
      retSolver = new LDLSolver<Complex>( solverNode, olasInfo );
      LOG_DBG(genSolver) << " GenerateSolver: Generated LDLSolver<COMPLEX>";
    }
  }
  break;


  case BaseSolver::LAPACK_LU:
    if ( mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX ) {
      EXCEPTION( "LAPACK_LU only works with a LAPACK_GBMATRIX!");
    }
    else {
      const StdMatrix &stdmat = dynamic_cast<const StdMatrix &>(mat);
      if ( stdmat.GetStorageType() != BaseMatrix::LAPACK_GBMATRIX ) {
        EXCEPTION( "LAPACK_LU only works with a LAPACK_GBMATRIX!");
      }
      else {
        retSolver = new Lapack_LU( solverNode, olasInfo );
        ASSERTMEM( retSolver, sizeof(Lapack_LU) );
        LOG_DBG(genSolver) << " GenerateSolver: Generated Lapack_LU solver";
      }
    }
    break;

  case BaseSolver::LAPACK_LL:
    if ( mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX ) {
      EXCEPTION( "LAPACK_LL only works with an SCRS_MATRIX!" );
    }
    else {
      const StdMatrix &stdmat = dynamic_cast<const StdMatrix &>(mat);
      if ( stdmat.GetStorageType() != BaseMatrix::SPARSE_SYM ) {
        EXCEPTION( "LAPACK_LL only works with an SCRS_MATRIX!" );
      }
      else {
        retSolver = new Lapack_LL( solverNode, olasInfo );
        ASSERTMEM( retSolver, sizeof(Lapack_LL) );
        LOG_DBG(genSolver) << " GenerateSolver: Generated Lapack_LL solver";
      }
    }
    break;

  case BaseSolver::EXTERNAL_SOLVER:
      if ( eType == BaseMatrix::DOUBLE ) {
        retSolver = new ExternalSolver<Double>( solverNode, olasInfo );
        ASSERTMEM( retSolver, sizeof(ExternalSolver<Double>) );
        LOG_DBG(genSolver) << " GenerateSolver: Generated real " << BaseSolver::solverType.ToString(solver);
      }
      else if ( eType == BaseMatrix::COMPLEX ) {
        retSolver = new ExternalSolver<Complex>( solverNode, olasInfo );
        ASSERTMEM( retSolver, sizeof(ExternalSolver<Complex>) );
        LOG_DBG(genSolver) << " GenerateSolver: Generated complex " << BaseSolver::solverType.ToString(solver);
      }
      break;

  case BaseSolver::PARDISO_SOLVER:

#ifdef USE_PARDISO

    LOG_DBG(genSolver) << "structure type of matrix: " << mat.GetStructureType();
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

    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new PardisoSolver<Double>( solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(PardisoSolver<Double>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated real Pardiso solver";
    }
    if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new PardisoSolver<Complex>( solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(PardisoSolver<Complex>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated complex Pardiso solver";
    }
#else

    EXCEPTION( "Compile with USE_PARDISO to enable interface to Pardiso "
               "library" );
#endif
    break;

  case BaseSolver::UMFPACK:

#ifdef USE_SUITESPARSE

    // Check suitability of matrix
    if ( mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX ) {
      EXCEPTION( "UMFPACKSolver only works with (S)CRS_Matrix class!" );
    }
    else {
      const StdMatrix &stdmat = dynamic_cast<const StdMatrix &>(mat);
      if ( stdmat.GetStorageType() != BaseMatrix::SPARSE_NONSYM &&
          stdmat.GetStorageType() != BaseMatrix::SPARSE_SYM  ) {
        EXCEPTION( "UMFPACKSolver only works with (S)CRS_Matrix class!" );
      }
    }

    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new UMFPACKSolver<Double>( solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(UMFPACKSolver<Double>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated real UMFPACK solver";
    }
    if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new UMFPACKSolver<Complex>( solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(UMFPACKSolver<Complex>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated complex UMFPACK solver";
    }
#else

    EXCEPTION( "Compile with USE_SUITESPARSE to enable interface to UMFPACK "
               "library" );
#endif
    break;

  case BaseSolver::SUPERLU:

#ifdef USE_SUPERLU

    // Check suitability of matrix
    if ( mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX ) {
      EXCEPTION( "SuperLUSolver only works with (S)CRS_Matrix class!" );
    }
    else {
      const StdMatrix &stdmat = dynamic_cast<const StdMatrix &>(mat);
      if ( stdmat.GetStorageType() != BaseMatrix::SPARSE_NONSYM &&
          stdmat.GetStorageType() != BaseMatrix::SPARSE_SYM  ) {
        EXCEPTION( "SuperLUSolver only works with (S)CRS_Matrix class!" );
      }
    }

    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new SuperLUSolver<Double>( solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(SuperLUSolver<Double>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated real SuperLU solver";
    }
    if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new SuperLUSolver<Complex>( solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(SuperLUSolver<Complex>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated complex SuperLU solver";
    }
#else

    EXCEPTION( "Compile with USE_SUPERLU to enable interface to SuperLU "
               "library" );
#endif
    break;

  case BaseSolver::LIS:
#ifdef USE_LIS
    // GetSolverCompatMatrixFormats() ensures proper non-sym type
    retSolver = new LISSolver(solverNode, olasInfo, eType);
#else
    throw Exception("Compile with USE_LIS to enable interface to LIS.");
#endif
  break;

  case BaseSolver::GINKGO:
#ifdef USE_GINKGO
    // GetSolverCompatMatrixFormats() ensures proper non-sym type
    retSolver = new GinkgoSolver(solverNode, olasInfo, eType);
#else
    throw Exception("Compile with USE_GINKGO to enable interface to ginkgo.");
#endif
  break;

  case BaseSolver::PETSC:

#ifdef USE_PETSC
   {
     // Check suitability of matrix
     if (mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX)
       EXCEPTION("PETSC only works with (S)CRS_Matrix class!");
     const StdMatrix &stdmat = dynamic_cast<const StdMatrix &>(mat);
      if(stdmat.GetStorageType() != BaseMatrix::SPARSE_NONSYM
         && stdmat.GetStorageType() != BaseMatrix::SPARSE_SYM  )
       EXCEPTION("PETSC only works with (S)CRS_Matrix class!");

     retSolver = new PETSCSolver(solverNode, olasInfo, eType);
     
     LOG_DBG(genSolver) << " GenerateSolver: Generated PETSC solver";
   }
 #else
   EXCEPTION("Compile with USE_PETSC to enable interface to PETSC.");
 #endif
   break;

  case BaseSolver::PHIST:

  #ifdef USE_PHIST_CG
    {
      // Check suitability of matrix
      if (mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX)
        EXCEPTION("PETSC only works with (S)CRS_Matrix class!");
      const StdMatrix &stdmat = dynamic_cast<const StdMatrix &>(mat);
       if(stdmat.GetStorageType() != BaseMatrix::SPARSE_NONSYM
          && stdmat.GetStorageType() != BaseMatrix::SPARSE_SYM  )
        EXCEPTION("PETSC only works with (S)CRS_Matrix class!");

      retSolver = new PhistLinearSolver(solverNode, olasInfo, eType);

      LOG_DBG(genSolver) << " GenerateSolver: Generated PETSC solver";
    }
  #else
    EXCEPTION("Compile with USE_PHIST_CG to enable interface to PHIST_CG.");
  #endif
    break;

  case BaseSolver::CHOLMOD:
#ifdef USE_SUITESPARSE
  {
    if(mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX || dynamic_cast<const StdMatrix &>(mat).GetStorageType() != BaseMatrix::SPARSE_SYM){
      EXCEPTION("CholMod only works with SCRS_Matrix class!");
    }
    if(eType == BaseMatrix::DOUBLE){
      retSolver = new CholMod<Double>(solverNode, olasInfo, eType);
      LOG_DBG(genSolver) << " GenerateSolver: Generated real CholMod solver";
    }else if(eType == BaseMatrix::COMPLEX){
      retSolver = new CholMod<Complex>(solverNode, olasInfo, eType);
      LOG_DBG(genSolver) << " GenerateSolver: Generated complex CholMod solver";
    }
  }
#else
    EXCEPTION("Compile with USE_SUITESPARSE to enable interface to CholMod");
#endif
    break;

  default:
    EXCEPTION("GenerateSolver: Request for unknown solver type!");
  }

  // Check for unsupported matrix entry type
  if (retSolver == NULL ) EXCEPTION("unhandled type " << eType);

  retSolver->PostInit();

  return retSolver;
}



// *******************************************************
//   Generate a direct solver object for coarse system AMG
// *******************************************************
BaseSolver* GenerateDirSolverObjectAMG( const BaseMatrix &mat,
                                        PtrParamNode  olasInfo ){

  BaseSolver *retSolver = NULL;
  BaseMatrix::EntryType eType = mat.GetEntryType();
  std::string solverStr = "";

  // Hardcoded ... do we need other direct solvers (LU, ...) ????
  //BaseSolver::SolverType solver = BaseSolver::PARDISO_SOLVER;

#ifdef USE_PARDISO
    LOG_DBG(genSolver) << "structure type of matrix: " << mat.GetStructureType();
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
    PtrParamNode solverNode;
    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new PardisoSolverPrimitive<Double>( solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(PardisoSolverPrimitive<Double>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated real Pardiso solver";
    }
    if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new PardisoSolverPrimitive<Complex>( solverNode, olasInfo );
      ASSERTMEM( retSolver, sizeof(PardisoSolverPrimitive<Complex>) );
      LOG_DBG(genSolver) << " GenerateSolver: Generated complex Pardiso solver";
    }
#else
    EXCEPTION( "Compile with USE_PARDISO to enable interface to Pardiso "
               "library" );
#endif

  // Check for unsupported matrix entry type
  if (retSolver == NULL ) EXCEPTION("unhandled type " << eType);

  return retSolver;
}


std::set<BaseMatrix::StorageType> 
GetSolverCompatMatrixFormats(BaseSolver::SolverType st) {
  std::set<BaseMatrix::StorageType> ret;
  switch(st) {
    case BaseSolver::NOSOLVER:
      break;
      
    // All iterative solvers can handle any sparse matrix format
    case BaseSolver::RICHARDSON: 
    case BaseSolver::CG: 
    case BaseSolver::LANCZOS: 
    case BaseSolver::QMR:
    case BaseSolver::GMRES:
    case BaseSolver::MINRES:
    case BaseSolver::SYMMLQ:
      break;

    case BaseSolver::LAPACK_LU:
    case BaseSolver::LAPACK_LL:
      ret.insert(BaseMatrix::LAPACK_GBMATRIX);
      break;

    case BaseSolver::PARDISO_SOLVER:
      ret.insert(BaseMatrix::SPARSE_SYM);
      ret.insert(BaseMatrix::SPARSE_NONSYM);
      break;

    case BaseSolver::UMFPACK:
      ret.insert(BaseMatrix::SPARSE_NONSYM);
      break;

    case BaseSolver::LU_SOLVER:
      ret.insert(BaseMatrix::SPARSE_NONSYM);
      break;

    case BaseSolver::CHOLMOD:
      ret.insert(BaseMatrix::SPARSE_SYM);
      break;

    case BaseSolver::GINKGO:
    case BaseSolver::LIS:
      // ginkgo and lis can don't know about symmetric storage saving
      ret.insert(BaseMatrix::SPARSE_NONSYM);
      break;

    case BaseSolver::PETSC:
      ret.insert(BaseMatrix::SPARSE_SYM);
      ret.insert(BaseMatrix::SPARSE_NONSYM);
      break;


    case BaseSolver::PHIST:
       ret.insert(BaseMatrix::SPARSE_SYM);
       ret.insert(BaseMatrix::SPARSE_NONSYM);
       break;

    case BaseSolver::SUPERLU:
      ret.insert(BaseMatrix::SPARSE_NONSYM);
      break;

    case BaseSolver::LDL_SOLVER:
      ret.insert(BaseMatrix::SPARSE_SYM);
      break;

    case BaseSolver::DIAGSOLVER:
      ret.insert(BaseMatrix::DIAG);
      break;

    // The external solver should not be limited by the matrix format within cfs
    case BaseSolver::EXTERNAL_SOLVER:
      break;
  }
  return ret;
}

//! Return if preconditioner is capable of solving a SBM system
bool IsSolverSBMCapable(BaseSolver::SolverType st)
{
  switch(st) {
    // All iterative solvers can handle any sparse matrix format
    case BaseSolver::NOSOLVER:
    case BaseSolver::RICHARDSON: 
    case BaseSolver::CG:
    case BaseSolver::LANCZOS: 
    case BaseSolver::QMR:
    case BaseSolver::GMRES:
    case BaseSolver::MINRES:
    case BaseSolver::SYMMLQ: 
      return true;

    case BaseSolver::LAPACK_LU:
    case BaseSolver::LAPACK_LL:
    case BaseSolver::PARDISO_SOLVER:  
    case BaseSolver::LU_SOLVER:
    case BaseSolver::CHOLMOD: 
    case BaseSolver::LDL_SOLVER:
    case BaseSolver::DIAGSOLVER:

    // check these
    case BaseSolver::UMFPACK:
    case BaseSolver::LIS:
    case BaseSolver::GINKGO:
    case BaseSolver::PETSC:
    case BaseSolver::SUPERLU:
    case BaseSolver::PHIST:
    case BaseSolver::EXTERNAL_SOLVER:
      return false;
  }
  assert(false);
  return false; // please compiler
}


}
