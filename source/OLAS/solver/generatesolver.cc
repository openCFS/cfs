// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <def_use_ilupack.hh>
#include <def_use_lapack.hh>
#include <def_use_pardiso.hh>
#include <def_use_cholmod.hh>

#include "DataInOut/ParamHandling/ParamNode.hh"

#include "generatesolver.hh"
#include "General/exception.hh"

#ifdef USE_LAPACK
#include "OLAS/external/lapack/lapacklu.hh"
#include "OLAS/external/lapack/lapackll.hh"
#endif

#ifdef USE_PARDISO
#include "OLAS/external/pardiso/pardisosolver.hh"
#endif

#ifdef USE_ILUPACK
#include "OLAS/external/ilupack/Ilupack.hh"
#endif

#ifdef USE_CHOLMOD
#include "OLAS/external/cholmod/CholMod.hh"
#endif

// include source code for templated solvers
#include "basesolver.hh"
#include "richardson.hh"
#include "cgsolver.hh"
#include "gmres.hh"
#include "minres.hh"
#include "lusolver.hh"
#include "ldlsolver.hh"
#include "diagsolver.hh"

namespace CoupledField {



// **********************************************************
//   Macro for generation of a solver object intended for
//   the case that the solver depends not only on the entry
//   type, but also on the block size of the matrix.
// *********************************************************
#define SOLVER_OBJ( entryType, solverObjType ) \
  if ( entryType == eType  ) { \
    retSolver = new solverObjType( params, report ); \
    (*cla) << " GenerateSolver: Generated " \
    << MACRO2STRING( solverObjType ) \
    << " solver" << std::endl;\
  }



// ****************************
//   Generate a solver object
// ****************************
BaseSolver* GenerateSolverObject( const BaseMatrix &mat,
                                  SolverType solver,
                                  ParamNode* xml,
                                  InfoNode*  olasInfo,
                                  OLAS_Params *params,
                                  OLAS_Report *report ){

  BaseSolver *retSolver = NULL;
  BaseMatrix::EntryType eType = mat.GetEntryType();
  bool eTypeUnknown = false;

  ParamNode* solver_xml = NULL;
  if(xml != NULL && xml->Has("solver")) solver_xml = xml->Get("solver");
  
  // Branch depending on desired solver
  switch( solver ) {

  case RICHARDSON:
    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new RichardsonSolver<Double>( params, report );
      ASSERTMEM( retSolver, sizeof(RichardsonSolver<Double>) );
      (*cla) << " GenerateSolver: Generated real Richardson solver"
              << std::endl;
    }
    else if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new RichardsonSolver<Complex>( params, report );
      ASSERTMEM( retSolver, sizeof(RichardsonSolver<Complex>) );
      (*cla) << " GenerateSolver: Generated complex Richardson solver"
              << std::endl;
    }
    else {
      eTypeUnknown = true;
    }
    break;

  case DIAGSOLVER:
    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new DiagSolver<Double>( params, report );
      ASSERTMEM( retSolver, sizeof(DiagSolver<Double>) );
      (*cla) << " GenerateSolver: Generated real Diagonal solver"
      << std::endl;
    }
    else if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new DiagSolver<Complex>( params, report );
      ASSERTMEM( retSolver, sizeof(DiagSolver<Complex>) );
      (*cla) << " GenerateSolver: Generated complex Diagonal solver"
      << std::endl;
    }
    else {
      eTypeUnknown = true;
    }
    break;

  case CG:
    if(eType == BaseMatrix::DOUBLE) {
      retSolver = new CGSolver<Double>(solver_xml, params, report );
      (*cla) << " GenerateSolver: Generated real CG solver" << std::endl;
    }
    if(eType == BaseMatrix::COMPLEX) {
      retSolver = new CGSolver<Complex>(solver_xml, params, report );
      (*cla) << " GenerateSolver: Generated complex CG solver" << std::endl;
    }
    break;

  case GMRES:
    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new GMRESSolver<Double>( params, report );
      ASSERTMEM( retSolver, sizeof(GMRESSolver<Double>) );
      (*cla) << " GenerateSolver: Generated real GMRES solver" << std::endl;
    }
    else if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new GMRESSolver<Complex>( params, report );
      ASSERTMEM( retSolver, sizeof(GMRESSolver<Complex>) );
      (*cla) << " GenerateSolver: Generated complex GMRES solver"
      << std::endl;
    }
    else {
      eTypeUnknown = true;
    }
    break;

  case MINRES:
    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new MINRESSolver<Double>( params, report );
      ASSERTMEM( retSolver, sizeof(MINRESSolver<Double>) );
      (*cla) << " GenerateSolver: Generated real MINRES solver" << std::endl;
    }
    else if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new MINRESSolver<Complex>( params, report );
      ASSERTMEM( retSolver, sizeof(MINRESSolver<Complex>) );
      (*cla) << " GenerateSolver: Generated complex MINRES solver"
      << std::endl;
    }
    else {
      eTypeUnknown = true;
    }
    break;

  case LU_SOLVER:
    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new LUSolver<Double>( params, report );
      ASSERTMEM( retSolver, sizeof(LUSolver<Double>) );
      std::string tmp;
      Enum2String( LU_SOLVER, tmp );
      (*cla) << " GenerateSolver: Generated real "
      << tmp << std::endl;
    }
    else if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new LUSolver<Complex>( params, report );
      ASSERTMEM( retSolver, sizeof(LUSolver<Complex>) );
      std::string tmp;
      Enum2String( LU_SOLVER, tmp );
      (*cla) << " GenerateSolver: Generated complex "
      << tmp << std::endl;
    }
    else {
      eTypeUnknown = true;
    }
    break;
  case LDL_SOLVER:
    // Get block size of matrix entries (we need an StdMatrix for this,
    // but, if it's not, an LDL_Solver makes no sense anyhow)
  {
    SOLVER_OBJ( BaseMatrix::DOUBLE,   LDLSolver<Double>  );
    SOLVER_OBJ( BaseMatrix::COMPLEX,  LDLSolver<Complex> );
  }
  break;


#ifdef USE_LAPACK
  case LAPACK_LU:
    if ( mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX ) {
      EXCEPTION( "LAPACK_LU only works with a LAPACK_GBMATRIX!");
    }
    else {
      const StdMatrix &stdmat = dynamic_cast<const StdMatrix &>(mat);
      if ( stdmat.GetStorageType() != BaseMatrix::LAPACK_GBMATRIX ) {
        EXCEPTION( "LAPACK_LU only works with a LAPACK_GBMATRIX!");
      }
      else {
        retSolver = new Lapack_LU( params, report );
        ASSERTMEM( retSolver, sizeof(Lapack_LU) );
        (*cla) << " GenerateSolver: Generated Lapack_LU solver"
        << std::endl;
      }
    }
    break;

  case LAPACK_LL:
    if ( mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX ) {
      EXCEPTION( "LAPACK_LL only works with an SCRS_MATRIX!" );
    }
    else {
      const StdMatrix &stdmat = dynamic_cast<const StdMatrix &>(mat);
      if ( stdmat.GetStorageType() != BaseMatrix::SPARSE_SYM ) {
        EXCEPTION( "LAPACK_LL only works with an SCRS_MATRIX!" );
      }
      else {
        retSolver = new Lapack_LL( params, report );
        ASSERTMEM( retSolver, sizeof(Lapack_LL) );
        (*cla) << " GenerateSolver: Generated Lapack_LL solver"
        << std::endl;
      }
    }
    break;
#else
  case LAPACK_LU:
    EXCEPTION( "Compile with USE_LAPACK to enable support for LAPACK_LU solver" );
    break;
  case LAPACK_LL:
    EXCEPTION( "Compile with USE_LAPACK to enable support for LAPACK_LL solver" );
    break;
#endif

  case PARDISO:

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

    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new PardisoSolver<Double>( params, report, solver_xml );
      ASSERTMEM( retSolver, sizeof(PardisoSolver<Double>) );
      (*cla) << " GenerateSolver: Generated real Pardiso solver"
      << std::endl;
    }
    if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new PardisoSolver<Complex>( params, report, solver_xml );
      ASSERTMEM( retSolver, sizeof(PardisoSolver<Complex>) );
      (*cla) << " GenerateSolver: Generated complex Pardiso solver"
      << std::endl;
    }
#else

    EXCEPTION( "Compile with USE_PARDISO to enable interface to Pardiso "
    "library" );
#endif
    break;


  case ILUPACK_SOLVER:

#ifdef USE_ILUPACK
  {
    // Check suitability of matrix
    if (mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX)
      EXCEPTION("Ilupack only works with (S)CRS_Matrix class!");

    const StdMatrix &stdmat = dynamic_cast<const StdMatrix &>(mat);
    if(stdmat.GetStorageType() != BaseMatrix::SPARSE_NONSYM
        && stdmat.GetStorageType() != BaseMatrix::SPARSE_SYM  )
      EXCEPTION("Ilupack only works with (S)CRS_Matrix class!");

    if(eType == BaseMatrix::DOUBLE) {
      retSolver = new Ilupack<Double>(solver_xml, olasInfo, eType);
      (*cla) << " GenerateSolver: Generated real ilupack solver" << std::endl;
    }
    if(eType == BaseMatrix::COMPLEX) {
      retSolver = new Ilupack<Complex>(solver_xml, olasInfo, eType);
      (*cla) << " GenerateSolver: Generated complex ilupack solver" << std::endl;
    }

  }
#else
  EXCEPTION("Compile with USE_ILUPACK to enable interface to ilupack");
#endif
  break;
  
  case CHOLMOD:
#ifdef USE_CHOLMOD
  {
    if(mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX || dynamic_cast<const StdMatrix &>(mat).GetStorageType() != BaseMatrix::SPARSE_SYM){
      EXCEPTION("CholMod only works with SCRS_Matrix class!");
    }
    if(eType == BaseMatrix::DOUBLE){
      retSolver = new CholMod<Double>(solver_xml, olasInfo, eType);
      (*cla) << " GenerateSolver: Generated real CholMod solver" << std::endl;
    }else if(eType == BaseMatrix::COMPLEX){
      retSolver = new CholMod<Complex>(solver_xml, olasInfo, eType);
      (*cla) << " GenerateSolver: Generated complex CholMod solver" << std::endl;
    }
  }
#else
    EXCEPTION("Compile with USE_CHOLMOD to enable interface to CholMod");
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

}
