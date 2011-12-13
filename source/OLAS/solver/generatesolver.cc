// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stddef.h>
#include <map>
#include <ostream>
#include <string>
#include <utility>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/basematrix.hh"
#include "MatVec/stdmatrix.hh"
#include "def_use_cholmod.hh"
#include "def_use_ilupack.hh"
#include "def_use_lapack.hh"
#include "def_use_pardiso.hh"
#include "generatesolver.hh"

#ifdef USE_LAPACK
#include "OLAS/external/lapack/lapackll.hh"
#include "OLAS/external/lapack/lapacklu.hh"
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
#include "cgsolver.hh"
#include "diagsolver.hh"
#include "gmres.hh"
#include "ldlsolver.hh"
#include "lusolver.hh"
#include "minres.hh"
#include "richardson.hh"

namespace CoupledField {



// **********************************************************
//   Macro for generation of a solver object intended for
//   the case that the solver depends not only on the entry
//   type, but also on the block size of the matrix.
// *********************************************************
#define SOLVER_OBJ( entryType, solverObjType ) \
  if ( entryType == eType  ) { \
    retSolver = new solverObjType( solverXML, olasInfo ); \
    (*cla) << " GenerateSolver: Generated " \
    << MACRO2STRING( solverObjType ) \
    << " solver" << std::endl;\
  }



// ****************************
//   Generate a solver object
// ****************************
BaseSolver* GenerateSolverObject( const BaseMatrix &mat,
                                  PtrParamNode xml,
                                  PtrParamNode  olasInfo ){

  BaseSolver *retSolver = NULL;
  BaseMatrix::EntryType eType = mat.GetEntryType();
//  bool eTypeUnknown = false; // TODO: Unused variable eTypeUnknown
  std::string solverStr = "";

  PtrParamNode solverXML;
  solverXML = xml->Get("solver", ParamNode::INSERT);
  
  EnumMap::iterator it, end;
  it = BaseSolver::solverType.map.begin();
  end = BaseSolver::solverType.map.end();
  
  for( ; it != end; it++ ) {
    if( solverXML->HasByVal("type", it->second) ) {
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
      retSolver = new RichardsonSolver<Double>( solverXML, olasInfo );
      ASSERTMEM( retSolver, sizeof(RichardsonSolver<Double>) );
      (*cla) << " GenerateSolver: Generated real Richardson solver"
              << std::endl;
    }
    else if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new RichardsonSolver<Complex>( solverXML, olasInfo );
      ASSERTMEM( retSolver, sizeof(RichardsonSolver<Complex>) );
      (*cla) << " GenerateSolver: Generated complex Richardson solver"
              << std::endl;
    }
//    else {
//      eTypeUnknown = true;
//    }
    break;
  case BaseSolver::DIAGSOLVER:
    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new DiagSolver<Double>(solverXML, olasInfo );
      ASSERTMEM( retSolver, sizeof(DiagSolver<Double>) );
      (*cla) << " GenerateSolver: Generated real Diagonal solver"
      << std::endl;
    }
    else if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new DiagSolver<Complex>(solverXML, olasInfo );
      ASSERTMEM( retSolver, sizeof(DiagSolver<Complex>) );
      (*cla) << " GenerateSolver: Generated complex Diagonal solver"
      << std::endl;
    }
//    else {
//      eTypeUnknown = true;
//    }
    break;

  case BaseSolver::CG:
    if(eType == BaseMatrix::DOUBLE) {
      retSolver = new CGSolver<Double>(solverXML, olasInfo );
      (*cla) << " GenerateSolver: Generated real CG solver" << std::endl;
    }
    if(eType == BaseMatrix::COMPLEX) {
      retSolver = new CGSolver<Complex>(solverXML, olasInfo );
      (*cla) << " GenerateSolver: Generated complex CG solver" << std::endl;
    }
    break;

  case BaseSolver::GMRES:
    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new GMRESSolver<Double>(solverXML, olasInfo );
      ASSERTMEM( retSolver, sizeof(GMRESSolver<Double>) );
      (*cla) << " GenerateSolver: Generated real GMRES solver" << std::endl;
    }
    else if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new GMRESSolver<Complex>(solverXML, olasInfo );
      ASSERTMEM( retSolver, sizeof(GMRESSolver<Complex>) );
      (*cla) << " GenerateSolver: Generated complex GMRES solver"
      << std::endl;
    }
//    else {
//      eTypeUnknown = true;
//    }
    break;

  case BaseSolver::MINRES:
    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new MINRESSolver<Double>( solverXML, olasInfo );
      ASSERTMEM( retSolver, sizeof(MINRESSolver<Double>) );
      (*cla) << " GenerateSolver: Generated real MINRES solver" << std::endl;
    }
    else if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new MINRESSolver<Complex>( solverXML, olasInfo );
      ASSERTMEM( retSolver, sizeof(MINRESSolver<Complex>) );
      (*cla) << " GenerateSolver: Generated complex MINRES solver"
      << std::endl;
    }
//    else {
//      eTypeUnknown = true;
//    }
    break;

  case BaseSolver::LU_SOLVER:
    if ( eType == BaseMatrix::DOUBLE ) {
      retSolver = new LUSolver<Double>( solverXML, olasInfo );
      ASSERTMEM( retSolver, sizeof(LUSolver<Double>) );
      (*cla) << " GenerateSolver: Generated real "
      << BaseSolver::solverType.ToString(solver) << std::endl;
    }
    else if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new LUSolver<Complex>( solverXML, olasInfo );
      ASSERTMEM( retSolver, sizeof(LUSolver<Complex>) );
      (*cla) << " GenerateSolver: Generated complex "
      << BaseSolver::solverType.ToString(solver) << std::endl;
    }
//    else {
//      eTypeUnknown = true;
//    }
    break;
  case BaseSolver::LDL_SOLVER:
    // Get block size of matrix entries (we need an StdMatrix for this,
    // but, if it's not, an LDL_Solver makes no sense anyhow)
  {
    SOLVER_OBJ( BaseMatrix::DOUBLE,   LDLSolver<Double>  );
    SOLVER_OBJ( BaseMatrix::COMPLEX,  LDLSolver<Complex> );
  }
  break;


#ifdef USE_LAPACK
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
        retSolver = new Lapack_LU( solverXML, olasInfo );
        ASSERTMEM( retSolver, sizeof(Lapack_LU) );
        (*cla) << " GenerateSolver: Generated Lapack_LU solver"
        << std::endl;
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
        retSolver = new Lapack_LL( solverXML, olasInfo );
        ASSERTMEM( retSolver, sizeof(Lapack_LL) );
        (*cla) << " GenerateSolver: Generated Lapack_LL solver"
        << std::endl;
      }
    }
    break;
#else
  case BaseSolver::LAPACK_LU:
    EXCEPTION( "Compile with USE_LAPACK to enable support for LAPACK_LU solver" );
    break;
  case BaseSolver::LAPACK_LL:
    EXCEPTION( "Compile with USE_LAPACK to enable support for LAPACK_LL solver" );
    break;
#endif

  case BaseSolver::PARDISO:

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
      retSolver = new PardisoSolver<Double>( solverXML, olasInfo );
      ASSERTMEM( retSolver, sizeof(PardisoSolver<Double>) );
      (*cla) << " GenerateSolver: Generated real Pardiso solver"
      << std::endl;
    }
    if ( eType == BaseMatrix::COMPLEX ) {
      retSolver = new PardisoSolver<Complex>( solverXML, olasInfo );
      ASSERTMEM( retSolver, sizeof(PardisoSolver<Complex>) );
      (*cla) << " GenerateSolver: Generated complex Pardiso solver"
      << std::endl;
    }
#else

    EXCEPTION( "Compile with USE_PARDISO to enable interface to Pardiso "
    "library" );
#endif
    break;


  case BaseSolver::ILUPACK:

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
      retSolver = new Ilupack<Double>(solverXML, olasInfo, eType);
      (*cla) << " GenerateSolver: Generated real ilupack solver" << std::endl;
    }
    if(eType == BaseMatrix::COMPLEX) {
      retSolver = new Ilupack<Complex>(solverXML, olasInfo, eType);
      (*cla) << " GenerateSolver: Generated complex ilupack solver" << std::endl;
    }

  }
#else
  EXCEPTION("Compile with USE_ILUPACK to enable interface to ilupack");
#endif
  break;
  
  case BaseSolver::CHOLMOD:
#ifdef USE_CHOLMOD
  {
    if(mat.GetStructureType() != BaseMatrix::SPARSE_MATRIX || dynamic_cast<const StdMatrix &>(mat).GetStorageType() != BaseMatrix::SPARSE_SYM){
      EXCEPTION("CholMod only works with SCRS_Matrix class!");
    }
    if(eType == BaseMatrix::DOUBLE){
      retSolver = new CholMod<Double>(solverXML, olasInfo, eType);
      (*cla) << " GenerateSolver: Generated real CholMod solver" << std::endl;
    }else if(eType == BaseMatrix::COMPLEX){
      retSolver = new CholMod<Complex>(solverXML, olasInfo, eType);
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
