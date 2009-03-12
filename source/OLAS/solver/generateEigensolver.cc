// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <def_use_arpack.hh>

#include "MatVec/basematrix.hh"

#include "generateEigensolver.hh"
#include "baseEigensolver.hh"

#ifdef USE_ARPACK
#include "OLAS/external/arpack/arpackEigensolver.hh"
#endif


namespace CoupledField {


  // *********************************
  //   Generate a EigenSolver object
  // *********************************
  BaseEigenSolver* GenerateEigenSolverObject( BaseMatrix &mat, 
                                              EigenSolverType solver,
                                              ParamNode* xml,
                                              OLAS_Params *params, 
                                              OLAS_Report *report ) {
    
    BaseEigenSolver *retSolver = NULL;
    // TODO: Check if this is still needed
    // BaseMatrix::EntryType eType = mat.GetEntryType();
    // COMPWARNING: unused variable bool eTypeUnknown = false;
    
    // Branch depending on desired EigenSolver
    switch( solver ) {

#ifdef USE_ARPACK
    case ARPACK:
      retSolver = new ArpackEigenSolver( xml, params, report );
      (*cla) << " GenerateEigenSolver: Generated ARPACK Eigensolver"
             << std::endl;
      break;
#endif

    case SUBSPACE:
      EXCEPTION( "GenerateEigenSolver: Subsapce algorithm not yet supported.\n" );
      break;

      default:
        EXCEPTION( "GenerateEigenSolver: Request for unknown solver type!" );
    }

    return retSolver;
  }

}
