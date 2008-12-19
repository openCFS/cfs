// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <def_use_arpack.hh>

#include "solver/generateEigensolver.hh"
#include "solver/baseEigensolver.hh"
#include "matvec/basematrix.hh"

#ifdef USE_ARPACK
#include "external/arpack/arpackEigensolver.cc"
#endif


namespace OLAS {


  // *********************************
  //   Generate a EigenSolver object
  // *********************************
  BaseEigenSolver* GenerateEigenSolverObject( BaseMatrix &mat, 
                                              EigenSolverType solver,
                                              ParamNode* xml,
                                              OLAS_Params *params, 
                                              OLAS_Report *report ) {
    
    BaseEigenSolver *retSolver = NULL;
    MatrixEntryType eType = mat.GetEntryType();
    // COMPWARNING: unused variable bool eTypeUnknown = false;
    
    // Branch depending on desired EigenSolver
    switch( solver ) {

#ifdef USE_ARPACK
    case ARPACK:
      retSolver = New ArpackEigenSolver( xml, params, report );
      (*cla) << " GenerateEigenSolver: Generated ARPACK Eigensolver"
             << std::endl;
      break;
#endif

    case SUBSPACE:
      (*error) << " GenerateEigenSolver: Subsapce algorithm not yet supported."
               << std::endl;
      Error( __FILE__, __LINE__ );
      break;

      default:
      Error( "GenerateEigenSolver: Request for unknown solver type!", __FILE__,
	     __LINE__ );
    }


    // Force instantiation of member functions of templated matrix classes
    // Note that the InstantiatePublicMethods() method itself should never
    // actually be called, but the compiler must not notice this. Thus, the
    // akward if condition.
    bool neverTrue = (eType == NOENTRYTYPE);
    if ( neverTrue ) {
      retSolver->InstantiatePublicMethods( mat );
    }

    return retSolver;
  }

}
