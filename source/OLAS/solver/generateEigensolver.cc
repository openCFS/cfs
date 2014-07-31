#include <def_use_arpack.hh>

#include "MatVec/BaseMatrix.hh"
#include "OLAS/algsys/SolStrategy.hh"

#include "generateEigensolver.hh"
#include "BaseEigenSolver.hh"

#ifdef USE_ARPACK
#include "OLAS/external/arpack/ArpackEigenSolver.hh"
#endif


namespace CoupledField {


  // *********************************
  //   Generate a EigenSolver object
  // *********************************
  BaseEigenSolver* GenerateEigenSolverObject( BaseMatrix &mat, 
                                              shared_ptr<SolStrategy> strat,
                                              PtrParamNode eSolverList,
                                              PtrParamNode solverList,
                                              PtrParamNode precondList,
                                              PtrParamNode  eigenInfo ) {
    
    BaseEigenSolver *retSolver = NULL;

    // Obtain current eigensolver id from strategy object and try
    // to find eigenSolver in xml list.
    std::string eSolverId = strat->GetEigenSolverId();
    ParamNodeList sNodes =  eSolverList->GetChildren();
    PtrParamNode eSolverXML;
    for( UInt i = 0; i < sNodes.GetSize(); ++i ) {
      if( sNodes[i]->Get("id")->As<std::string>() == eSolverId ) {
        eSolverXML = sNodes[i]; 
      }
    }

    if(! eSolverXML) {
      EXCEPTION("Solver with id '" << eSolverId << "' was not found!");
    }
  
      // Convert string to enum
    EnumMap::iterator it, end;
    it = BaseEigenSolver::eigenSolverType.map.begin();
    end = BaseEigenSolver::eigenSolverType.map.end();
  
    std::string solverStr;
    for( ; it != end; it++ ) {
      if( eSolverXML->GetName() ==  it->second) {
        if(solverStr != "")
          EXCEPTION("Two eigensolvers have been specified: " << solverStr
                    << " and " << (it->second))
        
        solverStr = it->second;
      }
    }
  
    if(solverStr == "") {
      EXCEPTION("Could not determine eigensolver!")
    }

    BaseEigenSolver::EigenSolverType solver = BaseEigenSolver::NOEIGENSOLVER;
    solver = BaseEigenSolver::eigenSolverType.Parse(solverStr);
    
    // Branch depending on desired EigenSolver
    switch( solver ) {

#ifdef USE_ARPACK
    case BaseEigenSolver::ARPACK:
      retSolver = new ArpackEigenSolver( strat, eSolverXML, solverList, precondList, eigenInfo );
      (*cla) << " GenerateEigenSolver: Generated ARPACK Eigensolver"
             << std::endl;
      break;
#endif

    case BaseEigenSolver::SUBSPACE:
      EXCEPTION( "GenerateEigenSolver: Subsapce algorithm not yet supported.\n" );
      break;

      default:
        EXCEPTION( "GenerateEigenSolver: Request for unknown solver type!" );
    }

    return retSolver;
  }

}
