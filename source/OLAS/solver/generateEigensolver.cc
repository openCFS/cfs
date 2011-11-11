#include <def_use_arpack.hh>

#include "MatVec/BaseMatrix.hh"

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
                                              PtrParamNode xml,
                                              PtrParamNode  eigenInfo ) {
    
    BaseEigenSolver *retSolver = NULL;
    // TODO: Check if this is still needed
    // BaseMatrix::EntryType eType = mat.GetEntryType();
    // COMPWARNING: unused variable bool eTypeUnknown = false;
    
    std::string solverStr;
    PtrParamNode solverXML;
    solverXML = xml->Get("eigenSolver", ParamNode::INSERT);
  
    EnumMap::iterator it, end;
    it = BaseEigenSolver::eigenSolverType.map.begin();
    end = BaseEigenSolver::eigenSolverType.map.end();
  
    for( ; it != end; it++ ) {
      if( solverXML->HasByVal("type", it->second) ) {
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
      retSolver = new ArpackEigenSolver( xml, eigenInfo, eigenInfo );
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
