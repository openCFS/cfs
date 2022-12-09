#include <def_use_arpack.hh>
#include <def_use_phist_ev.hh>
#include <def_use_pardiso.hh>
#include <def_use_feast.hh>
#include <def_build_quadeigensolver.hh>
#include <def_use_superlu.hh>

#include "MatVec/BaseMatrix.hh"
#include "OLAS/algsys/SolStrategy.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "generateEigensolver.hh"

#include "BaseEigenSolver.hh"



#ifdef USE_ARPACK
  #include "OLAS/external/arpack/ArpackEigenSolver.hh"
#endif
#ifdef USE_PHIST_EV
  #include "OLAS/external/phist/PhistEigenSolver.hh"
   #include "OLAS/external/phist/PhistCore.hh"
#endif

#ifdef USE_FEAST
  #include "OLAS/external/feast/FeastEigenSolver.hh"
#endif

#if defined(USE_ARPACK) && defined(USE_SUPERLU)
  #include "PALMSolver.hh"
#endif


#include "QuadraticEigenSolver.hh"


namespace CoupledField {

DEFINE_LOG(genEigSolver, "genEigSolver")

  // *********************************
  //   Generate a EigenSolver object
  // *********************************
  BaseEigenSolver* GenerateEigenSolverObject( shared_ptr<SolStrategy> strat,
                                              PtrParamNode eSolverList,
                                              PtrParamNode solverList,
                                              PtrParamNode precondList,
                                              PtrParamNode  eigenInfo,
                                              const std::string& opt_solver_string) {
    
    BaseEigenSolver* retSolver = NULL;

    // If an optional solver is given via string create it, otherwise create the one given in strat
    std::string eSolverId="";
    if(opt_solver_string.length()){
      eSolverId = opt_solver_string;
    }
    else{
      eSolverId = strat->GetEigenSolverId();
    }
    // Obtain current eigensolver id from strategy object and try
    // to find eigenSolver in xml list.
    ParamNodeList sNodes =  eSolverList->GetChildren();
    LOG_DBG(genEigSolver) << "GESO: id=" << eSolverId << " childs=" << sNodes.GetSize();
    PtrParamNode eSolverXML;
    for( UInt i = 0; i < sNodes.GetSize(); ++i ) {
      LOG_DBG(genEigSolver) << "GESO: test " << sNodes[i]->Get("id")->As<std::string>();
      if( sNodes[i]->Get("id")->As<std::string>() == eSolverId ) {
        eSolverXML = sNodes[i]; 
      }
    }

    if(!eSolverXML)
      EXCEPTION("Solver with id '" << eSolverId << "' was not found!");

    // Convert string to enum
    std::string solverStr;
    for(auto it = BaseEigenSolver::eigenSolverType.map.begin(); it != BaseEigenSolver::eigenSolverType.map.end(); it++ )
    {
      if(eSolverXML->GetName() == it->second)
      {
        if(solverStr != "")
          EXCEPTION("Two eigensolvers have been specified: " << solverStr << " and " << (it->second))
        solverStr = it->second;
      }
    }
  
    if(solverStr == "")
      EXCEPTION("Could not determine eigensolver!")

    BaseEigenSolver::EigenSolverType solver = BaseEigenSolver::NO_EIGENSOLVER;
    solver = BaseEigenSolver::eigenSolverType.Parse(solverStr);

    LOG_DBG(genEigSolver) << "GESO: solver -> " << solverStr;

    // Failsafe if the optional solver given was QUADRATIC to avoid a recursion.
    if(opt_solver_string.length()){
      if(solver==BaseEigenSolver::QUADRATIC){
        EXCEPTION("Use a different Solver for the GEVP, e.g. FEAST");
      }
    }

    // Branch depending on desired EigenSolver
    switch(solver)
    {
    case BaseEigenSolver::ARPACK:
      #ifdef USE_ARPACK
        retSolver = new ArpackEigenSolver( strat, eSolverXML, solverList, precondList, eigenInfo );
      #else
        EXCEPTION( "compiled without ARPACK: set USE_ARPACK=ON to use the ARPACK solver" );
      #endif
      break;

    case BaseEigenSolver::PHIST:
      #ifdef USE_PHIST_EV
        retSolver = new PhistEigenSolver( strat, eSolverXML, solverList, precondList, eigenInfo );
      #else
        EXCEPTION( "compiled without PHIST_EV: set USE_PHIST_EV=ON to use the PHIST_EV solver" );
      #endif
      break;

    case BaseEigenSolver::FEAST:
      #ifdef USE_FEAST
        retSolver = new FeastEigenSolver(strat, eSolverXML, solverList, precondList, eigenInfo);
      #else
        EXCEPTION( "compiled without FEAST: set USE_FEAST=ON to use the FEAST solver!" );
      #endif
        break;

    case BaseEigenSolver::PALM:
      #if defined(USE_ARPACK) && defined(USE_SUPERLU)
        retSolver = new PALMEigenSolver(strat, eSolverXML, solverList, precondList, eigenInfo);
      #else
        EXCEPTION( "compiled without PALM: set USE_ARPACK=ON and USE_SUPERLU=ON to use the PALM solver.");
      #endif
        break;

    case BaseEigenSolver::QUADRATIC:
      #if defined(USE_ARPACK) || defined(USE_FEAST) || defined(USE_PHIST_EV)
        retSolver = new QuadraticEigenSolver(strat, eSolverXML, solverList, precondList, eigenInfo,eSolverList);
      #else
        EXCEPTION( "compiled without Quadratic Eigenvalue solver: enable an eigenvalue solver like ARPACK, FEAST or PHIST_EV" );
      #endif
        break;

    case BaseEigenSolver::NO_EIGENSOLVER:
      assert(false);
    }

    retSolver->eigenSolverType_ = solver;
    return retSolver;
  }
}
