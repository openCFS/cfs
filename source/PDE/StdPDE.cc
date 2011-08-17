// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "StdPDE.hh"

#include "MatVec/vector.hh"
#include "Driver/stdSolveStep.hh"
#include "Driver/transientdriver.hh"

#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"

// headers for Paramhandling
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"
#include "DataInOut/ParamHandling/CFSOLASParams.hh"
#include "Driver/assemble.hh"

#include "OLAS/algsys/basesystem.hh"

namespace CoupledField {


  StdPDE::StdPDE(Grid *aptgrid, ParamNode* paramNode )
    : BasePDE( paramNode )
      {

    
    // =====================================================================
    // initialize variables
    // =====================================================================
    numPdeEquations_ = 0;
    numElems_ = 0;
    nonLin_ = false;
    nonLinMaterial_ = false;
    isHysteresis_ = false;
    TS_alg_ = NULL;
    effectiveMass_ = false;
    diagMass_ = false;
    firstTimeStepStatic_ = true;
    isAlwaysStatic_ = false;
    isaxi_ = false;
    isComplex_ = false;
    fracDamping_ = false;
    sol_ = NULL;     
    isIncrFormulation_ = false;
    ComputeRHSforHarm_ = false;
    solveStep_ = NULL;
    needSolPrev_ = false;
    solPrev_ = NULL;
    isSetInitialCondition_ = false;
    InitialCondition_=0.0;
    isInstationary_ = false;
    olasInfo_ = NULL; // set in the child-nodes
    // =====================================================================
    // set file pointers
    // =====================================================================
    ptgrid_     = aptgrid;
    assemble_   = NULL;
    algsys_     = NULL;
    
    // =====================================================================
    // set analysis parameters
    // =====================================================================
    couplingBCsCounter_ = 0;
    isIterCoupled_ = false;
    updateCouplingBCs_ = false;
    dim_ = ptgrid_->GetDim();
    iterCoupledCounter_ = 0;
    effectiveMass_ = false;

    // =====================================================================
    // various parameters
    // =====================================================================
    needsAlgsys_ = true;

    // check if we have an axi-symmetric setup
    isaxi_ = param->Get("domain")->Get("geometryType")->AsString() == "axi";


  }
  
  StdPDE::~StdPDE() 
  {
    
  }



  const Vector<Double>& StdPDE::getS1() const {
  
  
    if ( TS_alg_ != NULL ) {
      return TS_alg_->GetDeriv1();
    }
    else {
      EXCEPTION( pdename_ << ":getS1: No timestepping defined for this PDE" );

      // Only a dummy line for compiler
      return TS_alg_->GetDeriv1();      
    }
  }
  
  
  const Vector<Double>& StdPDE::getS2() const {
    
    
    if ( TS_alg_ != NULL ) {
      return TS_alg_->GetDeriv2();
    }
    else {
      EXCEPTION( pdename_ << ":getS2: No timestepping defined for this PDE" );

      // Only a dummy line for compiler
      return TS_alg_->GetDeriv2();
    }
  }

  const Vector<Double>& StdPDE::getOld1() const {

    if ( TS_alg_ != NULL ) {
      return TS_alg_->GetOld1();
    }
    else {
      EXCEPTION( pdename_ << ":getOld1: No timestepping defined for this PDE");

      // Only a dummy line for compiler
      return TS_alg_->GetOld1();
    }
  }
  
  // ======================================================
  // GRID SECTION (Meshing, ...) 
  // ======================================================

  // calculates L2-norm of RHS regarding dirichlet entries due to penalty 
  // formulation by setting them 0
  Double StdPDE::RhsL2Norm(Vector<Double>& actRHS)
  {

    StdVector<Integer> eqns;
    shared_ptr<BaseFeFunction> actFunction; 
    // Eliminate inhom. dirichlet node from RHS (due to penalty formulation)
    for ( UInt i = 0; i < idBcs_.GetSize(); i++ ) {

      // Get grip of current idBC
      InhomDirichletBc const & actBc = *idBcs_[i];

      // Get entity iterator
      EntityIterator it = actBc.entities->GetIterator();
      
      for ( it.Begin(); !it.IsEnd(); it++ ) {
        actFunction = GetFeFunction(actBc.result->resultType);
        actFunction->GetFeSpace()->GetEqns( eqns, it, actBc.dof );
        for(UInt iEqn = 0 ; iEqn < eqns.GetSize();iEqn){
          if ( eqns[iEqn] != 0 ) {
            actRHS[eqns[iEqn]] = 0.0;
          }
        }
      }
    } 
    

    return actRHS.NormL2();
  }


  void StdPDE::CreateMatrices_Solver() {
    
    
    // ==============================
    // create and initialize matrices 
    // ==============================

    // create algebraic system and intialize matrices
    SETPROFILE("Before CreatLinSys()");
    algsys_->CreateLinSys();
    SETPROFILE("After CreatLinSys()");
    algsys_->InitMatrix();
    
    // Check for analysistype
    if ( analysistype_ != EIGENFREQUENCY ) {
      
      // create solver and preconditioner
      SETPROFILE("Before CreateSolver()");
      algsys_->CreateSolver(olasInfo_);
      SETPROFILE("Before CreatePrecond()");
      algsys_->CreatePrecond();
      
    } else {
      // create eigenvalue solver
      algsys_->CreateEigenSolver(olasInfo_->Get("eigenSolver"));
    }
        
    // now reset AlgebraicSystem 
    algsys_->InitRHS();
  	algsys_->InitSol();

    
    
    
    SETPROFILE("-- Finished CreateMatrices_Solver--");
  }


  BaseSolveStep * StdPDE::GetSolveStep() {
    
    return solveStep_;
  }
  

  // ======================================================
  // ALGSYS SECTION (SOLVER, ...) 
  // ======================================================

  ParamNode* StdPDE::FindLinearSystem(const std::string& sysName) {


    ParamNode* pn = NULL;
    pn = param->Get("sequenceStep", "index", GenStr(sequenceStep_), false );
    if(pn != NULL) pn = pn->Get("linearSystems", false);
    if(pn != NULL) pn = pn->Get("system", "name", sysName, false);
    
    return pn;
  }

  Double StdPDE::GetFracDampMatrixCoeff(RegionIdType regionId) {
    

    Double coeff;

    // pre factor of fractional derivative (same for all algorithms)
    TransientDriver * driver = NULL; 
    driver = dynamic_cast<TransientDriver*>(domain->GetSingleDriver() );
    
    if( driver == NULL) {
      EXCEPTION( "Fractional damping only possible for transient simulation!" );
    }
    Double dt = driver->GetDeltaT();

    Double y;
    materials_[regionId]->GetScalar(y,FRACTIONAL_EXPONENT,Global::REAL);

    coeff = std::exp(-(y-1.0) * std::log(dt));

    // needed for formulation with only MASS and STIFFNESS matrix
    // pre factor of Newmark time stepping scheme
    Double beta = TS_alg_->GetNewmarkBeta();
    coeff *= 1.0 / (beta*dt*dt);

    /*
    (*debug) << std::endl << "Parameters in GetFracDampMatrixCoeff are:"
             << std::endl << "regionId: " << regionId
             << std::endl << "y    = " << y
             << std::endl << "dt   = " << dt
             << std::endl << "beta = " << beta
             << std::endl << "coeff= " << coeff << std::endl << std::endl;
    */ 
    return coeff;
  }
  

  // real valued method (for TRANSIENT and STATIC)
  void StdPDE::GetSolVecOfElement( Vector<Double>& elemSol,
                                   const EntityIterator& it,
                                   shared_ptr<ResultInfo> res ) {


    StdVector<Integer> eqns;
    shared_ptr<BaseFeFunction> aFct = GetFeFunction(res->resultType);
    aFct->GetFeSpace()->GetEqns( eqns, it );


    elemSol.Resize( eqns.GetSize() );
    elemSol.Init(0);
    NodeStoreSol<Double> * solhelp = 
      dynamic_cast<NodeStoreSol<Double>*>(sol_);
    Vector<Double> & sol = solhelp->GetAlgSysVector();
    
    for( UInt i = 0; i < eqns.GetSize(); i++ ) {
      if ( eqns[i] != 0 ) {
        elemSol[i] = sol[abs(eqns[i])-1];
      } else {
        elemSol[i] = 0.0;
      }
     }
  }


  // complex valued method (for HARMONIC)
  void StdPDE::GetSolVecOfElement( Vector<Complex>& elemSol,
                                   const EntityIterator& it,
                                   shared_ptr<ResultInfo> res ) {


    StdVector<Integer> eqns;
    shared_ptr<BaseFeFunction> aFct = GetFeFunction(res->resultType);
    aFct->GetFeSpace()->GetEqns( eqns, it );


    elemSol.Resize( eqns.GetSize() );
    elemSol.Init( Complex(0.0, 0.0) );
    NodeStoreSol<Complex> * solhelp = 
      dynamic_cast<NodeStoreSol<Complex>*>(sol_);
    Vector<Complex> & sol = solhelp->GetAlgSysVector();
    
    for( UInt i = 0; i < eqns.GetSize(); i++ ) {
      if ( eqns[i] != 0 ) {
        elemSol[i] = sol[abs(eqns[i])-1];
      } else {
        elemSol[i] = Complex(0.0, 0.0);
      }
     }
  }

  
  // real valued method (for TRANSIENT)
  void StdPDE::GetDerivSolVecOfElement(Vector<Double>& sol,
                                       const EntityIterator& it,
                                       shared_ptr<ResultInfo> res) {

    StdVector<Integer> eqns;
    shared_ptr<BaseFeFunction> aFct = GetFeFunction(res->resultType);
    aFct->GetFeSpace()->GetEqns( eqns, it );
    sol.Resize( eqns.GetSize() );
    sol.Init( 0.0 ); 
    
    if (  analysistype_ == TRANSIENT) {
      const Vector<Double> & sol_der1 = getS1();
        
      for( UInt i = 0; i < eqns.GetSize(); i++ ) {
        if ( eqns[i] != 0 ) {
          sol[i] = sol_der1[abs(eqns[i])-1];
        } else {
          sol[i] = 0.0;
        }
      }
    }
  }



  // complex valued method (for HARMONIC)
  void StdPDE::GetDerivSolVecOfElement(Vector<Complex>& sol,
                                       const EntityIterator& it,
                                       shared_ptr<ResultInfo> res) {


    StdVector<Integer> eqns;
    shared_ptr<BaseFeFunction> aFct = GetFeFunction(res->resultType);
    aFct->GetFeSpace()->GetEqns( eqns, it );
    
    sol.Resize( eqns.GetSize() );
    sol.Init( 0.0 );
    
    // we obtain from assemble: frequency =  2*PI*actFreq
    Double omega = solveStep_->GetActFreq() * 2 * PI;
    Complex jomega = Complex(0.0,omega);
    
    if ( analysistype_ == HARMONIC ) {
      NodeStoreSol<Complex> * solhelp = 
 	dynamic_cast<NodeStoreSol<Complex>*>(sol_);
      const Vector<Complex> & solAtNode = solhelp->GetAlgSysVector();

      for( UInt i = 0; i < eqns.GetSize(); i++ ) {
        if ( eqns[i] != 0 ) {
          sol[i] = jomega * solAtNode[abs(eqns[i])-1];
        } else {
          sol[i] = Complex(0.0, 0.0);
        }
      }
    }
  }



  // real valued method (for TRANSIENT)
  void StdPDE::GetDeriv2SolVecOfElement( Vector<Double>& sol,
                                         const EntityIterator& it,
                                         shared_ptr<ResultInfo> res) {


    StdVector<Integer> eqns;
    shared_ptr<BaseFeFunction> aFct = GetFeFunction(res->resultType);
    aFct->GetFeSpace()->GetEqns( eqns, it );
    
    sol.Resize( eqns.GetSize() );
    sol.Init( 0.0 );
    
    if ( analysistype_ == TRANSIENT ) {
      const Vector<Double> & sol_der2 = getS2();
      for( UInt i = 0; i < eqns.GetSize(); i++ ) {
        if ( eqns[i] != 0 ) {
          sol[i] = sol_der2[abs(eqns[i])-1];
        } else {
          sol[i] = 0.0;
        }
      }
    }
  }


  // real valued method (for HARMONIC)
  void StdPDE::GetDeriv2SolVecOfElement( Vector<Complex>& sol,
                                         const EntityIterator& it,
                                         shared_ptr<ResultInfo> res) {


    StdVector<Integer> eqns;
    shared_ptr<BaseFeFunction> aFct = GetFeFunction(res->resultType);
    aFct->GetFeSpace()->GetEqns( eqns, it );
    
    sol.Resize( eqns.GetSize() );
    sol.Init( 0.0 );
    
    // we obtain from assemble: frequency 
    Double omega = solveStep_->GetActFreq() * 2 * PI;
    
    if ( analysistype_ == HARMONIC ) {
      NodeStoreSol<Complex> * solhelp = 
 	dynamic_cast<NodeStoreSol<Complex>*>(sol_);
      Vector<Complex> & solAtNode = solhelp->GetAlgSysVector();

      for( UInt i = 0; i < eqns.GetSize(); i++ ) {
        if ( eqns[i] != 0 ) {
          sol[i] = - omega * omega *solAtNode[abs(eqns[i])-1];
        } else {
          sol[i] = Complex(0.0, 0.0);
        }
      }
    }
  }

  //stores an algsys_ vector into an StdVector
  void StdPDE::StoreAlgsysToVec(Vector<Double>& vec, Double * pt) {


    const UInt numElems = numPdeUnknowns_;

    vec.Resize(numElems);
    for (UInt i=0; i<numElems; i++) {
      vec[i] = pt[i];
    }
  }

  SingleVector *  StdPDE::GetPrevSolutionVector() {
    return solVecPrev_;
  }

  shared_ptr<ResultInfo> StdPDE::GetResultInfo( SolutionType solType ) {
    
    shared_ptr<ResultInfo> res;
    for( UInt i = 0; i < results_.GetSize(); i++ ) {
      if( results_[i]->resultType == solType) {
        res = results_[i];
        break;
      }
    }

    // check, if resultinfo was found
    if( !res.get() ) {
      EXCEPTION( "A result with solutionType '" << SolutionTypeEnum.ToString(solType)
                 << "' was not found for " << pdename_ );
    }

    return res;
  }


  // ******************
  //   ReadOlasParams
  // ******************
  void StdPDE::ReadOlasParams( std::string sysName ) {


    // Log to .las file
    (*cla) <<  " --- CFS: Setting parameters for linear system '"
           << sysName << " ---" << std::endl;

    // Set parameters for OLAS
    std::string amExpert = "no";
    param->Get( "override", amExpert, false);

    ParamNode * linSysNode = NULL;
    ParamNode * temp = param->Get("sequenceStep", "index", GenStr(sequenceStep_), false );
    if ( temp )
      temp = temp->Get("linearSystems", false);
    if ( temp ) {
      linSysNode = temp ->Get("system", "name", sysName, false );
    }
    
    CFSOLASParams::SetParams( sysName, linSysNode, olasParams_, 
                              analysistype_, assemble_,
                              (amExpert == "yes") );

  }

  //============================================================================================
  //FeFunction Methods
  //============================================================================================

  shared_ptr<BaseFeFunction> StdPDE::GetFeFunction( SolutionType solType ) {
    
    //TODO> We need to find a more failsafe way to store the entity names associated with a 
    //FeFunction
    shared_ptr<BaseFeFunction> feFct;
    SolutionType mySolType;
    if( feFunctions_.find(solType) == feFunctions_.end()){
      //ok so it could be that we are looking for a postProc Result
      if(postProcResults_.find(solType) != postProcResults_.end()){
        mySolType = postProcResults_[solType];
      }else{
          EXCEPTION( "A FeFunction descriptor with solutionType '" << SolutionTypeEnum.ToString(solType)
                      << "' was not found for " << pdename_ );
      }
    }else{
      mySolType = solType;
    }
    if(feFunctions_.find(mySolType) != feFunctions_.end()){
      feFct = feFunctions_[mySolType];
    }else{
      EXCEPTION("StdPDE::GetFeFunction: Could not find the corresponding FeFunction for the given entity name\n \
                         Did you specify all Regions, Surfregions and NamedNodes in the xml?");
    }
    return feFct;

  }

} // end of namespace
