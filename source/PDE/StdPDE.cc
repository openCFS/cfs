// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "StdPDE.hh"

#include "Utils/vector.hh"
#include "Driver/stdSolveStep.hh"
#include "Driver/transientdriver.hh"

#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"

// headers for Paramhandling
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/ParamHandling/CFSOLASParams.hh"
#include "Driver/assemble.hh"

namespace CoupledField {


  StdPDE::StdPDE(Grid *aptgrid, ParamNode* paramNode )
    : BasePDE( paramNode )
      {

    ENTER_FCN( "StdPDE::StdPDE");
    
    // =====================================================================
    // initialize variables
    // =====================================================================
    numPDENodes_ = 0;
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
  
  StdPDE::~StdPDE() {
    ENTER_FCN( "StdPDE::~StdPDE" );
  }



  const Vector<Double>& StdPDE::getS1() const {
  
    ENTER_FCN( "StdPDE::getS1" );
  
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
    
    ENTER_FCN( "StdPDE::getS2" );
    
    if ( TS_alg_ != NULL ) {
      return TS_alg_->GetDeriv2();
    }
    else {
      EXCEPTION( pdename_ << ":getS2: No timestepping defined for this PDE" );

      // Only a dummy line for compiler
      return TS_alg_->GetDeriv2();
    }
  }

  // ======================================================
  // GRID SECTION (Meshing, ...) 
  // ======================================================

  // calculates L2-norm of RHS regarding dirichlet entries due to penalty 
  // formulation by setting them 0
  Double StdPDE::RhsL2Norm(Vector<Double>& actRHS)
  {
    ENTER_FCN( "StdPDE::RhsL2Norm" );

    Integer eqnNr;
  
    // Eliminate inhom. dirichlet node from RHS (due to penalty formulation)
    for ( UInt i = 0; i < idBcs_.GetSize(); i++ ) {

      // Get grip of current idBC
      InhomDirichletBc const & actBc = *idBcs_[i];

      // Get entity iterator
      EntityIterator it = actBc.entities->GetIterator();
      
      for ( it.Begin(); !it.IsEnd(); it++ ) {
        eqnNr = eqnMap_->GetEqn( *actBc.result, it, actBc.dof );
        if ( eqnNr != 0 ) {
          actRHS[(eqnNr-1)] = 0.0;
        }
      }
    } 
    

    return actRHS.NormL2();
  }


  void StdPDE::CreateMatrices_Solver() {
    
    ENTER_FCN( "StdPDE::CreateMatrices_Solver" );
    
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
      algsys_->CreateSolver();
      SETPROFILE("Before CreatePrecond()");
      algsys_->CreatePrecond();
      
    } else {
      // create eigenvalue solver
      algsys_->CreateEigenSolver();
    }
        
    // now reset AlgebraicSystem 
    algsys_->InitRHS();
    algsys_->InitSol();
    SETPROFILE("-- Finished CreateMatrices_Solver--");
  }


  BaseSolveStep * StdPDE::GetSolveStep() {
    
    ENTER_FCN( "StdPDE::GetSolveStep" );
    return solveStep_;
  }
  

  // ======================================================
  // ALGSYS SECTION (SOLVER, ...) 
  // ======================================================

  Double StdPDE::GetFracDampMatrixCoeff(RegionIdType regionId) {
    
    ENTER_FCN( "StdPDE::GetFracDampMatrixCoeff" );

    Double coeff;

    // pre factor of fractional derivative (same for all algorithms)
    TransientDriver * driver = NULL; 
    driver = dynamic_cast<TransientDriver*>(domain->GetSingleDriver() );
    
    if( driver == NULL) {
      EXCEPTION( "Fractional damping only possible for transient simulation!" );
    }
    Double dt = driver->GetDeltaT();

    Double y;
    materials_[regionId]->GetScalar(y,FRACTIONAL_EXPONENT,REAL);

    coeff = std::exp(-(y-1.0) * std::log(dt));

    // needed for formulation with only MASS and STIFFNESS matrix
    // pre factor of Newmark time stepping scheme
    Double beta = TS_alg_->GetNewmarkBeta();
    coeff *= 1.0 / (beta*dt*dt);

#ifdef DEBUG
    (*debug) << std::endl << "Parameters in GetFracDampMatrixCoeff are:"
             << std::endl << "regionId: " << regionId
             << std::endl << "y    = " << y
             << std::endl << "dt   = " << dt
             << std::endl << "beta = " << beta
             << std::endl << "coeff= " << coeff << std::endl << std::endl;
#endif

    return coeff;
  }
  

  // real valued method (for TRANSIENT and STATIC)
  void StdPDE::GetSolVecOfElement( Vector<Double>& elemSol,
                                   const EntityIterator& it ) {

    ENTER_FCN( "StdPDE::GetSolVecOfElement" );

    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *results_[0], it );


    elemSol.Resize( eqns.GetSize() );
    elemSol.Init(0);
    NodeStoreSol<Double> * solhelp = 
      dynamic_cast<NodeStoreSol<Double>*>(sol_);
    Vector<Double> sol = solhelp->GetAlgSysVector();
    
    for( UInt i = 0; i < eqns.GetSize(); i++ ) {
      if ( eqns[i] != 0 ) {
        elemSol[i] = sol[(abs(eqns[i]-1))];
      } else {
        elemSol[i] = 0.0;
      }
     }
  }


  // complex valued method (for HARMONIC)
  void StdPDE::GetSolVecOfElement( Vector<Complex>& elemSol,
                                   const EntityIterator& it ) {

    ENTER_FCN( "StdPDE::GetSolVecOfElement" );

    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *results_[0], it );


    elemSol.Resize( eqns.GetSize() );
    elemSol.Init( Complex(0.0, 0.0) );
    NodeStoreSol<Complex> * solhelp = 
      dynamic_cast<NodeStoreSol<Complex>*>(sol_);
    Vector<Complex> sol = solhelp->GetAlgSysVector();
    
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
                                       const EntityIterator& it) {

    ENTER_FCN( "StdPDE::GetDerivSolVecOfElement" );
    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *results_[0], it );
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
                                        const EntityIterator& it) {

    ENTER_FCN( "StdPDE::GetDerivSolVecOfElement" );

    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *results_[0], it );
    
    sol.Resize( eqns.GetSize() );
    sol.Init( 0.0 );
    
    // we obtain from assemble: frequency =  2*PI*actFreq
    Double omega = solveStep_->GetActFreq() * 2 * PI;
    Complex jomega = Complex(0.0,omega);
    
    if ( analysistype_ == HARMONIC ) {
      NodeStoreSol<Complex> * solhelp = 
 	dynamic_cast<NodeStoreSol<Complex>*>(sol_);
      Vector<Complex> solAtNode = solhelp->GetAlgSysVector();

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
                                         const EntityIterator& it) {

    ENTER_FCN( "StdPDE::GetDeriv2SolVecOfElement" );

    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *results_[0], it );
    
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
                                         const EntityIterator& it) {

    ENTER_FCN( "StdPDE::GetDeriv2SolVecOfElement" );

    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *results_[0], it );
    
    sol.Resize( eqns.GetSize() );
    sol.Init( 0.0 );
    
    // we obtain from assemble: frequency 
    Double omega = solveStep_->GetActFreq() * 2 * PI;
    
    if ( analysistype_ == HARMONIC ) {
      NodeStoreSol<Complex> * solhelp = 
 	dynamic_cast<NodeStoreSol<Complex>*>(sol_);
      Vector<Complex> solAtNode = solhelp->GetAlgSysVector();

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

    ENTER_FCN( "StdPDE::StoreAlgsysToVec" );

    //const UInt numElems = numPDENodes_ * dofspernode_;
    UInt numElems = eqnMap_->GetNumEqns();

    vec.Resize(numElems);
    for (UInt i=0; i<numElems; i++) {
      vec[i] = pt[i];
    }
  }

  CFSVector *  StdPDE::GetSolutionVector() {
    ENTER_FCN( "StdPDE::GetSolutionVector" );
    return solVec_;
  }



  shared_ptr<ResultInfo> StdPDE::GetResultInfo( SolutionType solType ) {
    ENTER_FCN( "StdPDE::GetResult" );
    
    shared_ptr<ResultInfo> res;
    for( UInt i = 0; i < results_.GetSize(); i++ ) {
      if( results_[i]->resultType == solType) {
        res = results_[i];
        break;
      }
    }

    // check, if resultinfo was found
    if( !res.get() ) {
      std::string solString;
      Enum2String( solType, solString );
      EXCEPTION( "A result with solutionType '" << solString
                 << "' was not found for " << pdename_ );
    }

    return res;
  }


  // ******************
  //   ReadOlasParams
  // ******************
  void StdPDE::ReadOlasParams( std::string sysName ) {

    ENTER_FCN( "StdPDE::ReaOlasParams" );

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

} // end of namespace
