#include <fstream>
#include <iostream>
#include <string>

#include "stdSolveStep.hh"
#include "assemble.hh"

#include "Utils/nodestoresol.hh"
#include "PDE/StdPDE.hh"

namespace CoupledField {

  StdSolveStep::StdSolveStep(StdPDE & apde) 
    :BaseSolveStep(),
     PDE_(apde)
  {

    ENTER_FCN( "StdSolveStep::StdSolveStep" );
 
    pdename_      = PDE_.pdename_;
    numPDENodes_  = PDE_.numPDENodes_;
    numPDEElems_  = PDE_.numElems_;
    isaxi_        = PDE_.isaxi_;
    subdoms_      = PDE_.subdoms_;
    materialData_ = PDE_.materials_;
    ptgrid_       = PDE_.ptgrid_;
    algsys_       = PDE_.algsys_;
    sol_          = PDE_.sol_;
    eqnData_      = PDE_.eqnData_;
    assemble_     = PDE_.assemble_;
    outFile_      = PDE_.outFile_;

    TS_alg_        = PDE_.TS_alg_;

    lineSearch_       = PDE_.lineSearch_;
    nonLin_           = PDE_.nonLin_;
    isHyst_           = PDE_.isHysteresis_;
    incStopCrit_      = PDE_.incStopCrit_;
    residualStopCrit_ = PDE_.residualStopCrit_;
    nonLinMaxIter_    = PDE_.nonLinMaxIter_;
    nonLinMethod_     = PDE_.nonLinMethod_;
    nonLinLogging_    = PDE_.nonLinLogging_ ;
    nonLinPDEName_    = PDE_.nonLinPDEName_;

    isIterCoupled_       = PDE_.isIterCoupled_;
    firstTimeStepStatic_ =PDE_.firstTimeStepStatic_;
    iterCoupledCounter_ = &PDE_.iterCoupledCounter_;

    isIncrFormulation_  = PDE_.isIncrFormulation_;


    pdeId1_   = NO_PDE_ID;
    pdeId2_   = NO_PDE_ID;
    numReset_ = 0;

    startStep_ = 1;
  }

  
    //! Destructor
  StdSolveStep::~StdSolveStep() {
    ENTER_FCN( "StdSolveStep::~StdSolveStep" );

  }

  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

  // time is used for a series of static calculations
  // don't get confused with REAL transient simulations!
  void StdSolveStep::SolveStepStatic( const Boolean reset ) {

    ENTER_FCN( "StdSolveStep::SolveStepStatic" );


    if (PDE_.nonLin_) {
      StepStaticNonLin(reset);
    }
    else {
      StepStaticLin(reset);
    }
  }


  void StdSolveStep::StepStaticLin( const Boolean reset ) {

    ENTER_FCN( "StdSolveStep::StepStaticLin" );

    Integer job = 3; // only update BCs
    Double * ptSol;
    UInt size = 0;

    // If the geometry has changed or the system matrix
    // is calculated for the first time,
    // the matrices have to be reassembled and therfore
    // the preconditioner has to be recalculated
    
    if ( PDE_.geoUpdate_ == TRUE ||  PDE_.firstTimeStepStatic_ == TRUE) {
      PDE_.AssembleMatrices();
      PDE_.AssembleSprings( actTime_ );
      job = 1; // calc new preconditioner
    }
  
    // The RHS-sources have to be reassembled each time
    PDE_.AssembleSrcRHS( actTime_ );

    PDE_.SetBCs( actTime_ );

    PDE_.algsys_->ConstructEffectiveMatrix( matrix_factor_ );

    // Incorporate Boundary conitions and
    // recalc the prconditioner eventually
    PDE_.algsys_->BuildInDirichlet();

    SETPROFILE("Before SetupPrecond");
    PDE_.algsys_->SetupPrecond();
    SETPROFILE("After SetupPrecond / Before SetupSolver");

    PDE_.algsys_->SetupSolver( );
    SETPROFILE("After SetupSolver / Before Solve");

    // Solve problem
    PDE_.algsys_->Solve();
    SETPROFILE("After Solve");
    

    // Get the solution and store it
    size = PDE_.algsys_->GetSolutionVal(ptSol);
    PDE_.SaveSolution(ptSol,size);

    PDE_.firstTimeStepStatic_ = FALSE;
  }


  // ======================================================
  // Solve Step Transient SECTION  
  // ======================================================

  void StdSolveStep::PreStepTrans( const Boolean reset ) {

    ENTER_FCN( "StdSolveStep::PreStepTrans" );

    // due to coupling-pdes, the RHS has to be initialized BEFORE 
    // the coupling forces are assembled to the RHS
    PDE_.algsys_->InitRHS();

    if (PDE_.geoUpdate_) {
      PDE_.algsys_->InitRHS();
      PDE_.algsys_->InitSol();
      PDE_.algsys_->InitMatrix();
      PDE_.SetReassemble();  
    }
  }


  void StdSolveStep::SolveStepTrans( const Boolean reset ) {

    ENTER_FCN( "StdSolveStep::SolveStepTrans" );

    // /DELETE
    // if (PDE_.laststepcalc_ == kstep && kstep != 1) {
    //       PDE_.recalc_ = TRUE;
    //     }
    //     else {
    //       PDE_.laststepcalc_= kstep;
    //     }

    if (PDE_.nonLin_) {
      StepTransNonLin(reset);
    }
    else {
      StepTransLin(reset);
    }
  }


  //! \todo delete job parameter or replace it by a 
  //! meaningfull attribute
  void StdSolveStep::StepTransLin( const Boolean reset ) {

    ENTER_FCN( "StdSolveStep::StepTransLin" );

    Double * ptsol;
    Integer job;
    UInt length = 0;

    //account for RHS
    PDE_.AssembleSrcRHS( actTime_ );

    if ( PDE_.isIterCoupled_ == FALSE || PDE_.iterCoupledCounter_ == 0 ) {        
      Vector<Double> & solHelp = 
        dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
      PDE_.TS_alg_->Predictor(solHelp);
    }

    if ( actStep_ == 1 || actStep_ == startStep_) {
      job = 3;
      // why is the first statement checking for 'pdeIsCoupled'?
      if ( PDE_.isIterCoupled_ == FALSE || PDE_.iterCoupledCounter_ == 0 
           || PDE_.geoUpdate_ == TRUE ) {
        job = 1;
        PDE_.AssembleMatrices();
        PDE_.AssembleSprings( actTime_ );
        PDE_.algsys_->ConstructEffectiveMatrix(matrix_factor_);
      }  
    }
    else if (reset) {
      job = 1;
	  
      PDE_.algsys_->InitMatrix();
      PDE_.AssembleSprings( actTime_ );
      PDE_.algsys_->ConstructEffectiveMatrix(matrix_factor_);
    }
    else {
      job = 3;

      // The following section is only an experiment up to now
      if ( PDE_.geoUpdate_ == TRUE && PDE_.isIterCoupled_ == TRUE ) {
        if (PDE_.isIncrFormulation_) {
          Error( "Incremental formulation and geoUpdate are currently not "
                 "working together", __FILE__, __LINE__ );
        }
        job = 1;
        PDE_.AssembleMatrices();
        PDE_.AssembleSprings( actTime_ );
        PDE_.algsys_->ConstructEffectiveMatrix(matrix_factor_);      
      }
    }

    if (PDE_.isIncrFormulation_) {
      Vector<Double> & solvector =
        dynamic_cast<NodeStoreSol<Double>*>(PDE_.sol_)->GetAlgSysVector();

      solvector *= -(1.0);
      PDE_.algsys_->UpdateRHS(SYSTEM,solvector.GetPointer());
    }

    PDE_.TS_alg_->UpdateRHS();

    PDE_.SetBCs( actTime_ );
    PDE_.algsys_->BuildInDirichlet();

    if ( job == 1 ) {
      SETPROFILE("Before SetupPrecond");
      PDE_.algsys_->SetupPrecond( );
      SETPROFILE("After SetupPrecond / Before SetupSolver");
      PDE_.algsys_->SetupSolver( );
      SETPROFILE("After SetupSolver / Before Solve");  
    }

    PDE_.algsys_->Solve();
    SETPROFILE("After Solve");

    if ( PDE_.isIncrFormulation_ ) {

      //what a heuristic!!!!!
      Double relaxVal = 0.5;

      if (PDE_.iterCoupledCounter_ == 0)
        relaxVal = 0.1;
      if (PDE_.iterCoupledCounter_ == 1)
        relaxVal = 0.25;
      
      //we do not want to have iterations!!
      if (PDE_.incStopCrit_ > 1.0 ) {
        relaxVal = 1.0;
        if (PDE_.iterCoupledCounter_ == 1) {
          (*error) << "Mechanic-Acoustic-Coupling with Norm > 1 should have "
                   << "just a forward coupling from Mechanic to Acoustic";
          Error( __FILE__, __LINE__ );
        }
      }

      length = algsys_->GetSolutionVal(ptsol);
      solIncr_.Replace(length, ptsol, FALSE);

      if (PDE_.iterCoupledCounter_ == 0) {
        actSol_ = solIncr_ * relaxVal;
      } else {
        actSol_ += solIncr_ * relaxVal;
      }

      PDE_.SaveSolution(actSol_.GetPointer(), length);
    }
    else {
      length = algsys_->GetSolutionVal(ptsol);
      PDE_.SaveSolution(ptsol,length);
    }


    Vector<Double> & solHelp = 
      dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());

    //std::cerr << "Doing Corrector step of PDE " << PDE_.GetName() << std::endl;
    PDE_.TS_alg_->Corrector(solHelp);

     if ( PDE_.isIterCoupled_ ) {
      PDE_.iterCoupledCounter_++;
    }


  }


  void StdSolveStep::SolveStepTrans4Slice( const Boolean reset ) {

    ENTER_FCN( "StdSolveStep::SolveStepTrans4Slice" );

    Double * ptsol;
    UInt length = 0;

    if ( reset ) {
      numReset_ += 1;
    }

    //account for RHS
    PDE_.AssembleSrcRHS( actTime_ );

    //perform predictor step
    Vector<Double> & solHelp1 = 
      dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
    PDE_.TS_alg_->Predictor(solHelp1);

    if ( actStep_ == 1) {
      PDE_.AssembleMatrices();
      PDE_.AssembleSprings( actTime_ );
      PDE_.algsys_->ConstructEffectiveMatrix(matrix_factor_);
    }
    else if (reset) {
      PDE_.algsys_->InitMatrix(SYSTEM);
      PDE_.algsys_->ConstructEffectiveMatrix(matrix_factor_);
    }

    //update RHS due to time stepping
    PDE_.TS_alg_->UpdateRHS();

    if ( numReset_ == 0) {
      PDE_.SetBCs( actTime_ );
      PDE_.algsys_->BuildInDirichlet();
    }

    if ( actStep_ == 1 || reset == TRUE ) {
      PDE_.algsys_->SetupPrecond( );
      PDE_.algsys_->SetupSolver( );
    }

    PDE_.algsys_->Solve();

    //get the solution
    length = algsys_->GetSolutionVal(ptsol);
    PDE_.SaveSolution(ptsol,length);
   
    //perform corrector step 
    Vector<Double> & solHelp2 = 
      dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
    PDE_.TS_alg_->Corrector(solHelp2);
        
  }


  void StdSolveStep::PostStepTrans( ) {

    ENTER_FCN( "StdSolveStep::PostStepTrans" );

//     std::string type_str;
//     Enum2String( PDE_.GetSolutionVector()->GetEntryType(), type_str);
//     std::cout << " In StdSolveStep::PostStepTrans the solution vector of"
//               << pdename_ << "computation is of type: "
//               << type_str << std::endl << std::endl;

    if ( PDE_.GetFracDamping() ) {
      Vector<Double> & solHelp = 
        dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());

      // Following method is essential for fractional damping model
      PDE_.TS_alg_->AdvanceTimestep(solHelp);
    }
  
  }



  // ======================================================
  // Solve Step Harmonic  SECTION  
  // ======================================================w

  void StdSolveStep::PreStepHarmonic( const Boolean reset ) {

    ENTER_FCN( "StdSolveStep::PreStepHarmonic" );

  
    PDE_.SetFrequency( actFreq_ );
    PDE_.algsys_->InitRHS();
    if (reset) {
      PDE_.algsys_->InitMatrix();
      PDE_.SetReassemble();
    }
  }


  void StdSolveStep::SolveStepHarmonic(const Boolean reset ) {

    ENTER_FCN( "StdSolveStep::SolveStepHarmonic" );

    if ( PDE_.nonLin_ ) {
      StepHarmonicNonLin( reset );
    }
    else {
      StepHarmonicLin( reset );
    }
  }


  void StdSolveStep::StepHarmonicLin( const Boolean reset ) {

    ENTER_FCN( "StdSolveStep::StepHarmonicLin" );

    Integer job;
    Complex * ptSol = NULL;
    UInt length = 0;

    if ( reset ) {
      PDE_.AssembleMatrices( );
      PDE_.AssembleSpecial( );
    }


    //this has to be done each time!
    PDE_.AssembleSrcRHS( actFreq_ );


    if (PDE_.ComputeRHSforHarm_)
      {
        // Evaluating RHS with nodal srcs for harmonic flownoise problems
        PDE_.ComputeRHS(actFreq_);
      }
    
    
    if ( reset ) {

      //account for bcs
      PDE_.SetBCs( actFreq_ );

      // calc new preconditioner
      job = 1;
    }
    else {
      job = 3;
    }

    PDE_.algsys_->ConstructEffectiveMatrix( matrix_factor_ );
    PDE_.algsys_->BuildInDirichlet();
 
    if (job == 1) {
      PDE_.algsys_->SetupPrecond();
      PDE_.algsys_->SetupSolver();
    }

    PDE_.algsys_->Solve();

    length =  PDE_.algsys_->GetSolutionVal(ptSol);
    PDE_.SaveSolution(ptSol,length);
  }

  // ======================================================
  // METHODS FOR EIGENVALUE COMPUTATION
  // ======================================================
  
  UInt StdSolveStep::CalcEigenFrequencies( Vector<Double> & frequencies,
                                           Vector<Double> & errBounds,
                                           UInt numFreq, Double shift ) {
    ENTER_FCN( "StdSolveStep::CalcEigenFrequencies" );
    
    // If the geometry has changed or the system matrix
    // is calculated for the first time,
    // the matrices have to be reassembled and therfore
    // the preconditioner has to be recalculated
    PDE_.algsys_->InitRHS();
    PDE_.algsys_->InitSol();
    PDE_.algsys_->InitMatrix();
    
    PDE_.AssembleMatrices();

    // Setup solver
    PDE_.algsys_->SetupEigenSolver( numFreq, shift, false );

    // Calculate eigenfrequencies
    const Double * val = NULL;
    const Double * err = NULL;
    UInt numConverged = PDE_.algsys_->CalcEigenFrequencies( val, err);

    // Copy eigenvalues into vector
    frequencies.Resize( numConverged );
    errBounds.Resize( numConverged );
    for ( UInt i = 0; i < numConverged; i++ ) {
      frequencies[i] = val[i];
      errBounds[i] = err[i];
    }

    return numConverged;
  }

  UInt StdSolveStep::CalcEigenFrequencies( Vector<Complex> & frequencies,
                                           Vector<Double> & errBounds,
                                           UInt numFreq, Double shift ) {
    ENTER_FCN( "StdSolveStep::CalcEigenFrequencies<Complex>" );
    
    // If the geometry has changed or the system matrix
    // is calculated for the first time,
    // the matrices have to be reassembled and therfore
    // the preconditioner has to be recalculated
    PDE_.algsys_->InitRHS();
    PDE_.algsys_->InitSol();
    PDE_.algsys_->InitMatrix();
    
    PDE_.AssembleMatrices();
    
    // Setup solver
    PDE_.algsys_->SetupEigenSolver( numFreq, shift, true );
    
    // Calculate eigenfrequencies
    const  Complex* val = NULL;
    const Double * err = NULL;
    UInt numConverged = PDE_.algsys_->CalcEigenFrequencies( val, err);
    
    // Copy eigenvalues into vector
    frequencies.Resize( numConverged );
    errBounds.Resize( numConverged );
    for ( UInt i = 0; i < numConverged; i++ ) {
      frequencies[i] = val[i];
      errBounds[i] = err[i];
    }
    
    return numConverged;
  }

  void StdSolveStep::CalcEigenMode( UInt numMode ) {
    
    ENTER_FCN( "StdSolveStep::CalcEigenMode" );

    Double * ptSol = NULL;
    UInt size = 0;
    
    algsys_->CalcEigenMode( numMode );

    // Get the solution and store it
    size = PDE_.algsys_->GetSolutionVal(ptSol);
    PDE_.SaveSolution(ptSol,size);
  }


  

  // ======================================================
  // METHODS FOR NONLINEAR ANALYSIS
  // ======================================================

  // sets excitation coil and returns L2Norm of them
  Double StdSolveStep::SetLinRHS( Double loadFactor)
  {
    ENTER_FCN( "StdSolveStep::SetLinRHS" );

    Double RhsLinL2Norm;  


    // to incorporate loads
    assemble_->AssembleSrcRHS(actTime_); 

    // Stores rhs vector into extForces and returns that L2-norm
 
    Double *solPtr;
    algsys_->GetRHSVal( solPtr );
    StoreAlgsysToVec(RhsLinVal_, solPtr );

    RhsLinVal_ *= loadFactor;

    RhsLinL2Norm = RhsLinVal_.NormL2();
 
    // If extForcesL2Norm is 0, no residual norm can be calculated
    if (!RhsLinL2Norm) {
      Warning("Zero external force vector!! ", __FILE__,__LINE__);
    }

    return RhsLinL2Norm;
  }


  //stores an algsys_ vector into a StdVector and returns that L2-norm
  void StdSolveStep::StoreAlgsysToVec(Vector<Double>& vec, Double * pt) {

    ENTER_FCN( "StdSolveStep::StoreAlgsysToVec" );

    //const UInt numElems = numPDENodes_ * dofspernode_;
    UInt numElems = eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN();
    vec.Resize(numElems);
      
    for (UInt i=0; i<numElems; i++) {
      vec[i] = pt[i];
    }
  }


  void StdSolveStep::SetTimeStep( Double dt ) {
    ENTER_FCN( "StdSolveStep::SetTimeStep") ;

    // Check if timestepping is eneabled
    if ( TS_alg_ != NULL )
      TS_alg_->Init( matrix_factor_, dt );
  }
  

  Double StdSolveStep::LineSearch(Vector<Double>& solIncrement, Vector<Double>& actSol, 
                                  Double& etaLineSearch, Boolean trans)
  {
    ENTER_FCN( "StdSolveStep::LineSearch" );

    Vector<Double> solOld(actSol);
    const UInt nrEtas = 5;
    const Double eta[nrEtas] = {1, 0.5, 0.25, 0.125, 0.1};
    Double etaOpt;
    Double residualL2NormOpt = 1e15;

    for(UInt i=0; i<nrEtas; i++)
      {
        actSol = solIncrement * eta[i];
        actSol += solOld;

        //store new solution
        sol_->SetAlgSysVector(actSol);

        // recalculate RHS with new values to get new residual (f^(k+1))========
        algsys_->InitRHS(RhsLinVal_.GetPointer());

        if(trans) {
          assemble_->AssembleNLRHS( actTime_ );
          TS_alg_->UpdateRHS(actSol);
        }
        else {
          assemble_->AssembleNLRHS();
        }


        // =====================================================================
        // calculation of error norms
        // =====================================================================
        Vector<Double> actRHS;
        Double *rhsPtr;
        algsys_->GetRHSVal( rhsPtr );
        StoreAlgsysToVec(actRHS, rhsPtr );

        // calculation of residual error =======================================
        Double residualL2Norm = RhsL2Norm(actRHS); // L2Norm of  ( f_i^(k+1) - f_a )

        //        std::cout << "LineSearch: eta=" << eta[i] << "  res=" << residualL2Norm << std::endl;

        if (residualL2Norm < residualL2NormOpt)
          {
            residualL2NormOpt = residualL2Norm;
            etaOpt = eta[i];
          }
      }

    etaLineSearch = etaOpt;
  
    actSol =  solIncrement * etaOpt;
    actSol += solOld;

    return residualL2NormOpt;
  }


  // returns that L2-norm of an algsys vector
  Double StdSolveStep::AlgsysL2Norm(Double * pt)
  {
    ENTER_FCN( "StdSolveStep::AlgsysL2Norm" );

    //const UInt numElems = numPDENodes_ * dofspernode_;
    UInt numElems = eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN();
    Double quadSum = 0;
  
    for (UInt i=0; i<numElems; i++)   
      quadSum += pt[i]*pt[i];

    return sqrt(quadSum);
  }
  

  void StdSolveStep::WriteClaNlNorms(const UInt iterationCounter, 
                                     const Double residualL2Norm, 
                                     const Double extForcesL2Norm,
                                     const Double residualErr, 
                                     const Double solIncrL2Norm, 
                                     const Double actSolL2Norm, 
                                     const Double incrementalErr)
  {
    ENTER_FCN( "StdSolveStep::WriteClaNlNorms" );
  
    *cla << std::endl << " ======================================================= " 
         << std::endl
         << " NONLINEAR ITERATION " << iterationCounter << std::endl
         << " ======================================================= " << std::endl;
    *cla << " === Residual norm:          " << residualL2Norm << std::endl;
    *cla << "     Norm of ext. forces:    " << extForcesL2Norm << std::endl;
    *cla << "     Residual error          " << residualErr << std::endl;
  
    *cla << " === Incremental sol L2Norm: " << solIncrL2Norm << std::endl;
    *cla << "     Actual solution L2Norm: " << actSolL2Norm << std::endl;
    *cla << "     Incremental error       " << incrementalErr << std::endl;      
  }



  // ======================================================
  // get StdPDE-methods  
  // ======================================================


  void StdSolveStep::SetBCs(const Double time )
  {
    PDE_.SetBCs(time);
  }

  void StdSolveStep::GetElemCoords(const StdVector< UInt > connect,
                                   Matrix< Double > &coordMat)
  {
    PDE_.GetElemCoords(connect, coordMat);
  }

  void StdSolveStep::GetSolVecOfElement(Vector<Double>& sol, 
                                        StdVector<UInt>& connect_PDE)
  {
    PDE_.GetSolVecOfElement(sol, connect_PDE);
  }
 
  void StdSolveStep::GetDerivSolVecOfElement(Vector<Double>& sol, 
                                             StdVector<UInt>& connect_PDE)
  {
    PDE_.GetDerivSolVecOfElement(sol, connect_PDE);
  }

  void StdSolveStep::GetDeriv2SolVecOfElement(Vector<Double>& sol, 
                                              StdVector<UInt>& connect_PDE)
  {
    PDE_.GetDeriv2SolVecOfElement(sol, connect_PDE);
  }

  Double StdSolveStep::RhsL2Norm(Vector<Double>& actRHS) {
    return PDE_.RhsL2Norm(actRHS);
  }

  void StdSolveStep::SetPDEId( const PdeIdType pdeId ) {
    ENTER_FCN( "StdSolveStep::SetPDEId" );
    
    pdeId1_ = pdeId;
  }

  void StdSolveStep::TransformSol4Slice(UInt & shiftFactor, UInt & nodeShift,
                                        UInt & elemgrid, Double &  meshsize, const UInt flag) {
    ENTER_FCN( "StdSolveStep::TransformSol4Slice" );
    
    PDE_.TransformSol4Slice(shiftFactor, nodeShift,
                            elemgrid, meshsize, flag); 
  }

  void StdSolveStep::SaveNodes(const UInt shiftFactor, const Double timeStep,
                               const UInt numShift, const Integer nodeShift, 
                               const UInt maxnumelemz_) {

    ENTER_FCN( "StdSolveStep::SaveNodes" );
    
    PDE_.SaveNodes(shiftFactor, timeStep, numShift, nodeShift, maxnumelemz_);
  }


} // end of namespace

