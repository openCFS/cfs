#include <fstream>
#include <iostream>
#include <string>

#include "solveStepAcousticMechBubble.hh"
#include "assemble.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Forms/forms_header.hh"
#include "General/environment.hh"
#include "PDE/timestepping.hh"
#include "PDE/SinglePDE.hh"
#include "Domain/domain.hh"

#include "CoupledPDE/DirectCoupledPDE.hh"


namespace CoupledField {

  SolveStepAcousticMechBubble::SolveStepAcousticMechBubble(DirectCoupledPDE& apde,
                                                   BubbleDynType bubbleDynType) 
    : StdSolveStep(apde) 
  {
    ENTER_FCN( "SolveStepAcousticMechBubble::SolveStepAcousticMechBubble" );
     std::cerr << ":SolveStepAcousticMechBubble" << std::endl;
  
    directPDE_ = &apde;
    
    // determine pointers to acoustic and mechanc PDE
    acouPDE_ = domain->GetSinglePDE("acoustic");
    mechPDE_ = domain->GetSinglePDE("mechanic");


    bubbleDynType_ = bubbleDynType;
  
    if( bubbleDynType_ != NOBUBBLETYPE ) {
      SetupBubbleDynamics();
    }
  

  }

  SolveStepAcousticMechBubble::~SolveStepAcousticMechBubble() {
    ENTER_FCN( "SolveStepAcousticMechBubble::~SolveStepAcousticMechBubble" );
  }
 
  // ======================================================
  // Solve Step Transient SECTION  
  // ======================================================


  void SolveStepAcousticMechBubble::SolveStepTrans( const Boolean reset ) {

    ENTER_FCN( "SolveStepAcousticMechBubble::SolveStepTrans" );

    StepTransBubble(reset);

  }


  void SolveStepAcousticMechBubble::StepTransBubble( const Boolean reset ) {

    ENTER_FCN( "SolveStepAcousticMechBubble::StepTransBubble" );

    //    StdVector<SinglePDE*>& singlePDEs = directPDE_->GetSinglePDEs();
//     UInt numPDENodes = singlePDEs[0]->getPDE_numPDENodes() +
//       singlePDEs[1]->getPDE_numPDENodes();

    UInt numeq = directPDE_->GetTotalUnknowns();


    Double * ptsol;
    UInt job;

    Double oldpressure=0;
    Double oldpressureDeriv=0;

    Double nlRhsNorm;

    Vector<Double> newSol(numeq);
    Vector<Double> oldSol(numeq);

    Vector<Double> nlRhsNew(numeq);
    Vector<Double> nlRhsOld(numeq);
    nlRhsOld.Init();

    //check for assembling of matrices
    job = 3;

    if (actStep_ == 1) {
      job = 1;
      PDE_.AssembleMatrices();
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
    }  

    //ACHTUNG: Wurde angepasst
    //perform predictor and update RHS according to time stepping algorithm
    //NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
    //Vector<Double> & sol4pred= solhelp->GetAlgSysVector();

 
   Vector<Double> & sol4pred = 
      dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
    TS_alg_->Predictor(sol4pred);


 
    //set old solution  
    oldSol = sol4pred;


    //do nonlinear iteration
    UInt iterationCounter=0;
    Boolean performOneMoreStep;




     do {

    //  std::cerr<<"schleife "<<  iterationCounter<<std::endl;



      iterationCounter++;
      // for every time step write out number of iteration loops to standard out
      //    if (iterationCounter == 1)
      // std::cout << std::endl << "Time step:   "  << kstep 
      //  << "  ,Iterations: " << iterationCounter;
      //else 
      // std::cout        << "  " << iterationCounter;
    
#ifdef DEBUG
      *debug << std::endl
             << "====================================================== "
             << std::endl
             <<   "Nonlinear Acoustics: Perform internal loop no. "
             << iterationCounter << std::endl;      
#endif
    
      // Set RHS to zero
      algsys_->InitRHS();
    
      //compute source term due to bubble ODE
        ComputeBubbleRHS();



      algsys_->GetRHSVal( ptsol );
      for (UInt i=0; i<nlRhsNew.GetSize(); i++) {
        nlRhsNew[i] = ptsol[i];
      }
     
      //      StoreAlgsysToVec(nlRhsNew, ptsol );       


      //compute norm
      Double normRhs=0.0;
      nlRhsNorm = 0.0;
      for (UInt i=0; i<nlRhsNew.GetSize(); i++) {
        nlRhsNorm += (nlRhsNew[i]-nlRhsOld[i])*(nlRhsNew[i]-nlRhsOld[i]);
        normRhs   += nlRhsNew[i]*nlRhsNew[i];
      }
     
      if (normRhs >0) {
        nlRhsNorm = sqrt(nlRhsNorm) / sqrt(normRhs);
      }
      else {
        nlRhsNorm = sqrt(nlRhsNorm);
      }
 
      nlRhsOld = nlRhsNew;

      //account for source terms RHS
      PDE_.AssembleSrcRHS(actTime_);
    
    
      //acoust for time stepping
      TS_alg_->UpdateRHS();
    


      SetBCs(actTime_);

      if ( iterationCounter>1 ) 
        job = 3;

      algsys_->BuildInDirichlet();
      if ( job == 1 ) {
        algsys_->SetupSolver();
        algsys_->SetupPrecond();
      }
        
      //solve the PDE system
      algsys_->Solve();
    
      //ACHTUNG: Wurde angepasst
      //get the solution
      UInt length = algsys_->GetSolutionVal(ptsol);
      PDE_.SaveSolution(ptsol,length);
      
      // store new solution in newSol
      //StoreAlgsysToVec(newSol, ptsol );
      Vector<Double> & solHelp = 
        dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
       TS_alg_->Corrector(solHelp);

       newSol = solHelp;
      // perform corrector step, if effective mass formulation is used,
      //   we need the Corrector step before we store newsol to sol_,
      //   because newsol is second time derivative at first!
      //TS_alg_->Corrector(newSol);

          newSol = solHelp;
      // compute L2-Norm of error between last incremental solution and
      //   actual incremental solution
      Double solIncrL2Norm=0;
      for (UInt i=0; i<newSol.GetSize(); i++){
        solIncrL2Norm += (newSol[i]-oldSol[i])*(newSol[i]-oldSol[i]);
      }
      solIncrL2Norm = sqrt(solIncrL2Norm);
      Double actSolL2Norm = newSol.NormL2();
    
      Double incrementalErr;
      if (actSolL2Norm != 0.0)
        incrementalErr = solIncrL2Norm / actSolL2Norm;
      else
        incrementalErr = solIncrL2Norm;



      // output of norms and data
      if ( nonLinLogging_ == TRUE ) {
	Info->WriteNonLinIter(pdename_, iterationCounter, nlRhsNorm,
			      incrementalErr );
      }
    


      //store new solution to old one;
      oldSol = newSol;
    


      // ***********************************************************
      //   Compute bubble ODE
      // ***********************************************************


      Matrix<Double> ptCoord;
      BaseFE         * ptElem;
      StdVector<UInt> connect, Eqns;  
      UInt numEl=0;
      
      StdVector<RegionIdType> acouRegions = *(acouPDE_->getSDsPDE());

      BaseNodeStoreSol * acouSol = acouPDE_->getPDESolution();

      for (UInt actSD=0; actSD<acouRegions.GetSize(); actSD++) {
        StdVector<Elem*> elemssd;
        ptgrid_->GetVolElems(elemssd,acouRegions[actSD]);
      

	for (UInt j=0; j < elemssd.GetSize(); j++) {
	
	  ptElem  = elemssd[j]->ptElem;
	  connect = elemssd[j]->connect;

	  Vector<Double> elPressure;
	  Vector<Double> elPressureDeriv;
	
	  //we now have also in coupled computation a pressure formulation =============
	
	  //get the pressure
	  acouSol->GetElemSolution(elPressure,connect);

	
	  //get 1st derivative of pressure
	  acouPDE_->GetDerivSolVecOfElement(elPressureDeriv,connect);
	

	  //compute average values
	  Double pressure=0;
	  Double pressureDeriv=0;
	  for (UInt elnode=0; elnode<elPressure.GetSize(); elnode++) {
	    pressure      += elPressure[elnode];
	    pressureDeriv += elPressureDeriv[elnode];
	  }
	
	  pressure      /= elPressure.GetSize();
	  pressureDeriv /= elPressureDeriv.GetSize();


	
	  //set to ODE



	  //Dimensionless case ++++++++++++++++++++++++++++++++++++++++++++++++


	  pressureNoDim = pressure / pStatic_ ;
	  presDerivNoDim = pressureDeriv * initRadius_ / pStatic_ / ( sqrt ( pStatic_ / density_));
	  ptBubble_[numEl]->SetP(pressureNoDim);
	  ptBubble_[numEl]->SetDpdt(presDerivNoDim);

	  // In case of explicit Euler watch out suggested stepsize
	  Double dt = TS_alg_->GetTimeStep();
	  Double steptime = actTime_ - dt;

	  if ( hTry_ > dt){
	    hTry_ = dt / 3.0;
	  }

	  steptime   = steptime / initRadius_ * (sqrt(pStatic_/ density_));
	  tNoDim_    = actTime_ / initRadius_ * (sqrt(pStatic_/ density_));
	  hTry_      = hTry_ / initRadius_ * (sqrt(pStatic_/ density_));

	  //get the current values
	  bubbleValues_[0] = radius_[numEl];
	  bubbleValues_[1] = velocity_[numEl];
	  yNoDim_[0] = bubbleValues_[0] / initRadius_;
	  yNoDim_[1] = bubbleValues_[1] / (sqrt( pStatic_/ density_));

	  //set numEl to ODE-Solver
	  ptODESolver_->SetNumEl(numEl);
	
	  ptODESolver_->Solve(steptime, tNoDim_, yNoDim_, *ptBubble_[numEl], hTry_,
			      0, dt);
	
	  //set the new values
	  radiusWork_[numEl]   = yNoDim_[0] * initRadius_;
	  velocityWork_[numEl] = yNoDim_[1] * sqrt( pStatic_/ density_);
	  //Dimensionless case ++++++++++++++++++++++++++++++++++++++++++++++++++++



	  //Normal case++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	  //          ptBubble_[numEl]->SetP(pressure);
	  //          ptBubble_[numEl]->SetDpdt(pressureDeriv);
        
          // In case of explicit Euler watch out suggested stepsize
	  //          Double dt = TS_alg_->GetTimeStep();
	  //          Double steptime = actTime_ - dt;

	  //	  if ( hTry_ > dt){
	  //	    hTry_ = dt / 3.0;
	  //	  }
        
          //get the current values
	  //          bubbleValues_[0] = radius_[numEl];
	  //          bubbleValues_[1] = velocity_[numEl];
               
        
          //set numEl to ODE-Solver
	  //          ptODESolver_->SetNumEl(numEl);
      	  //          ptODESolver_->Solve(steptime, actTime_, bubbleValues_, *ptBubble_[numEl], hTry_,
	  //                              0, dt);
        
          //set the new values
	  //          radiusWork_[numEl]   = bubbleValues_[0]; 
	  //          velocityWork_[numEl] = bubbleValues_[1]; 
	  //Normal case+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          numEl++;
        }  
      }

      //check if we have reached the accuracy!
      performOneMoreStep =    (incrementalErr > incStopCrit_);// || (nlRhsNorm > incStopCrit_) ;



       } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);



    radius_   = radiusWork_;
    velocity_ = velocityWork_;

  }

  // ======================================================================
  //   Setup structure for bubbly dynamics
  // ======================================================================
  void SolveStepAcousticMechBubble::SetupBubbleDynamics(){  

    ENTER_FCN( "SolveStepAcousticMechBubble::SetupBubbleDynamics" );

    UInt numPDEElems = acouPDE_->getPDE_numElems();

    ptBubble_.Resize(numPDEElems);
    bubbleValues_.Resize(2);
    radius_.Resize(numPDEElems);  
    velocity_.Resize(numPDEElems);  
    velocity_.Resize(numPDEElems);
    radiusWork_.Resize(numPDEElems);  
    velocityWork_.Resize(numPDEElems);

    hTry_=1e-10;
    yNoDim_.Resize(2);
    // =============================================
    //  Query ParamHandler for material parameters
    // =============================================

    // First make sure that there is a section 'bubbleDynamic'
    // We rely on params using a validating Schema parser here.
    StdVector<Double> auxVec;
    params->GetList( "initRadius", auxVec, "bubbleDynamic" );
    if ( auxVec.GetSize() != 1 ) {
      Error( "Cannot find initRadius! Assuming that section 'bubbleDynamic' "
             "is missing in xml-file", __FILE__, __LINE__ );
    }

    //Double initRadius = auxVec[0];
    initRadius_ = auxVec[0];
    //    Double initVel = 0.0;
    params->Get( "initVel", initVel_, "bubbleDynamic" );
    //    Double density = 0.0;
    params->Get( "density", density_, "bubbleDynamic" );
    //    Double sonicVel = 0.0;
    params->Get( "sonicVel", sonicVel_, "bubbleDynamic" );
    //    Double pStatic = 0.0;
    params->Get( "pStatic", pStatic_, "bubbleDynamic" );
    //    Double pVapour = 0.0;
    params->Get( "pVapour", pVapour_, "bubbleDynamic" );
    //    Double surfaceTension = 0.0;
    params->Get( "surfaceTension", surfaceTension_, "bubbleDynamic" );
    //    Double polytrop = 0.0;
    params->Get( "polytrop", polytrop_, "bubbleDynamic" );
    //    Double viscosity = 0.0;
    params->Get( "viscosity", viscosity_, "bubbleDynamic" );

    //Set bubbledensity
    params->Get( "bubbleNumDensity", bubbleDensity_, "acoustic", "bubbles" );


    Info->PrintF ("", "The following paramters and methods were used\n");
    Info->PrintF ("", "to compute the bubble behaviour:\n\n");
    Info->PrintF ("", "Initial Radius of the bubbles:\t %10.6e\n", initRadius_);
    Info->PrintF ("", "Initial Velocity of the bubble wall:\t %10.6e\n", initVel_);
    Info->PrintF ("", "Density of bubbles per unit volume:\t %10.6e\n", bubbleDensity_);
    Info->PrintF ("", "Density of the fluid:\t %10.6e\n", density_);
    Info->PrintF ("", "Sound velocity of the fluid:\t %10.6e\n", sonicVel_);
    Info->PrintF ("", "Static pressure:\t %10.6e\n", pStatic_);
    Info->PrintF ("", "Vapour pressure of the fluid:\t %10.6e\n", pVapour_);
    Info->PrintF ("", "Surface tension of the fluid:\t %10.6e\n", surfaceTension_);
    Info->PrintF ("", "Polytropic exponent of the fluid:\t %10.6e\n", polytrop_);
    Info->PrintF ("", "Viscosity of the fluid:\t %10.6e\n\n", viscosity_);

    Info->PrintF ("", "\nAnfangsschrittweite:\t %10.6e\n", hTry_);



    //Set the nonlinear parameters for the iteration between bubble and pressure
      // perform logging?
      nonLinLogging_ = params->IsSet( "logging", "acoustic", "nonLinear" );
      // incremental stopping criterion
      params->Get("incStopCrit", incStopCrit_, "acoustic", "nonLinear" );
      // residual stopping criterion
      params->Get("resStopCrit", residualStopCrit_, "acoustic", "nonLinear" );
      // maximal number of NL-iterations
      params->Get("maxNumIters", nonLinMaxIter_, "acoustic", "nonLinear");

    Info->PrintF ("", "Incrementelle Abbruchbedingung:\t %10.6e\n" , incStopCrit_);
    Info->PrintF ("", "Maximale Anzahl von Iterationen:\t %10.6e\n", nonLinMaxIter_);

    for (UInt el=0; el<numPDEElems; el++) {
      // Choice which bubbledynamical method is used and 
      // creation of pointer to choosen class
      // so far method is chosen in cfs program call with -bubbletyp
      switch(bubbleDynType_){
      case KELLERMIKSIS:
	ptBubble_[el] = new KellerMiksis(initRadius_,density_, sonicVel_, pStatic_, 
					 pVapour_, surfaceTension_, polytrop_, viscosity_);
	if (el == 0){
	  Info->PrintF( pdename_, "Using Keller-Miksis-Bubble-Model\n");
	}
        break;
     

      case GILMORE:
	//ptBubble_[el] = new Gilmore(initRadius_,density_, sonicVel_, pStatic_, 
	//			pVapour_, surfaceTension_, polytrop_, viscosity_);
// 	if (el == 0){
// 	  Info->PrintF( pdename_, "Using Gilmore-Bubble-Model\n");
// 	}
	ptBubble_[el] = new Gilmoredimlos(initRadius_,density_, sonicVel_, pStatic_, 
				      pVapour_, surfaceTension_, polytrop_, viscosity_);
	if (el == 0){
	  Info->PrintF( pdename_, "Using dimensionless Gilmore-Bubble-Model\n");
	}
	break;

      default:
	Error("No bubblemethod specified ",__FILE__,__LINE__);
      }
  
      radius_[el]       = initRadius_;
      velocity_[el]     = initVel_;
      radiusWork_[el]   = initRadius_;
      velocityWork_[el] = initVel_;
   
    }

    // Generate ODE solver object
    ptODESolver_ = new ODESolver_RKF45;
    Info->PrintF( pdename_, "Using the Runge-Kutta-Method to solve bubbledynamics\n");
    // ptODESolver_ = new ODESolver_Rosenbrock;
    // Info->PrintF( pdename_, "Using the Rosenbrock-Method to solve bubbledynamics\n");
    // Info->PrintF( pdename_, "which computation of the jacobian\n");
    // Info->PrintF( pdename_, "which approximation of the jacobian\n");



  }


  void SolveStepAcousticMechBubble::ComputeBubbleRHS(){  
    ENTER_FCN( "SolveStepAcousticMechBubble::ComputeBubbleRHS" );

    NodeEQN * ptEQN = acouPDE_->getPDE_eqnData();
    StdVector<RegionIdType> *subdoms = acouPDE_->getSDsPDE();
    PdeIdType pdeId = acouPDE_->GetPDEId();

    //rhs-integrator
    Double dummy=1.0;
    BaseForm *rhsForm = new VolumeSrcInt(dummy, isaxi_);        

    UInt numEl = 0;
    Double fluidDensity;


    for (UInt actDom=0; actDom <  subdoms->GetSize(); actDom++) {      
      StdVector<Elem*> elemssd;
      ptgrid_->GetVolElems(elemssd, (*subdoms)[actDom]);
    

      for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++)    {              
        BaseFE * ptEl = elemssd[actEl]->ptElem;
        StdVector<UInt> connecth = elemssd[actEl]->connect;

        Matrix<Double> ptCoord;
        GetElemCoords(connecth, ptCoord);
        
        Double beta2 = 1;

	bubbleValues_[0] = radiusWork_[numEl];
	bubbleValues_[1] = velocityWork_[numEl];
	
	// New rhs for cavitational fluid
	// rho0 * 4/3 * pi* n * ( 3 * R^2 * d^2R/dt^2 + 6 * R * (dR/dt)^2) 
	if (actStep_ == 1) 
	  beta2 =density_ * density_   *  4.0/3.0*PI*bubbleDensity_
	    *6*bubbleValues_[0]*bubbleValues_[1]*bubbleValues_[1]; 
	else {
	  StdVector<Double> dydt(2);
	  dydt[0]=0;
	  dydt[1]=0;

	  // dimensionless**************************************
	  yNoDim_[0] = bubbleValues_[0] / initRadius_;
	  yNoDim_[1] = bubbleValues_[1] / (sqrt( pStatic_/ density_));
          tNoDim_    = actTime_ / initRadius_ * (sqrt(pStatic_/ density_));
	  ptBubble_[numEl]->CompDeriv(tNoDim_, yNoDim_, dydt);
	  dydt[1] = dydt[1] * pStatic_ / ( density_ * initRadius_ );
	  beta2 =density_ * density_   *  4.0/3.0*PI*bubbleDensity_*
	    (6*bubbleValues_[0]*bubbleValues_[1]*bubbleValues_[1]
	     + 3*bubbleValues_[0]*bubbleValues_[0]*dydt[1] ); 
	  // dimensionless**************************************

	  //Normal case+++++++++++++++++++++++++++++++++++++++
	  //	  ptBubble_[numEl]->CompDeriv(actTime_, bubbleValues_, dydt);
	  //	  beta2 = density_ * density_   * 4.0/3.0*PI*bubbleDensity_*
	  //	    (6*bubbleValues_[0]*bubbleValues_[1]*bubbleValues_[1]
	  //	     + 3*bubbleValues_[0]*bubbleValues_[0]*dydt[1] ); 
	  //Normal case+++++++++++++++++++++++++++++++++++++++
	}


	rhsForm->SetFactor(beta2);            
	rhsForm->SetElemPtr(ptEl);
      
	Vector<Double> elemVec;
	rhsForm->CalcElemVector(ptCoord, elemVec);

	StdVector<Integer> connect_PDE;
	ptEQN->Node2EQN(connecth, connect_PDE);
	algsys_->SetElementRHS(&elemVec[0], pdeId, connect_PDE.GetPointer(), 
			       connect_PDE.GetSize());

	numEl++;
      }
    }

    delete rhsForm;
  

  }


  //  return pointer to vector with first derivative of solution
  const StdVector<Double>& SolveStepAcousticMechBubble::GetResultData(std::string resultType) {

    ENTER_FCN( "SolveStepAcousticMechBubble::WriteResults" );

    if ( resultType == "bubbleRadius" ) {
      return radius_;
    }
    else if ( resultType == "bubbleVelocity" ) {
      return velocity_;
    }
    else {
      Error("SolveStepAcousticBubble: Result-Type not possible",
            __FILE__, __LINE__ );
      
      // only to satisfy compiler
      return velocity_;
    }

  }

  

} // end of namespace

