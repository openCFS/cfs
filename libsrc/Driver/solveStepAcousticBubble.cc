#include <fstream>
#include <iostream>
#include <string>

#include "solveStepAcousticBubble.hh"
#include "assemble.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Forms/forms_header.hh"
#include "General/environment.hh"
#include "PDE/timestepping.hh"
#include "PDE/StdPDE.hh"

namespace CoupledField {

  SolveStepAcousticBubble::SolveStepAcousticBubble(StdPDE& apde,
                                                   BubbleDynType bubbleDynType) 
    : StdSolveStep(apde) 
  {
    ENTER_FCN( "SolveStepAcousticBubble::SolveStepAcousticBubble" );
    std::cerr << ":SolveStepAcousticBubble" << std::endl;
  

    bubbleDynType_ = bubbleDynType;
  
    if( bubbleDynType_ != NOBUBBLETYPE ) {
      SetupBubbleDynamics();
    }
  

  }

  SolveStepAcousticBubble::~SolveStepAcousticBubble() {
    ENTER_FCN( "SolveStepAcousticBubble::~SolveStepAcousticBubble" );
  }
 
  // ======================================================
  // Solve Step Transient SECTION  
  // ======================================================


  void SolveStepAcousticBubble::SolveStepTrans( const Boolean reset ) {

    ENTER_FCN( "SolveStepAcousticBubble::SolveStepTrans" );

    StepTransBubble(reset);

  }


  void SolveStepAcousticBubble::StepTransBubble( const Boolean reset ) {

    ENTER_FCN( "SolveStepAcousticBubble::StepTransBubble" );

    Double * ptsol;
    UInt job;


    Double oldpressure=0;
    Double oldpressureDeriv=0;

    Double nlRhsNorm;

    Vector<Double> newSol(numPDENodes_);
    Vector<Double> oldSol(numPDENodes_);

    Vector<Double> nlRhsNew(numPDENodes_);
    Vector<Double> nlRhsOld(numPDENodes_);
    nlRhsOld.Init();

    //check for assembling of matrices
    job = 3;

    if (actStep_ == 1) {
      job = 1;
      assemble_->AssembleMatrices();
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
      StoreAlgsysToVec(nlRhsNew, ptsol );       

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
      assemble_->AssembleSrcRHS(actTime_);
    
    
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
      StoreAlgsysToVec(newSol, ptsol );
    
      // perform corrector step, if effective mass formulation is used,
      //   we need the Corrector step before we store newsol to sol_,
      //   because newsol is second time derivative at first!
      TS_alg_->Corrector(newSol);
   
      // compute L2-Norm of error between last incremental solution and
      //   actual incremental solution
      Double solIncrL2Norm=0;
      for (UInt i=0; i<newSol.GetSize(); i++)
        solIncrL2Norm += (newSol[i]-oldSol[i])*(newSol[i]-oldSol[i]);
    
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
    
      //compute bubble ODE
    
      Matrix<Double> ptCoord;
      BaseFE         * ptElem;
      StdVector<UInt> connect, Eqns;  
      UInt numEl=0;
    
      for (UInt actSD=0; actSD<subdoms_.GetSize(); actSD++) {
        StdVector<Elem*> elemssd;
        ptgrid_->GetVolElems(elemssd,subdoms_[actSD]);
      

	for (UInt j=0; j < elemssd.GetSize(); j++) {
	
	  ptElem  = elemssd[j]->ptElem;
	  connect = elemssd[j]->connect;

		


	  Vector<Double> elPressure;
	  Vector<Double> elPressureDeriv;
	
	  //we now have also in coupled computation a pressure formulation =============
	
	  //get the pressure
	  sol_->GetElemSolution(elPressure,connect);

	
	  //get 1st derivative of pressure
	  GetDerivSolVecOfElement(elPressureDeriv,connect);
	

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
	  // std::cerr<<"numEl "<<numEl<<"    presserNoDim  "<< pressureNoDim<<std::endl;
	  ptBubble_[numEl]->SetP(pressureNoDim);
	  //		  std::cerr<<"hier "<<std::endl;
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
	
	  //	if(numEl==90){
	  //(*cla)<< actTime_ << "    " << numEl << "   " <<  bubbleValues_[0]  << "    " <<  bubbleValues_[1]  << "   "
	  //<<  radiusWork_[numEl]   << "    " << velocityWork_[numEl];
	  //}

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
                
          //      if(numEl==90){
          //(*cla)<< actTime_ << "    " << numEl << "   " <<  bubbleValues_[0]  << "    " <<  bubbleValues_[1]  << "   "
          //<<  radiusWork_[numEl]   << "    " << velocityWork_[numEl];
          //}
        
          //set numEl to ODE-Solver
	  //          ptODESolver_->SetNumEl(numEl);
        
	  //          ptODESolver_->Solve(steptime, actTime_, bubbleValues_, *ptBubble_[numEl], hTry_,
	  //                              0, dt);
        
          //set the new values
	  //          radiusWork_[numEl]   = bubbleValues_[0]; 
	  //          velocityWork_[numEl] = bubbleValues_[1]; 


          //      if(numEl==90){
          //        (*cla)<< "       " << numEl << "     " << actTime_ << "    " <<  bubbleValues_[0]  << "    " <<  bubbleValues_[1]  << "   "
          //      <<  radiusWork_[numEl]   << "    " << velocityWork_[numEl] << std::endl;
          //
          //}
	  //Normal case+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          numEl++;
        }  
      }

      //check if we have reached the accuracy!
      performOneMoreStep =    (incrementalErr > incStopCrit_) || (nlRhsNorm > incStopCrit_) ;

    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);



    radius_   = radiusWork_;
    velocity_ = velocityWork_;

  }

  // ======================================================================
  //   Setup structure for bubbly dynamics
  // ======================================================================
  void SolveStepAcousticBubble::SetupBubbleDynamics(){  

    ENTER_FCN( "SolveStepAcousticBubble::SetupBubbleDynamics" );

    ptBubble_.Resize(numPDEElems_);
    bubbleValues_.Resize(2);
    radius_.Resize(numPDEElems_);  
    velocity_.Resize(numPDEElems_);
    radiusWork_.Resize(numPDEElems_);  
    velocityWork_.Resize(numPDEElems_);

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



    for (UInt el=0; el<numPDEElems_; el++) {
      // Choice which bubbledynamical method is used and 
      // creation of pointer to choosen class
      // so far method is chosen in cfs program call with -bubbletyp
      switch(bubbleDynType_){
      case KELLERMIKSIS:
	ptBubble_[el] = new KellerMiksis(initRadius_,density_, sonicVel_, pStatic_, 
					 pVapour_, surfaceTension_, polytrop_, viscosity_);
	//	Info->PrintF( pdename_, "Using Keller-Miksis-Bubble-Model\n");
        break;
     

      case GILMORE:
	//ptBubble_[el] = new Gilmore(initRadius_,density_, sonicVel_, pStatic_, 
	//			pVapour_, surfaceTension_, polytrop_, viscosity_);
	//	Info->PrintF( pdename_, "Using Gilmore-Bubble-Model\n");
	ptBubble_[el] = new Gilmoredimlos(initRadius_,density_, sonicVel_, pStatic_, 
				      pVapour_, surfaceTension_, polytrop_, viscosity_);
	Info->PrintF( pdename_, "Using dimensionless Gilmore-Bubble-Model\n");
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
    // ptODESolver_ = new ODESolver_Rosenbrock;

    //Set bubbledensity
    params->Get( "bubbleNumDensity", bubbleDensity_, "acoustic", "bubbles" );
  }


  void SolveStepAcousticBubble::ComputeBubbleRHS(){  
    ENTER_FCN( "SolveStepAcousticBubble::ComputeBubbleRHS" );

    //rhs-integrator
    Double dummy=1.0;
    BaseForm *rhsForm = new VolumeSrcInt(dummy, isaxi_);        

    UInt numEl = 0;
    for (UInt actDom=0; actDom <  subdoms_.GetSize(); actDom++) {      
      StdVector<Elem*> elemssd;
      ptgrid_->GetVolElems(elemssd, subdoms_[actDom]);
    
    
      for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++)    {              
        BaseFE * ptEl = elemssd[actEl]->ptElem;
        StdVector<UInt> connecth = elemssd[actEl]->connect;

        Matrix<Double> ptCoord;
        GetElemCoords(connecth, ptCoord);
        
        Double beta2 = 1;
                

	bubbleValues_[0] = radiusWork_[numEl];
	bubbleValues_[1] = velocityWork_[numEl];

	//output r,v
	//     if (numEl == 91)
	//	{
	//	  Info->PrintF("","%e  %e  %e\n",actTime_,radius_[numEl],velocity_[numEl]);
	//	  (*data)<< actTime_ << "    " << radius_[numEl] << "    " << velocity_[numEl]<< std::endl;
	//	}	 

	if (actStep_ == 1) 
	  beta2 = 4*PI*bubbleDensity_*6*bubbleValues_[0]
	    *bubbleValues_[1]*bubbleValues_[1]; 
	else {
	  StdVector<Double> dydt(2);


	  // dimensionless**************************************
	  yNoDim_[0] = bubbleValues_[0] / initRadius_;
	  yNoDim_[1] = bubbleValues_[1] / (sqrt( pStatic_/ density_));
	  tNoDim_    = actTime_ / initRadius_ * (sqrt(pStatic_/ density_));
	  ptBubble_[numEl]->CompDeriv(tNoDim_, yNoDim_, dydt);
	  dydt[1] = dydt[1] * pStatic_ / ( density_ * initRadius_ );
	  beta2 = 4*PI*bubbleDensity_*
	    (6*bubbleValues_[0]*bubbleValues_[1]*bubbleValues_[1]
	     + 3*bubbleValues_[0]*bubbleValues_[0]*dydt[1] ); 
	  // dimensionless**************************************

	  //Normal case+++++++++++++++++++++++++++++++++++++++
	  //	  ptBubble_[numEl]->CompDeriv(actTime_, bubbleValues_, dydt);
	  //	  beta2 = 4*PI*bubbleDensity_*
	  //	    (6*bubbleValues_[0]*bubbleValues_[1]*bubbleValues_[1]
	  //	     + 3*bubbleValues_[0]*bubbleValues_[0]*dydt[1] ); 
	  //Normal case+++++++++++++++++++++++++++++++++++++++
	}


	rhsForm->SetFactor(beta2);            
	rhsForm->SetElemPtr(ptEl);
      
	Vector<Double> elemVec;
	rhsForm->CalcElemVector(ptCoord, elemVec);

	// map connect to PDE node numbers
	StdVector<Integer> connect_PDE;
	eqnData_->Node2EQN(connecth, connect_PDE);

	algsys_->SetElementRHS(&elemVec[0], pdeId1_, connect_PDE.GetPointer(), 
			       connect_PDE.GetSize());
	numEl++;
      }
    }

    delete rhsForm;
  

  }


  //  return pointer to vector with first derivative of solution
  const StdVector<Double>& SolveStepAcousticBubble::GetResultData(std::string resultType) {

    ENTER_FCN( "SolveStepAcousticBubble::WriteResults" );

    if ( resultType == "bubbleRadius" ) {
      return radius_;
    }
    else if ( resultType == "bubbleVelocity" ) {
      return velocity_;
    }
    else {
      Error("SolveStepAcousticBubble: Result-Type not possible");
      
      // only to satisfy compiler
      return velocity_;
    }

  }

  

} // end of namespace

