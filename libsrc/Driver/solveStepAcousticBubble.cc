#include <fstream>
#include <iostream>
#include <string>

#include "solveStepAcousticBubble.hh"
#include "assemble.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Forms/forms_header.hh"
#include "General/environment.hh"
#include "PDE/timestepping.hh"

namespace CoupledField {

SolveStepAcousticBubble::SolveStepAcousticBubble(StdPDE& apde,
						 BubbleDynType bubbleDynType) 
  : StdSolveStep(apde) 
{
  ENTER_FCN( "SolveStepAcousticBubble::SolveStepAcousticBubble" );
  std::cerr << ":SolveStepAcousticBubble" << std::endl;
  
  // setup the structure for bubble dynamics
  // query if cavitation model should be used will be added after setup of new
  // xml file 
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


void SolveStepAcousticBubble::SolveStepTrans( const Integer kstep, const Double asteptime, 
					      const Integer level, const Boolean reset ) {

  ENTER_FCN( "SolveStepAcousticBubble::SolveStepTrans" );

  lasttimecalc_ = asteptime;
  laststepcalc_ = kstep;

  StepTransBubble(kstep, asteptime, level, reset);

}


void SolveStepAcousticBubble::StepTransBubble(const Integer kstep, const Double asteptime,
					      const Integer level, const Boolean reset) {

  ENTER_FCN( "SolveStepAcousticBubble::StepTransBubble" );

  Double * ptsol;
  Integer job;


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

  if (laststepcalc_ == 1) {
    job = 1;
    assemble_->AssembleMatrices(level);
    algsys_->ConstructEffectiveMatrix(matrix_factor_);
  }  

  //perform predictor and update RHS according to time stepping algorithm
  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
  Vector<Double> & sol4pred= solhelp->GetAlgSysVector();
  TS_alg_->Predictor(sol4pred);

  //set old solution  
  oldSol = sol4pred;


  //do nonlinear iteration
  Integer iterationCounter=0;
  Boolean performOneMoreStep;
  do {
    iterationCounter++;

    // for every time step write out number of iteration loops to standard out
    //    if (iterationCounter == 1)
    // std::cout << std::endl << "Time step:   "  << kstep 
    //	<< "  ,Iterations: " << iterationCounter;
    //else 
    // std::cout	<< "  " << iterationCounter;
    
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
    for (Integer i=0; i<nlRhsNew.GetSize(); i++) {
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
    assemble_->AssembleSrcRHS(level,lasttimecalc_);
    
    
    //acoust for time stepping
    TS_alg_->UpdateRHS();
    
    
    SetBCs(level,lasttimecalc_);

    if ( iterationCounter>1 ) 
      job = 3;

    algsys_->BuildInDirichlet();
    if ( job == 1 ) {
      algsys_->SetupSolver(job);
      algsys_->SetupPrecond(job);
    }
	
    //solve the PDE system
    algsys_->Solve();
    
    // store new solution in newSol
    algsys_->GetSolutionVal( ptsol );
    StoreAlgsysToVec(newSol, ptsol );
    
    // perform corrector step, if effective mass formulation is used,
    //   we need the Corrector step before we store newsol to sol_,
    //   because newsol is second time derivative at first!
    TS_alg_->Corrector(newSol);
    
    //put new solution to sol_
    sol_->SetAlgSysVector(newSol);  


    
    // compute L2-Norm of error between last incremental solution and
    //   actual incremental solution
    Double solIncrL2Norm=0;
    for (Integer i=0; i<newSol.GetSize(); i++)
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
    StdVector<Integer> connect, Eqns;  
    Integer numEl=0;
    
    for (Integer actSD=0; actSD<subdoms_.GetSize(); actSD++) {
      StdVector<Elem*> elemssd;
      ptgrid_->GetElemSD(elemssd,subdoms_[actSD],actlevel_);
      
      for (Integer j=0; j < elemssd.GetSize(); j++) {
	
	ptElem  = elemssd[j]->ptElem;
	connect = elemssd[j]->connect;



        // Output of nlRHSNew for all nodes of desired element
	//if (numEl == 90){
	// Integer eqnNr = 0;
	// Integer eqnDof = 0;
	// Integer dof = 1;
	// std::cerr << lasttimecalc_ << "   "  << numEl; 
	// for ( Integer i=0; i < connect.GetSize(); i++){
	//   eqnData_ -> Node2EQN(connect[i],dof, eqnNr, eqnDof);
	//   std::cerr << "             " << connect[i];
	//   if (eqnNr != 0){
	//     std::cerr << "    " << nlRhsNew[ (eqnNr - 1)];
	//   }
	// }
	// std::cerr << "   " << std::endl;
	//}
		


	Vector<Double> elPressure;
	Vector<Double> elPressureDeriv;
	
	//========= we now have also in coupled computation a pressure formulation =============
	
	//get the pressure
	sol_->GetElemSolution(elPressure,connect);

	
	//get 1st derivative of pressure
	GetDerivSolVecOfElement(elPressureDeriv,connect);
	

	//compute average values
	Double pressure=0;
	Double pressureDeriv=0;
	for (Integer elnode=0; elnode<elPressure.GetSize(); elnode++) {
	  pressure      += elPressure[elnode];
	  pressureDeriv += elPressureDeriv[elnode];
	}
	
	pressure      /= elPressure.GetSize();
	pressureDeriv /= elPressureDeriv.GetSize();

	if(numEl==90){
	  (*cla)<< lasttimecalc_ << "    " << "numEl " <<  "   "
		<< pressure << "    " << pressureDeriv  << "   "
		<< pressure-oldpressure << "    " 
		<< pressureDeriv-oldpressureDeriv << std::endl;
	  oldpressure = pressure;
	  oldpressureDeriv = pressureDeriv;
	}
        
	
	//set to ODE

	//Dimensionless pressure to bubble
	//	pressure = pressure / 1e5;
	//pressureDeriv = pressureDeriv * 10e-6 /1e5 /(sqrt(1e5/998));

	ptBubble_[numEl]->SetP(pressure);
	ptBubble_[numEl]->SetDpdt(pressureDeriv);
	
	// In case of explicit Euler watch out suggested stepsize
	Double dt = TS_alg_->GetTimeStep();
	Double steptime = lasttimecalc_ - dt;
	
	//get the current values


	bubbleValues_[0] = radius_[numEl];
	bubbleValues_[1] = velocity_[numEl];
	
	//Dimensionless radius etc. to bubble
	//bubbleValues_[0] = bubbleValues_[0] / 10e-6;
	//bubbleValues_[1] = bubbleValues_[1] / (sqrt(1e5/998));
	
	//	if(numEl==90){
	//(*cla)<< lasttimecalc_ << "    " << numEl << "   " <<  bubbleValues_[0]  << "    " <<  bubbleValues_[1]  << "   "
	//<<  radiusWork_[numEl]   << "    " << velocityWork_[numEl];
	//}


	
	//set numEl to ODE-Solver
	ptODESolver_->SetNumEl(numEl);
	
	ptODESolver_->Solve(steptime, lasttimecalc_, bubbleValues_, *ptBubble_[numEl], dt / 3.0,
			    0, dt);
	
	//set the new values
	radiusWork_[numEl]   = bubbleValues_[0]; //* (sqrt(1e5/998)) ;
	velocityWork_[numEl] = bubbleValues_[1]; //* 1e5 / 10e-6 / 998;

	//	if(numEl==90){
	//  	  (*cla)<< "       " << numEl << "     " << lasttimecalc_ << "    " <<  bubbleValues_[0]  << "    " <<  bubbleValues_[1]  << "   "
	//	<<  radiusWork_[numEl]   << "    " << velocityWork_[numEl] << std::endl;
	//
	//}
	
	numEl++;
      }  
    }

    //check if we have reached the accuracy!
    performOneMoreStep =    (incrementalErr > incStopCrit_) || (nlRhsNorm > incStopCrit_) ;

  } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);


  radius_   = radiusWork_;
  velocity_ = velocityWork_;

  (*data) << lasttimecalc_ << "   " << "90 " <<  "   " << radius_[90] << "    " << velocity_[90] << std::endl;
  //   std::cerr << lasttimecalc_ << "   " << "0 " <<  "   " << radius_[0] << "    " << velocity_[0] << std::endl;
  
  // if ( iterationCounter = nonLinMaxIter_ ){
  //    Error("Nonlinear Iteration not converged", __FILE__, __LINE__ );
  //  }
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

    Double initRadius = auxVec[0];
    Double initVel = 0.0;
    params->Get( "initVel", initVel, "bubbleDynamic" );
    Double density = 0.0;
    params->Get( "density", density, "bubbleDynamic" );
    Double sonicVel = 0.0;
    params->Get( "sonicVel", sonicVel, "bubbleDynamic" );
    Double pStatic = 0.0;
    params->Get( "pStatic", pStatic, "bubbleDynamic" );
    Double pVapour = 0.0;
    params->Get( "pVapour", pVapour, "bubbleDynamic" );
    Double surfaceTension = 0.0;
    params->Get( "surfaceTension", surfaceTension, "bubbleDynamic" );
    Double polytrop = 0.0;
    params->Get( "polytrop", polytrop, "bubbleDynamic" );
    Double viscosity = 0.0;
    params->Get( "viscosity", viscosity, "bubbleDynamic" );


  for (Integer el=0; el<numPDEElems_; el++) {
    // Choice which bubbledynamical method is used and 
    // creation of pointer to choosen class
    // so far method is chosen in cfs program call with -bubbletyp
    switch(bubbleDynType_){
    case KELLERMIKSIS:
	ptBubble_[el] = new KellerMiksis(initRadius,density, sonicVel, pStatic, 
				     pVapour, surfaceTension, polytrop, viscosity);
	//	Info->PrintF( pdename_, "Using Keller-Miksis-Bubble-Model\n");

	break;
     
    case GILMORE:
 	ptBubble_[el] = new Gilmore(initRadius,density, sonicVel, pStatic, 
				pVapour, surfaceTension, polytrop, viscosity);
	//	Info->PrintF( pdename_, "Using Gilmore-Bubble-Model\n");
	
	break;

    default:
      Error("No bubblemethod specified ",__FILE__,__LINE__);
    }
  

    radius_[el]       = initRadius;
    velocity_[el]     = initVel;
    radiusWork_[el]   = initRadius;
    velocityWork_[el] = initVel;
    
  }

    // Generate ODE solver object
    ptODESolver_ = new ODESolver_RKF45;

    //Set bubbledensity
    params->Get( "bubbleNumDensity", bubbleDensity_, "acoustic", "bubbles" );
}


void SolveStepAcousticBubble::ComputeBubbleRHS(){  
  ENTER_FCN( "SolveStepAcousticBubble::ComputeBubbleRHS" );


  //rhs-integrator
  Double dummy=1.0;
  BaseForm *rhsForm = new VolumeSrcInt(dummy, isaxi_);        

  Integer numEl = 0;
  for (Integer actDom=0; actDom <  subdoms_.GetSize(); actDom++) {      
    StdVector<Elem*> elemssd;
    ptgrid_->GetElemSD(elemssd, subdoms_[actDom], actlevel_);
    
    
    for (Integer actEl=0; actEl< elemssd.GetSize(); actEl++)    {              
      BaseFE * ptEl = elemssd[actEl]->ptElem;
      StdVector<Integer> connecth = elemssd[actEl]->connect;

      Matrix<Double> ptCoord;
      GetElemCoords(connecth, ptCoord, actlevel_);
        
      Double beta2 = 1;
                
      bubbleValues_[0] = radiusWork_[numEl];
      bubbleValues_[1] = velocityWork_[numEl];

      //output r,v
      //     if (numEl == 91)
      //	{
	  //	  Info->PrintF("","%e  %e  %e\n",lasttimecalc_,radius_[numEl],velocity_[numEl]);
      //	  (*data)<< lasttimecalc_ << "    " << radius_[numEl] << "    " << velocity_[numEl]<< std::endl;
      //	}	 

      if (laststepcalc_ == 1) 
        beta2 = 4*PI*bubbleDensity_*6*bubbleValues_[0]
          *bubbleValues_[1]*bubbleValues_[1]; 
      else {
        StdVector<Double> dydt(2);
        ptBubble_[numEl]->CompDeriv(lasttimecalc_, bubbleValues_, dydt);

        beta2 = 4*PI*bubbleDensity_*
          (6*bubbleValues_[0]*bubbleValues_[1]*bubbleValues_[1]
           + 3*bubbleValues_[0]*bubbleValues_[0]*dydt[1] ); 
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
  }

}

  

} // end of namespace

