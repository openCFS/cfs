#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "Forms/forms_header.hh"
#include "newmarkFracDamp.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "basePDE.hh"
#include "DataInOut/MaterialData.hh"
#include "Utils/mathfunctions.hh"

namespace CoupledField {

NewmarkFracDamp::NewmarkFracDamp(std::string apdename, BaseSystem * algebraicsystem,
								 NodeEQN * ptEQN, Grid * aptgrid,
								 BasePDE * aptBasePDE,
								 StdVector<std::string> subdomainList,
								 StdVector<std::string> adampingList, 
								 Integer afracMemory, 
								 InterpolType ainType, Boolean isaxi)
  :TimeStepping(apdename, algebraicsystem, ptEQN)
{
  ENTER_FCN( "NewmarkFracDamp::NewmarkFracDamp" );

  pdename_     = apdename;
  ptgrid_      = aptgrid;
  ptBasePDE_   = aptBasePDE;

  dampingList_ = adampingList;
  subdoms_     = subdomainList;
  fracMemory_  = afracMemory;
  inType_      = ainType;
  isaxi_       = isaxi;
  
  alpha_ = 0.0;
  beta_  = 0.25;
  gamma_ = 0.5;

  //check if integration parameters are defined in conf-file
  std::string analysis;
  params->Get( "type", analysis, "analysis" );
  if(analysis != "paramIdent")
    Info->PrintF( pdename_,"Newmark: Using defaults for alpha, beta and gamma!" );


  Integer numEQNs = ptEQN_->GetNumEQNs();
  Integer dofs = ptEQN_->GetNumDofsPerEQN();

  // get the memory
  solderiv1_.Resize(numEQNs * dofs);  
  solderiv1_.Init();
  solderiv2_.Resize(numEQNs * dofs);  
  solderiv2_.Init();
  

  solpred_.Resize(numEQNs * dofs); 
  solpred_.Init();
  solderiv1pred_.Resize(numEQNs * dofs); 
  solderiv1pred_.Init();
  
  if ( fracMemory_ < 2 )
    Error("Attenuation model needs frac_memory larger than 1",__FILE__,__LINE__);
  else {
	// create vector decribing storing of values when using interpolation
	solMemoryVal_.resize(fracMemory_);

	// get the memory
	solMemory_ = new Vector<Double>[fracMemory_];
	for (Integer i=0; i<fracMemory_; i++) {
	  solMemory_[i].Resize(numEQNs * dofs);  
	  solMemory_[i].Init();
		
	  solMemoryVal_[i] = NOTUSED;
	}
  }

}

NewmarkFracDamp::~NewmarkFracDamp()
{
  ENTER_FCN( "NewmarkFracDamp::~NewmarkFracDamp" );

}

void NewmarkFracDamp::Init(Double * matrix_factors, Double dt)
{
  ENTER_FCN( "NewmarkFracDamp::Init" );

  dt_ = dt;
  CalcParameters(dt_);

  matrix_factors[0] = 1.0;       // factor for stiffness matrix
  matrix_factors[1] = 0.0;       // factor for damping matrix
  matrix_factors[2] = 0.0;       // factor for convection matrix
  matrix_factors[3] = 1.0*a2_;   // factor for mass matrix

}


void NewmarkFracDamp::Predictor(Vector<Double>& solold)
{
  ENTER_FCN( "NewmarkFracDamp::Predictor" );

  laststepcalc_ = ptBasePDE_->GetTimeStepCounter();

  // determine number of terms over which BDF is calculated
  if (laststepcalc_ == 1)  // assumes first nstep = 1 (see transientdriver.cc)!
	calclimit_ = 0;        // no stored values, no BDF has to be computed!

  else if ( laststepcalc_>1 && laststepcalc_<=solMemoryVal_.size() ) {
	calclimit_ = laststepcalc_;
  }
  else {
	calclimit_ = solMemoryVal_.size(); 
	Info->PrintF(pdename_, "solMemoryVal_ size is: %d",  calclimit_);
  }

  solpred_ = solold + solderiv1_*dt_ + solderiv2_*a0_;
  solderiv1pred_ = solderiv1_ + solderiv2_*a1_;
}


void NewmarkFracDamp::UpdateRHS()
{
  ENTER_FCN( "NewmarkFracDamp::UpdateRHS" );

  Vector<Double> coeffMass, coeffDamp;

  // mass part
  coeffMass = solpred_*a2_;
  algsys_->UpdateRHS(MASS,coeffMass.GetPointer());

  // additional term for fractional damping part
  //  if ( dampType_ == FRACTIONAL ) {

  Matrix<Double> elemmat;  
  Matrix<Double> ptCoord;
  BaseFE         * ptElem;
  StdVector<Integer> connecth;
  Vector<Double> rhsvec;
  Vector<Double> elemsol;
  
  Integer level = 0, no_int;
  
  Double density, compressibility, c0, factor;

  MaterialData *mymaterialData;
  mymaterialData = ptBasePDE_->getPDEMaterialData();
  
  for ( Integer actSD=0; actSD < subdoms_.GetSize(); actSD++ ) {

	density         = mymaterialData[actSD].GetDensity();
	compressibility = mymaterialData[actSD].GetCompressibility();
	c0 = sqrt(compressibility/density);

	alpha0_ = mymaterialData[actSD].GetDampingAlfa();
	y_      = mymaterialData[actSD].GetDampingBeta();

	// factor: pre factor for in Newmark timestepping scheme
	// coeff_: weight factors in BDF
	if ( dampingList_[actSD] == "fractional_gl" ) {
	  factor = a2_ * density * 2.0 * alpha0_ / c0 / sin((y_-1.0)*PI/2.0)
		       * exp(-(y_-1.0)*log(dt_));
	  GLWeights(calclimit_, y_);
	}
	else if ( dampingList_[actSD] == "fractional_blank" ) {
	  factor = a2_ * density * 2.0 * alpha0_ / c0 / sin((y_-1.0)*PI/2.0)
		       * exp(-(y_-1.0)*log(dt_)) * exp(-gammaln(1.0- (y_- 1.0)) );
	  BlankWeights(calclimit_, y_, TRUE);
	}

	// get elements belonging to actual subdomain	
	StdVector<Elem*> elemssd;
	ptgrid_->GetElemSD(elemssd,subdoms_[actSD],level);
	
	for (Integer el=0; el < elemssd.GetSize(); el++) {
	  // pointer to element
	  ptElem=elemssd[el]->ptElem;
	  BaseForm * bilinear_mass = new MassInt(ptElem, factor, isaxi_);
	  
	  connecth=elemssd[el]->connect;
	  ptgrid_->GetCoordNodesElemMat(connecth, ptCoord, level);
	  
	  bilinear_mass->CalcElementMatrix(ptCoord, elemmat);
	  
	  // map connect to PDE node numbers
	  StdVector<Integer> connect_PDE;
	  ptEQN_->Node2EQN(connecth, connect_PDE);
	 
	  // compute BDF
	  GetElemSolution(solpred_, elemsol, connect_PDE);
	  rhsvec = -elemsol * coeff_[0];

	  no_int = 0; 
	  for (Integer i=1; i<=calclimit_; i++) {

		switch (solMemoryVal_[i-1]) {
		case TRUEVAL:
		  GetElemSolution(solMemory_[i-1], elemsol, connect_PDE);
		  rhsvec += elemsol * coeff_[i];
		  break;
		case LIN1PT:
		  no_int++;
		  GetElemSolution(solMemory_[i-1-no_int], elemsol, connect_PDE);
		  rhsvec +=  solMemory_[i-1-no_int] * 0.5 * coeff_[i];
		  GetElemSolution(solMemory_[i-1-no_int+1], elemsol, connect_PDE);
		  rhsvec += solMemory_[i-1-no_int+1] * 0.5 * coeff_[i];
		  break;
		default:
		  break;
		}
	  }

	  rhsvec = elemmat * elemsol;
	  
	  //assemble to RHS
	  algsys_->SetElementRHS(&rhsvec[0], connect_PDE.GetPointer(),
							 connect_PDE.GetSize());
	  
	  delete bilinear_mass;	  
	}

	PrintSolMemoryVal(actSD);
  }
  //  }

}


void NewmarkFracDamp::Corrector(Vector<Double>& solnew)
{
  ENTER_FCN( "NewmarkFracDamp::Corrector" );

  solderiv2_ = (solnew - solpred_) * a2_;
  solderiv1_ = solderiv1pred_ + solderiv2_*a3_;

  // store function values, reverse storing!!!!
  // solMemory_[0]=p_n, solMemory_[1]=p_(n-1)...solMemory_[fracMemory_ - 1]=p_(n-fracMemory_+1)
  //  if ( dampType_ == FRACTIONAL ) {

  solMemory_[0] = solnew - solpred_;

  if (inType_ == NOTUSED) {
	for ( Integer i=fracMemory_-1; i>=1; i-- ) {
	  solMemory_[i] = solMemory_[i-1];
	  solMemoryVal_[i] = solMemoryVal_[i-1];
	}
	solMemoryVal_[0] = TRUEVAL;
  }

  // Linear interpolation of one solution value	
  else if (inType_ == LIN1PT) {

	if (laststepcalc_ <= fracMemory_) {
	  for (Integer i=fracMemory_-1; i>=1; i--) {
		solMemory_[i] = solMemory_[i-1];
		solMemoryVal_[i] = solMemoryVal_[i-1];
	  }
	  solMemoryVal_[0] = TRUEVAL;
	}
	else if ( (laststepcalc_ > fracMemory_) && (laststepcalc_ < 2*fracMemory_) ) {
	  for (Integer i=(2*fracMemory_-laststepcalc_-2); i>=1; i--)
		solMemory_[i] = solMemory_[i-1];
		
	  solMemoryVal_[0] = TRUEVAL;
	  solMemoryVal_.insert( (solMemoryVal_.begin()+2*fracMemory_-laststepcalc_-1),LIN1PT);
	  Info->PrintF(pdename_, "Middle section, timestep: ",laststepcalc_ );
	}
	else if (laststepcalc_ >= 2*fracMemory_ ) {
	  if ( (laststepcalc_%2)==1 ) {
		for (Integer i=fracMemory_-1; i>=1; i--)
		  solMemory_[i] = solMemory_[i-1];

		solMemoryVal_.pop_back();
		solMemoryVal_.pop_back();
		solMemoryVal_.insert(solMemoryVal_.begin(), TRUEVAL);
	  }
	  else if ( (laststepcalc_%2)==0 )
		solMemoryVal_.insert( (solMemoryVal_.begin()+1), LIN1PT);
	}
  }
  //  }

}

void NewmarkFracDamp::CalcParameters(Double dt)
{
  ENTER_FCN( "NewmarkFracDamp::CalcParameters" );

  //for predictors
  a0_ = 0.5*(1-2.0*beta_)*dt_*dt_;
  a1_ = (1-gamma_)*dt_;

  //for correctors, matrices
  a2_ = 1.0/(beta_*dt_*dt_);
  a3_ = gamma_*dt_;

  //for RHS, Matrices
  a4_ = gamma_ / (beta_*dt_);
}


void NewmarkFracDamp::GetElemSolution ( const Vector<Double>& sol, 
									Vector<Double>& elemsol,
									const StdVector<Integer> & connectPDE )
{
  ENTER_FCN( "NewmarkFracDamp::GetElemSolution" );
   
  elemsol.Resize(connectPDE.GetSize());
  for (Integer eqn=0; eqn<connectPDE.GetSize(); eqn++) 
    elemsol[eqn] = sol[connectPDE[eqn]-1];
      
}

void NewmarkFracDamp::GLWeights(Integer memory, Double y )
{
  ENTER_FCN( "NewmarkFracDamp::GLWeights" );
   
  // reserve memory for weights of BDF, order of derivative is y-1
  coeff_.resize(memory+1);
   
  coeff_[0] = 1.0; // Index 0

  for (Integer i=1; i <= memory; i++)
    coeff_[i] = coeff_[i-1] * ((i-1) - (y_-1))/i;

}


void NewmarkFracDamp::BlankWeights(Integer memory, Double y, Boolean full)
{
  ENTER_FCN( "NewmarkFracDamp::BlankWeights" );
 
  Double pot;
  pot = 1.0- (y - 1.0);

  // reserve memory for weights of BDF, order of derivative is y-1
  coeff_.resize(memory+1);
      
  coeff_[0] = 1.0/ pot; // Index 0
 
  if (full) {
	for (Integer i=1; i<memory; i++) // Index 1 .. (memory-1)
	  coeff_[i]= ( exp(pot* log(i-1))- 2* exp(pot* log(i))+ exp(pot* log(i+1)) )
		/ pot;
	// Index memory
	coeff_[memory]= exp(-(y_-1.0)*log(fracMemory_))+ 
	  ( exp(pot* log(fracMemory_ -1))- exp(pot* log(fracMemory_)) ) / pot;
  }
  else {
	for (Integer i=1; i<=memory; i++)    // Index 1 .. memory
	  coeff_[i]= ( exp(pot* log(i-1))- 2* exp(pot* log(i))+ exp(pot* log(i+1)) )
		/ pot;
  }
}

void NewmarkFracDamp::PrintSolMemoryVal(Integer actSD)
{
  ENTER_FCN( "NewmarkFracDamp::PrintSolMemoryVal" );

  Info->PrintF(pdename_, "Region %d, solMemory_ of timestep %d is:", actSD, laststepcalc_ );
  std::string msg;
  for (Integer i=0; i<solMemoryVal_.size(); i++) {
	switch (solMemoryVal_[i]) {
	case NOTUSED:
	  msg += " N";
	  break;
	case TRUEVAL:
	  msg += " T";
	  break;
	case LIN1PT:
	  msg += " L";
	  break;
	}
  }
  Info->PrintF(pdename_, msg.c_str() );

}

} // end of namespace
