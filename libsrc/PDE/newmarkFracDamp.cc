#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

#include "Forms/forms_header.hh"
#include "newmarkFracDamp.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "basePDE.hh"
#include "DataInOut/MaterialData.hh"
#include "Utils/mathfunctions.hh"

namespace CoupledField {

NewmarkFracDamp::NewmarkFracDamp(std::string apdename,
								 BaseSystem * algebraicsystem,
								 NodeEQN * ptEQN, Grid * aptgrid,
								 BasePDE * aptBasePDE,
								 StdVector<std::string> asubdomainList,
								 StdVector<std::string> adampingList, 
								 Integer afracMemory, 
								 InterpolType ainType, Boolean isaxi)
  :TimeStepping(apdename, algebraicsystem, ptEQN)
{
  ENTER_FCN( "NewmarkFracDamp::NewmarkFracDamp" );

  pdename_     = apdename;
  ptgrid_      = aptgrid;
  ptBasePDE_   = aptBasePDE;

  subdoms_     = asubdomainList;
  dampingList_ = adampingList;
  fracMemory_  = afracMemory;
  inType_      = ainType;
  isaxi_       = isaxi;
  
  alpha_ = 0.0;
  beta_  = 0.25;
  gamma_ = 0.5;

  //check if integration parameters are defined
  std::string analysis;
  params->Get( "type", analysis, "analysis" );
  if(analysis != "paramIdent")
    Info->PrintF( pdename_,"Newmark: Using defaults for alpha, beta and gamma!\n" );


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
  
  if ( fracMemory_ <= 1 )
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
  //  matrix_factors[1] = 1.0*a2_;
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
  else if ( laststepcalc_>1 && laststepcalc_<=solMemoryVal_.size() )
	calclimit_ = laststepcalc_ - 1;
  else
	calclimit_ = solMemoryVal_.size(); 

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
  Matrix<Double> elemmat;         // element matrix
  Matrix<Double> ptCoord;
  BaseFE         * ptElem;
  StdVector<Integer> connecth;
  Vector<Double> rhsAssemble, rhsvec, elemsol; 
  
  Integer level = 0, noInt;
  
  Double density, compressibility, c0, alpha0, y, factor;

  MaterialData *mymaterialData;
  mymaterialData = ptBasePDE_->getPDEMaterialData();

  if ( dampingList_.GetSize() != subdoms_.GetSize() )
  	Error("Mismatch between dampingList_ and subdoms_!", __FILE__, __LINE__);  
  for ( Integer actSD=0; actSD < subdoms_.GetSize(); actSD++ ) {

	density         = mymaterialData[actSD].GetDensity();
	compressibility = mymaterialData[actSD].GetCompressibility();
	c0 = sqrt(compressibility/density);
	alpha0 = mymaterialData[actSD].GetDampingAlfa();
	y      = mymaterialData[actSD].GetDampingBeta();

	// factor: pre factor in Newmark timestepping scheme
	// coeff_: weight factors in BDF
	if ( dampingList_[actSD] == "fractional_gl" ) {
	  factor = density * 2.0 * alpha0 / c0 / sin((y-1.0)*PI/2.0);
	  factor *= exp(-(y-1.0)*log(dt_));
	  factor *= a2_;
	  GLWeights(calclimit_, y);
	}
	else if ( dampingList_[actSD] == "fractional_blank" ) {
	  factor =  density * 2.0 * alpha0 / c0 / sin((y-1.0)*PI/2.0);
	  factor *= exp(-(y-1.0)*log(dt_)) * exp(-gammaln(1.0- (y- 1.0)) );
	  factor *= a2_;
	  BlankWeights(calclimit_, y, TRUE);
	}

// 	// output coeff_
// 	std::cerr << "Weight factors of BDF are:" << std::endl;
// 	for (Integer k=0; k<coeff_.size(); k++)
// 	  std::cerr << "          w(" << k << ") = " << coeff_[k] << std::endl;

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

#ifdef DEBUG
	  // output matrix with which BDF is computed
	  (*debug) << "Mass Matrix, damping part RHS, of Element " << el << std::endl;
	  (*debug) << elemmat << std::endl;
#endif
	  // compute BDF
	  GetElemSolution(solpred_, elemsol, connect_PDE);
	  rhsvec = -elemsol * coeff_[0];
	  noInt = 0; // counts how many values have allready been interpolated
	  for (Integer i=1; i<=calclimit_; i++) {

		if ( solMemoryVal_[i-1] == TRUEVAL ) {
		  GetElemSolution(solMemory_[i-1], elemsol, connect_PDE);
		  rhsvec += elemsol * coeff_[i];
		}
		else if ( solMemoryVal_[i-1] == LIN1PT ) {
		  noInt++;
		  GetElemSolution(solMemory_[i-1-noInt], elemsol, connect_PDE);
		  rhsvec +=  solMemory_[i-1-noInt] * 0.5 * coeff_[i];
		  GetElemSolution(solMemory_[i-1-noInt+1], elemsol, connect_PDE);
		  rhsvec += solMemory_[i-1-noInt+1] * 0.5 * coeff_[i];
		}
		else
		  Error("Something went wrong in calculation of BDF!",__FILE__,__LINE__);
	  }

	  rhsAssemble = elemmat * rhsvec;

#ifdef DEBUG
	  (*debug) <<  "BDF vector of timestep " << laststepcalc_ << std::endl;
	  (*debug) << rhsvec << std::endl;
#endif
	  
	  //assemble to RHS
	  algsys_->SetElementRHS(&rhsAssemble[0], connect_PDE.GetPointer(),
							 connect_PDE.GetSize());
	  
	  delete bilinear_mass;	  
	}
	Info->PrintF(pdename_, "timestep: %d, BDF over %d calculated\n"
				 , laststepcalc_, calclimit_ );
	PrintSolMemoryVal();
  }
}


void NewmarkFracDamp::Corrector(Vector<Double>& solnew)
{
  ENTER_FCN( "NewmarkFracDamp::Corrector" );

  solderiv2_ = (solnew - solpred_) * a2_;
  solderiv1_ = solderiv1pred_ + solderiv2_*a3_;

  // store function values, reverse storing!!!!
  // solMemory_[0]=p_n, solMemory_[1]=p_(n-1)...solMemory_[fracMemory_ - 1]=p_(n-fracMemory_+1)
  // no interpolation
  if (inType_ == NOTUSED) {
	for ( Integer i=fracMemory_-1; i>=1; i-- ) {
	  solMemory_[i] = solMemory_[i-1];
	  solMemoryVal_[i] = solMemoryVal_[i-1];
	}
	solMemory_[0] = solnew - solpred_;
	solMemoryVal_[0] = TRUEVAL;
  }
  // linear interpolation of one solution value	
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
	  Info->PrintF(pdename_, "Middle section, timestep: %d\n",laststepcalc_ );
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
	solMemory_[0] = solnew - solpred_;
  }
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

  for (Integer i=1; i <= memory; i++) // Index 1 .. memory
    coeff_[i] = coeff_[i-1] * ((i-1) - (y-1))/i;
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
	coeff_[memory]= exp(-(y-1.0)*log(fracMemory_))+ 
	  ( exp(pot* log(fracMemory_ -1))- exp(pot* log(fracMemory_)) ) / pot;
  }
  else {
	for (Integer i=1; i<=memory; i++) // Index 1 .. memory
	  coeff_[i]= ( exp(pot* log(i-1))- 2* exp(pot* log(i))+ exp(pot* log(i+1)) )
		/ pot;
  }
}

void NewmarkFracDamp::PrintSolMemoryVal()
{
  ENTER_FCN( "NewmarkFracDamp::PrintSolMemoryVal" );

  std::string msg="solMemory_ is:";
  for (Integer i=0; i<solMemoryVal_.size(); i++) {
	if ( solMemoryVal_[i] == NOTUSED )
	  msg += " N";
	else if ( solMemoryVal_[i] == TRUEVAL )
	  msg += " T";
	else if ( solMemoryVal_[i] ==  LIN1PT )
	  msg += " L";
  }
  msg += "\n";
  Info->PrintF(pdename_, msg.c_str() );
}

} // end of namespace
