#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "Forms/forms_header.hh"
#include "newmarkdamp.hh"

namespace CoupledField
{

NewmarkDamp :: NewmarkDamp(std::string apdename, BaseSystem * algebraicsystem, NodeEQN * ptEQN, 
			   Grid * aptgrid, StdVector<std::string> adampingList, 
			   StdVector<std::string> subdomainList, Integer adamping, 
			   Integer afrac_memory, Double damp, Boolean isaxi)
:TimeStepping(apdename, algebraicsystem, ptEQN)
{
  ENTER_FCN( "NewmarkDamp::NewmarkDamp" );

  ptgrid_      = aptgrid;
  dampingList_ = adampingList;
  damping_     = adamping;
  frac_memory_ = afrac_memory;
  y_           = damp;
  subdoms_     = subdomainList;
  isaxi_       = isaxi;
  
  alpha_ = 0.0;
  beta_  = 0.25;
  gamma_ = 0.5;

  //check if integration parameters are defined in conf-file
  conf->ifget("alpha_NM",alpha_,pdename_); 
  conf->ifget("beta_NM",beta_,pdename_); 
  conf->ifget("gamma_NM",gamma_,pdename_);

  Integer numEQNs = ptEQN_->GetNumEQNs();
  Integer dofs = ptEQN_->GetNumDofsPerEQN();

  //get the memory
  solderiv1_.Resize(numEQNs * dofs);  
  solderiv1_.Init();
  solderiv2_.Resize(numEQNs * dofs);  
  solderiv2_.Init();
  

  solpred_.Resize(numEQNs * dofs); 
  solpred_.Init();
  solderiv1pred_.Resize(numEQNs * dofs); 
  solderiv1pred_.Init();
  
  if (frac_memory_==0)
    Error("Attenuation model needs frac_memory larger than 0",__FILE__,__LINE__);

  //get the memory
  solfrac_ = new Vector<Double>[frac_memory_];
  for (Integer i=0; i<frac_memory_; i++)
    {
      solfrac_[i].Resize(numEQNs * dofs);  
      solfrac_[i].Init();
    }
}

NewmarkDamp :: ~NewmarkDamp()
{
  ENTER_FCN( "NewmarkDamp::~NewmarkDamp" );

}

void NewmarkDamp::Init(Double * matrix_factors, Double dt)
{
  ENTER_FCN( "NewmarkDamp::Init" );

  dt_ = dt;
  CalcParameters(dt_);

  matrix_factors[0] = 1.0;       // factor for stiffness matrix

  matrix_factors[1] = 0.0; 

  //  if (damping_)
  //    matrix_factors[1] = 1/pow(dt_,(y_+1.0)); // factor for damping matrix

  matrix_factors[2] = 0.0;       // factor for convection matrix
  matrix_factors[3] = 1.0*a2_;   // factor for mass matrix

}


void NewmarkDamp::Predictor(Vector<Double>& solold)
{
  ENTER_FCN( "NewmarkDamp::Predictor" );

  solpred_ = solold + solderiv1_*dt_ + solderiv2_*a0_;
  solderiv1pred_ = solderiv1_ + solderiv2_*a1_;
}


void NewmarkDamp::UpdateRHS()
{
  ENTER_FCN( "NewmarkDamp::UpdateRHS" );

  Vector<Double> coeffMass, coeffDamp;

  // mass part
  coeffMass = solpred_*a2_;
  algsys_->UpdateRHS(MASS,coeffMass.GetPointer());


  // damping part
  //  Double coeff;
  //  coeff = 1;


//   for (Integer i=0; i<frac_memory_; i++)
//     {
//       coeff *= (i - y_-1)/(i+1);
//       coeffDamp = solfrac_[i] * coeff/pow(dt_,(y_+1));
//       algsys_->UpdateRHS(DAMPING,coeffDamp.get());
//     }


  Matrix<Double> elemmat;  
  Matrix<Double> ptCoord;
  BaseFE         * ptElem;
  
  StdVector<Integer> connecth;
  Vector<double> rhsvec;

  Integer level = 0;
  for (Integer sd=0; sd<subdoms_.GetSize(); sd++)
    {
      if ( dampingList_[sd] == "fractional") {
	// compute factor for subdomain
	Double factor = 1.0;
	
	StdVector<Elem*> elemssd;
	ptgrid_->GetElemSD(elemssd,subdoms_[sd],level);
	
	for (Integer el=0; el < elemssd.GetSize(); el++) {  
	  ptElem=elemssd[el]->ptElem;
	  BaseForm * bilinear_mass = new MassInt(ptElem, factor, isaxi_);
	  
	  connecth=elemssd[el]->connect;
	  ptgrid_->GetCoordNodesElemMat(connecth, ptCoord, level);
	  
	  bilinear_mass->CalcElementMatrix(ptCoord, elemmat);
	  
	  // map connect to PDE node numbers
	  StdVector<Integer> connect_PDE;
	  ptEQN_->Node2EQN(connecth, connect_PDE);

	  Vector<Double> elemsol;
	  Integer idx = 0;
	  GetElemSolution(solfrac_[idx], elemsol, connect_PDE);
	  
	  rhsvec = elemmat * elemsol;
	  
	  //assemble to RHS
	  algsys_->SetElementRHS(&rhsvec[0], connect_PDE.GetPointer(), connect_PDE.GetSize());
	  
	  delete bilinear_mass;	  
	}  
      }
    }

   
}


void NewmarkDamp::Corrector(Vector<Double>& solnew)
{
  ENTER_FCN( "NewmarkDamp::Corrector" );

  solderiv2_ = (solnew - solpred_) * a2_;
  solderiv1_ = solderiv1pred_ + solderiv2_*a3_;

  //solfrac[0]: p_n; solfrac[1]: p_(n-1); ..
  for (Integer i=frac_memory_-1; i>0; i--)
    solfrac_[i] = solfrac_[i-1];
      
  solfrac_[0] = solnew;
}

void NewmarkDamp :: CalcParameters(Double dt)
{
  ENTER_FCN( "NewmarkDamp::CalcParameters" );

  //for predictors
  a0_ = 0.5*(1-2.0*beta_)*dt_*dt_;
  a1_ = (1-gamma_)*dt_;

  //for correctors, matrices
  a2_ = 1.0/(beta_*dt_*dt_);
  a3_ = gamma_*dt_;

  //for RHS, Matrices
  a4_ = gamma_ / (beta_*dt_);
}


void NewmarkDamp :: GetElemSolution (const Vector<Double>& sol, Vector<Double>& elemsol, const StdVector<Integer> & connectPDE)
{
  ENTER_FCN( "NewmarkDamp::GetElemSol" );

  for (Integer eqn=0; eqn<connectPDE.GetSize(); eqn++) 
    elemsol[eqn] = sol[connectPDE[eqn]-1];
      
}

} // end of namespace
