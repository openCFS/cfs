#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "actimeerror.hh"

namespace CoupledField
{

AcousticTimeErrorEstimator::AcousticTimeErrorEstimator(BasePDE * aptPDE)
:TimeErrorEstimator(aptPDE)
{
#ifdef TRACE
  (*trace) << "entering AcousticTimeErrorEstimator::AcousticTimeErrorEstimator" << std::endl;
#endif

  Integer size=ptPDE_->getSize();  
  thirddersol_.Resize(size);
  thirddersol_.Init();

  counter_=0;

  // read tolerance for relative error
  conf->get("tolrelerr",tol_,"Acoustic");
  conf->get("theta", theta_, "Acoustic");

  // read parameters for coarsing strategy
  conf->get("coarsbeta",beta_,"Acoustic");
  conf->get("numrepeat", numrepeat_,"Acoustic");
}

void AcousticTimeErrorEstimator::CalcError(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering AcousticTimeErrorEstimator::CalcError" << std::endl;
#endif

  // calculation of third derivative of solution
  Vector<Double> thirddersolA=(ptPDE_->getS2()-ptPDE_->getS2old())/dt; 
 
  // calculation error of displacement
  Double gamma=ptPDE_->getGamma();
  Double beta=ptPDE_->getBeta();
  
  Vector<Double> errdisplacement;
  Double aux=std::pow(dt,3)/12;
  errdisplacement=(thirddersolA*(1-12*beta)+thirddersol_)*aux; 

  Vector<Double> errvelocity;
  aux=std::pow(dt,3)/6;
  errvelocity=(thirddersolA*(2-6*gamma)+thirddersol_)*aux;

  Double help1=ptPDE_->getAlgSys()->CalcEnergyNorm(0,0,2,errdisplacement.get());
  Double help2=ptPDE_->getAlgSys()->CalcEnergyNorm(0,0,5,errvelocity.get());
  Double error=sqrt(help1+help2);

   help1=ptPDE_->getAlgSys()->CalcEnergyNorm(0,0,2,(ptPDE_->getS()).get());
   help2=ptPDE_->getAlgSys()->CalcEnergyNorm(0,0,5,(ptPDE_->getS1()).get());
   Double normsol=sqrt(help1+help2);

  relativeerror_=error/normsol;
 
  thirddersol_=thirddersolA*2 - thirddersol_;
}

void AcousticTimeErrorEstimator::ChangeStep(Double & dt)
{
#ifdef TRACE
  (*trace) << "entering AcousticTimeErrorEstimator::ChangeStep" << std::endl;
#endif

  Double help1=theta_*tol_/relativeerror_;
  Double help=std::exp(0.333333333*std::log(help1));
  dt*=help;
}

Boolean AcousticTimeErrorEstimator::TestError(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering AcousticTimeErrorEstimator::TestError" << std::endl;
#endif

 CalcError(dt);

 std::cout << "relativeerror" << relativeerror_ << " tolerance " << tol_ << std::endl;

 if (relativeerror_<=beta_*tol_)
   {  counter_++;
      if ( counter_== numrepeat_ )
       { counter_=0; return TRUE;}
   }   
 else if (relativeerror_<= tol_)
       { counter_=0; return FALSE; }
      else { counter_=0; return TRUE;}
}

} // end of namespace
