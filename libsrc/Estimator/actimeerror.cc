#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

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

  thirddersoldt_.Resize(size);
  thirddersoldt_.Init();

  counter_=0;
  maxrelativeerror_=0;
  maxnormsol_=0;
  maxnorml2sol_=0;

  Calc3DerFromEquation_=FALSE;

  // read tolerance for relative error
  conf->get("tolrelerr",tol_,"Acoustic");
  conf->get("theta", theta_, "Acoustic");

  // read parameters for coarsing strategy
  conf->get("coarsbeta",beta_,"Acoustic");
  conf->get("numrepeat", numrepeat_,"Acoustic");

  // read bounds for timestep
  conf->get("maxdt",maxdt_,"Acoustic");
  conf->get("mindt",mindt_,"Acoustic"); 
}

void AcousticTimeErrorEstimator::CalcError(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering AcousticTimeErrorEstimator::CalcError" << std::endl;
#endif

  Double gamma=ptPDE_->getGamma();
  Double beta=ptPDE_->getBeta();

  Vector<Double> thirddersolAdt=(ptPDE_->getS2()-ptPDE_->getS2old());
  Integer size=ptPDE_->getSize();
  Double help=sqrt(thirddersolAdt*thirddersolAdt/size)*dt*dt*0.0833;
 
  Vector<Double> solution=ptPDE_->getS();
  Double norm=solution.norm_2();

  if (norm>maxnorml2sol_) maxnorml2sol_=norm;

  Double relativeerror;
  if (maxnorml2sol_ > 1e-10) relativeerror=help/norm; 
  else relativeerror=0.0;   

  maxrelativeerror_=relativeerror;
}

/*
void AcousticTimeErrorEstimator::CalcError(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering AcousticTimeErrorEstimator::CalcError" << std::endl;
#endif

  Double gamma=ptPDE_->getGamma();
  Double beta=ptPDE_->getBeta();

  // calculation of third derivative of solution
  Vector<Double> thirddersolAdt=(ptPDE_->getS2()-ptPDE_->getS2old()); 

  Vector<Double> help=ptPDE_->getS();
  Double norml2sol=help.norm_2();
 
  if (norml2sol>maxnorml2sol_) maxnorml2sol_=norml2sol;
 
  Double N=ptPDE_->getSize();
  Double errorl2=sqrt(thirddersolAdt*thirddersolAdt/N)*dt*dt*0.0833;

  std::cout << " Error " << errorl2 << "  " << norml2sol << std::endl;

  // calculation error of displacement
  Vector<Double> errdisplacement;
  Double aux=pow(dt,2)/12;
  errdisplacement=(thirddersolAdt*(1-12*beta)+thirddersoldt_)*aux; 

  Vector<Double> errvelocity;
  aux=pow(dt,1)/6;
  errvelocity=(thirddersolAdt*(2-6*gamma)+thirddersoldt_)*aux;

  Double help1=ptPDE_->getAlgSys()->CalcEnergyNorm(0,0,2,errdisplacement.get());
  Double help2=ptPDE_->getAlgSys()->CalcEnergyNorm(0,0,5,errvelocity.get());
  Double error=sqrt(help1+help2);

  help1=ptPDE_->getAlgSys()->CalcEnergyNorm(0,0,2,(ptPDE_->getS()).get());
  help2=ptPDE_->getAlgSys()->CalcEnergyNorm(0,0,5,(ptPDE_->getS1()).get());
  Double normsol=sqrt(help1+help2);

  if (normsol>maxnormsol_) maxnormsol_=normsol;

  Double relativeerror; 
  if (errorl2>1e-10) relativeerror=errorl2/maxnorml2sol_; /// !!!!! 
  else relativeerror=0.0;

//  if (relativeerror_>maxrelativeerror_) maxrelativeerror_=relativeerror_;
  maxrelativeerror_=relativeerror;

  std::cout << " ratio " << maxrelativeerror_ << std::endl;

  if (Calc3DerFromEquation_) {
                                  ptPDE_->CalcThirdDerivateFromEquation(thirddersoldt_);

  std::cout << " thirdder " << thirddersoldt_ << std::endl;

                                  Calc3DerFromEquation_=FALSE;
                             }
  else  thirddersoldt_=thirddersolAdt*2 - thirddersoldt_;

}
*/

void AcousticTimeErrorEstimator::ChangeStep(Double & dt)
{
#ifdef TRACE
  (*trace) << "entering AcousticTimeErrorEstimator::ChangeStep" << std::endl;
#endif


if (maxrelativeerror_<1e-10) dt*=3;
else
{
  Double help1=beta_*tol_/maxrelativeerror_;
  Double help=exp(0.333333333*log(help1));
  dt*=help;
}

std::cout << maxrelativeerror_ << " tol " << tol_ << std::endl;

  if (dt > maxdt_) dt=maxdt_;
  if (dt < mindt_) dt=mindt_;
}

Boolean AcousticTimeErrorEstimator::TestError(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering AcousticTimeErrorEstimator::TestError" << std::endl;
#endif

  CalcError(dt);

  if (InfoPrint)
  (*infofile) << " Ratio " << maxrelativeerror_ << " " << tol_ << std::endl;

  std::cout <<  " Ratio " << maxrelativeerror_ << " " << tol_ << std::endl;

  if (maxrelativeerror_>tol_) { counter_=0; return TRUE;} 
  else {
        if (maxrelativeerror_> beta_*tol_)
          { counter_=0; return FALSE;}
        else
          { counter_++; 
	    if (counter_== numrepeat_)
             {
               counter_=0;
	       return TRUE;
             }
            return FALSE;
          }
       }   
}

void AcousticTimeErrorEstimator::CalcThirdDer()
{
  Vector<Double> help=ptPDE_->getS2()-ptPDE_->getS2old();
  thirddersoldt_=help*2-thirddersoldt_;
}

} // end of namespace
