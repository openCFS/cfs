#include <fstream>
#include <iostream>
#include <string>

#include "transientdriver.hh"
#include "actimeerror.hh"

namespace CoupledField
{

TransientDriver :: TransientDriver(Domain<Point2D> * adomain)
:BaseDriver(adomain)
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::TransientDriver" << std::endl;
#endif

  // read step data from input file
  ptdomain_->GetInFile()->ReadStepData(numstep_,isavebegin_,isaveend_,isaveincr_,firstdt_);

}

TransientDriver :: ~TransientDriver()
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::~TransientDriver" << std::endl;
#endif

}


void TransientDriver :: SetupMatricesPDE(Integer pdenumber, Integer type)
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::SetUpMatricesPDE" << std::endl;
#endif
 
  ptdomain_->GetPDE(pdenumber)->SetupMatrices(type);

}

void TransientDriver :: SolveProblem()
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::SolveProblem" << std::endl;
#endif

  Integer level=0;
  Integer pdenumber  = 0;
  Integer nsys = 0;
  Double steptime=firstdt_;
  Integer stepsave=isavebegin_-1;

  Double dt=firstdt_;
  Boolean updatesysmat=FALSE;

  ptdomain_->GetPDE(pdenumber)->CalcParameters(dt);
  ptdomain_->GetPDE(pdenumber)->SetMatrixFactors();

  Integer nstep;
  for (nstep = 0; nstep<numstep_; nstep++)
    {
      ptdomain_->GetPDE(pdenumber)->SolveStepTrans(ptdomain_->GetBCs(), nstep, steptime, level, updatesysmat);

   // writing results in output-file
    if (nstep == stepsave && (nstep < isaveend_))
      { 
        ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
        stepsave+=isaveincr_;
      }

   steptime+=dt;
   }
}

void TransientDriver :: SolveProblemAdapt()
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::SolveProblemAdapt" << std::endl;
#endif

  Integer level=0;
  Integer pdenumber  = 0;
  Integer nsys = 0;
  Double steptime=firstdt_;
  Integer stepsave=isavebegin_-1;

//  TimeErrorEstimator ** ptTimeError=new TimeErrorEstimation * [pdenumber];
  TimeErrorEstimator * ptTimeError;
 
  ptTimeError=ptdomain_->GetPDE(pdenumber)->CreatePtTimeError();

  Double dt=firstdt_;
  Boolean resetsysmat=FALSE;

  ptdomain_->GetPDE(pdenumber)->CalcParameters(dt);
  ptdomain_->GetPDE(pdenumber)->SetMatrixFactors();

  Integer nstep;
  for (nstep = 0; nstep<numstep_; nstep++)
    {

      ptdomain_->GetPDE(pdenumber)->SolveStepTrans(ptdomain_->GetBCs(), nstep, steptime, level, resetsysmat);

   // writing results in output-file
    if (nstep == stepsave && (nstep < isaveend_))
      {
        ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
        stepsave+=isaveincr_;
      }

   if (ptTimeError->TestError(dt))
      {
         ptTimeError->ChangeStep(dt);

         ptdomain_->GetPDE(pdenumber)->CalcParameters(dt);
         ptdomain_->GetPDE(pdenumber)->SetMatrixFactors();

         resetsysmat=TRUE;
         std::cout << "We have change step" << stepsave << " " << dt << std::endl; 
      }

   steptime+=dt;
   }
}
}
