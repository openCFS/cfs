#include <fstream>
#include <iostream>
#include <string>

#include "transientdriver.hh"

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
 
  ptdomain_->GetPDE(pdenumber)->SetupMatrices(ptdomain_->GetAlgSys(),type);

}

void TransientDriver :: SolveProblem()
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::SolveProblem" << std::endl;
#endif

  Integer level=0;
  Integer pdenumber  = 0;
  Integer nsys = 0;
  Double steptime=0.0;
  Integer stepsave=isavebegin_-1;

  std::cout << "steps" << isavebegin_ << " " << isaveend_ << " " << isaveincr_ << std::endl;

  Integer nstep;

  for (nstep = 0; nstep<numstep_; nstep++)
    {
      steptime += firstdt_;

//      ptdomain->GetPDE(pdenumber)->CalcParamForNewmarkMethod(firstdt_);

      ptdomain_->GetPDE(pdenumber)->SolveStepTrans(ptdomain_->GetAlgSys(), ptdomain_->GetBCs(), nstep, steptime, level);

   // writing results in output-file
    if (nstep == stepsave && (nstep < isaveend_))
      { 
        ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
        stepsave+=isaveincr_;
      }
    }
}

}
