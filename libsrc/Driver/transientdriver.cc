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
  ptdomain->GetInFile()->ReadStepData(numstep,isavebegin,isaveincr,isaveend,firstdt);
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
 
  ptdomain->GetPDE(pdenumber)->SetupMatrices(ptdomain->GetAlgSys(),type);

}


void TransientDriver :: SolveProblem()
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::SolveProblem" << std::endl;
#endif
  Integer level=0;
  Integer pdenumber  = 0;
  Integer nsys = 0;

  Integer nstep;
  Double steptime=0.0;

  for (nstep = 0; nstep<numstep; nstep++)
    {
      steptime += firstdt;

//      ptdomain->GetPDE(pdenumber)->CalcParamForNewmarkMethod(firstdt);
      ptdomain->GetPDE(pdenumber)->SolveStepTrans(ptdomain->GetAlgSys(), ptdomain->GetBCs(), nstep, steptime, level);

      ptdomain->GetPDE(pdenumber)->WriteResultsInFile();
    }
}

}
