#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <math.h>

#include "transientdriver.hh"
#include <Utils/vector.hh>

#include <DataInOut/GMV/outGMV.hh>

#include <CoupledPDE/basecoupledpde.hh>
#include <CoupledPDE/itercoupledpde.hh>

#ifndef NEWBASEPDE
#include <PDE/basepde.hh>
#else
#include <PDE/newBasePDE.hh>
#endif //#ifndef NEWBASEPDE

namespace CoupledField
{

TransientDriver :: TransientDriver(Domain * adomain)
:BaseDriver(adomain)
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::TransientDriver" << std::endl;
#endif

  // get time steps information from conf-file
  conf->get("numsteps",numstep_);
  conf->get("firstdt", firstdt_);
  conf->get("stepsavebeg",isavebegin_);
  conf->get("stepsaveend",isaveend_);
  conf->get("stepsaveincr",isaveincr_);

  ptMeshes_=NULL;
}

TransientDriver :: ~TransientDriver()
{
#ifdef TRACE
  (*trace) << "entering TransientDriver::~TransientDriver" << std::endl;
#endif

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


  if (ptdomain_->GetNumPDE() <= 1) 
    {
      ptdomain_->GetPDE(pdenumber)->InitTimeStepping(dt);
      
      ptdomain_->PrintGrid(level);

      if (PrintGridOnly)
	exit(0);
      
      ptdomain_->GetPDE(pdenumber)->WriteGeneralPDEdefines();
      
      Integer nstep;
      for (nstep = 0; nstep<numstep_; nstep++)
	{
	  ptdomain_->GetPDE(pdenumber)->PreStepTrans(level, updatesysmat);
	  ptdomain_->GetPDE(pdenumber)->SolveStepTrans(nstep, steptime, level, updatesysmat);
	  ptdomain_->GetPDE(pdenumber)->PostStepTrans(level);
	  
	  // writing results in output-file
	  if (nstep == stepsave && (nstep < isaveend_))
	    { 
	      ptdomain_->GetPDE(pdenumber)->PostProcess(level);
	      ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
	      stepsave+=isaveincr_;
	    }
	  
	  steptime+=dt;
	}
    }
  else
    {
      ptdomain_->GetCoupledPDE() -> InitTimeStepping(dt);
      ptdomain_->GetCoupledPDE() -> SetTimeSteppingParams(numstep_, firstdt_, isavebegin_, 
							  isaveend_, isaveincr_); 
      
      ptdomain_->PrintGrid(level);
      
      if (PrintGridOnly)
	exit(0);
      
      ptdomain_->GetCoupledPDE() -> WriteGeneralPDEdefines();
      
      Integer nstep;
      for (nstep = 0; nstep<numstep_; nstep++)
	{
	  //	  ptdomain_->GetCoupledPDE()->PreStepTrans(level, updatesysmat);
	  ptdomain_->GetCoupledPDE()->SolveStepTrans(nstep, steptime, level, updatesysmat);
	  //	  ptdomain_->GetCoupledPDE()->PostStepTrans(level);
	  
	  // writing results in output-file
	  if (nstep == stepsave && (nstep < isaveend_))
	    { 
	      ptdomain_->GetCoupledPDE()->PostProcess(level);
	      ptdomain_->GetCoupledPDE()->WriteResultsInFile();
	      stepsave+=isaveincr_;
	    }
	  steptime+=dt;
	}
    }
}

} // end of namespace
