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
#include <DataInOut/WriteInfo.hh>

#include <PDE/basePDE.hh>

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

  if(isavebegin_ <= 0)
    Error("Value of stepsavebegin must be positive and nonzero! ",__FILE__,__LINE__);
   
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
  
  Integer level    = 0;
  Integer pdenumber= 0;
  Integer nsys     = 0;
  Double  steptime = firstdt_;
  Integer stepsave = isavebegin_;

  Double  dt = firstdt_;
  Boolean updatesysmat=FALSE;
  BasePDE * actPDE = ptdomain_->GetPDE(pdenumber);
  

  if (ptdomain_->GetNumPDE() <= 1) 
    {
      actPDE->InitTimeStepping(dt);
      
      ptdomain_->PrintGrid(level);

      if (PrintGridOnly)
	exit(0);
      
      actPDE->WriteGeneralPDEdefines();
      
      Integer nstep;
      for (nstep = 1; nstep <= numstep_; nstep++)
	{
	  Info->WriteTimeStep(actPDE->GetName(), nstep, steptime);

	  actPDE->PreStepTrans(nstep, steptime, level, updatesysmat);
	  actPDE->SolveStepTrans(nstep, steptime, level, updatesysmat);
	  actPDE->PostStepTrans(nstep,steptime,level);
	  
	  // writing results in output-file
	  if (nstep == stepsave && (nstep <= isaveend_))
	    { 
	      actPDE->PostProcess(level);
	      actPDE->WriteResultsInFile();
	      stepsave+=isaveincr_;
	    }
	  
	  steptime+=dt;	 
	}
    }
  else
    {
      BaseCoupledPDE * actCoupledPDE = ptdomain_->GetCoupledPDE();

      actCoupledPDE -> InitTimeStepping(dt);
      actCoupledPDE -> SetTimeSteppingParams(numstep_, firstdt_, isavebegin_, 
					     isaveend_, isaveincr_); 
      
      ptdomain_->PrintGrid(level);
      
      if (PrintGridOnly)
	exit(0);
      
      actCoupledPDE -> WriteGeneralPDEdefines();
      
      for (Integer nstep = 1; nstep <= numstep_; nstep++)
	{
	  Info->WriteTimeStep(actCoupledPDE->GetName(), nstep, steptime);

	  actCoupledPDE->InitStepTransCoupled(steptime);
	  
	  //  actCoupledPDE->PreStepTrans(nstep, steptime, level, updatesysmat);
	  actCoupledPDE->SolveStepTrans(nstep, steptime, level, updatesysmat);
	  // actCoupledPDE->PostStepTrans(level);
	  
	  // writing results in output-file
	  if (nstep == stepsave && (nstep <= isaveend_))
	    { 
	      actCoupledPDE->PostProcess(level);
	      actCoupledPDE->WriteResultsInFile();
	      stepsave+=isaveincr_;
	    }
	  steptime+=dt;
	}
    }
}

} // end of namespace
