#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <math.h>

#include "transientdriver.hh"
#include "Utils/vector.hh"

#include "DataInOut/GMV/outGMV.hh"

#include "CoupledPDE/basecoupledpde.hh"
#include "CoupledPDE/itercoupledpde.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

#include <PDE/basePDE.hh>

namespace CoupledField
{


  // ===============
  //   Constructor
  // ===============
  TransientDriver :: TransientDriver(Domain * adomain, 
				     Integer stepOffset,
				     Double timeOffset, 
				     std::string driverTag,
				     Boolean isPartOfSequence) 
    : SingleDriver(adomain, stepOffset, timeOffset, 
		   driverTag, isPartOfSequence)
  {
    ENTER_FCN( "TransientDriver::TransientDriver" );

#ifndef XMLPARAMS
    // get time steps information from conf-file
    conf->get("numsteps",numstep_);
    conf->get("firstdt", firstdt_);
    conf->get("stepsavebeg",isavebegin_);
    conf->get("stepsaveend",isaveend_);
    conf->get("stepsaveincr",isaveincr_);

#else

    // vecotrs for accessing parameters
    StdVector<std::string> keyVec, attrVec, valVec;

    attrVec = "tag";
    valVec = driverTag_;
   

    // Get time stepping information from parameter object
    keyVec = "transient", "numSteps";
    params->Get(keyVec, attrVec, valVec, numstep_);
    
    keyVec = "transient", "firstDt";
    params->Get(keyVec, attrVec, valVec, firstdt_);

    keyVec = "transient", "stepSaveBeg";
    params->Get(keyVec, attrVec, valVec, isavebegin_);

    keyVec = "transient", "stepSaveEnd";
    params->Get(keyVec, attrVec, valVec, isaveend_);

    keyVec = "transient", "stepSaveInc";
    params->Get(keyVec, attrVec, valVec, isaveincr_);  


#endif

    // Make consistency check. In fact in the XML case the Schema should catch
    // this error. But one can never be sure.
    if(isavebegin_ <= 0)
      {
	Error( "Value of stepsavebegin must be positive and nonzero!",
	       __FILE__, __LINE__ );
      }
   
    ptMeshes_ = NULL;

  }


  // ==============
  //   Destructor
  // ==============
  TransientDriver :: ~TransientDriver()
  {
    ENTER_FCN( "TransientDriver::~TransientDriver" );
  }


  // =================
  //   Solve Problem
  // =================
  void TransientDriver :: SolveProblem()
  {
    ENTER_FCN( "TransientDriver::SolveProblem" );
  
    Integer level     = 0;
    Integer pdenumber = 0;
    Integer nsys      = 0;
    Double  steptime  = firstdt_;
    Integer stepsave  = isavebegin_;

    Double  dt = firstdt_;
    Boolean updatesysmat=FALSE;

    // if driver is not part of multiSequence Driver, get list
    // of pdes which have to be solved and intialize them
    if (isPartOfSequence_ == FALSE){     
      GetMyPDEs();
      Info->StartProgress ("Starting to solve problem", FALSE);
    }
    
    // Solve problem
    if (pdes_.GetSize() <= 1) 
      {
	
	// branch for single PDE
	pdes_[0]->SetTimeStep(dt);

	// if multiSequence is performed, the ms-driver
      // writes out the grid one time
	if (! isPartOfSequence_)
	  ptdomain_->PrintGrid(level);

	if (PrintGridOnly) exit(0);
      
	pdes_[0]->WriteGeneralPDEdefines();
      
	Integer nstep;
	for (nstep = 1; nstep <= numstep_; nstep++)
	  {
	    Info->WriteTimeStep(pdes_[0]->GetName(), nstep+stepOffset_, 
				steptime+timeOffset_);

	    pdes_[0]->PreStepTrans(nstep, steptime, level, updatesysmat);
	    pdes_[0]->SolveStepTrans(nstep, steptime, level, updatesysmat);
	    pdes_[0]->PostStepTrans(nstep,steptime,level);
	  
	    // writing results in output-file
	    if (nstep == stepsave && (nstep <= isaveend_))
	      { 
		pdes_[0]->PostProcess(level);
		pdes_[0]->WriteResultsInFile(stepOffset_, timeOffset_);
		stepsave+=isaveincr_;
	      }
	  
	    steptime+=dt;	 
	  }
      }
    else
      {
	BaseCoupledPDE * actCoupledPDE = ptdomain_->GetCoupledPDE();
	
	actCoupledPDE->SetTimeStep(firstdt_);

	// if multiSequence is performed, the ms-driver
	// writes out the grid one time
	if (! isPartOfSequence_)
	  ptdomain_->PrintGrid(level);
      
	if (PrintGridOnly) exit(0);
      
	actCoupledPDE -> WriteGeneralPDEdefines();

	// define which PDEs participate in solving process
	ptdomain_->GetCoupledPDE()->DefineSolvingPDEs(pdes_);

	for (Integer nstep = 1; nstep <= numstep_; nstep++)
	  {
	    Info->WriteTimeStep(actCoupledPDE->GetName(), nstep, steptime);

	    actCoupledPDE->InitStepTransCoupled(steptime);
	  
	    // actCoupledPDE->PreStepTrans(nstep,steptime,level,updatesysmat);
	    actCoupledPDE->SolveStepTrans(nstep, steptime, level,updatesysmat);
	    // actCoupledPDE->PostStepTrans(level);
	    
	    // writing results in output-file
	    if (nstep == stepsave && (nstep <= isaveend_))
	      { 
		actCoupledPDE->PostProcess(level);
		actCoupledPDE->WriteResultsInFile(stepOffset_, timeOffset_);
		stepsave+=isaveincr_;
	      }
	    steptime+=dt;
	  }
      }
  }

} // end of namespace
