#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <math.h>

#include "transientdriver.hh"
#include "Utils/vector.hh"

#include "DataInOut/GMV/outGMV.hh"

#include "CoupledPDE/itercoupledpde.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

#include <PDE/basePDE.hh>

namespace CoupledField {


// ===============
//   Constructor
// ===============
TransientDriver::TransientDriver(Domain * adomain, 
								 Integer stepOffset,
								 Double timeOffset, 
								 std::string driverTag,
								 Boolean isPartOfSequence) 
  : SingleDriver(adomain, stepOffset, timeOffset, 
				 driverTag, isPartOfSequence)
{
  ENTER_FCN( "TransientDriver::TransientDriver" );
  
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

  // Make consistency check. In fact in the XML case the Schema should catch
  // this error. But one can never be sure.
  if(isavebegin_ <= 0) 	{
	Error( "Value of stepsavebegin must be positive and nonzero!",
	       __FILE__, __LINE__ );
  }
  
  ptMeshes_ = NULL;
  
}

// ==============
//   Destructor
// ==============
TransientDriver::~TransientDriver()
{
  ENTER_FCN( "TransientDriver::~TransientDriver" );
}


// =================
//   Solve Problem
// =================
void TransientDriver::SolveProblem()
{
  ENTER_FCN( "TransientDriver::SolveProblem" );
  
  Integer level     = 0;
  Double  steptime  = firstdt_;
  Integer stepsave  = isavebegin_;
  
  Double  dt = firstdt_;
  Boolean updatesysmat=FALSE;
  
  // if driver is not part of multiSequence Driver, get list
  // of pdes which have to be solved and intialize them
  if (isPartOfSequence_ == FALSE) {     
    GetMyPDEs();
    Info->StartProgress ("Starting to solve problem", FALSE);
  }

  Double timeStepPercent = (double)numstep_/10;
  Double percentCounter = timeStepPercent;
  
  // Solve problem
  ptPDE_->SetTimeStep(dt);
  
  // if multiSequence is performed, the ms-driver
  // writes out the grid one time
  if (! isPartOfSequence_)
    ptdomain_->PrintGrid(level);
  
  if (PrintGridOnly) 
    exit(0);
  
  ptPDE_->WriteGeneralPDEdefines();
  
  Integer nstep;
  for (nstep = 1; nstep <= numstep_; nstep++) {
    
    if ( numstep_ <= 50 )
      Info->WriteTimeStep(ptPDE_->GetName(), nstep+stepOffset_, 
			  steptime+timeOffset_);
    else if ( (numstep_ > 50) && (numstep_ <= 500) ) {
      if ( (nstep%10) == 0 )		
	Info->WriteTimeStep(ptPDE_->GetName(), nstep+stepOffset_, 
			    steptime+timeOffset_);
    }
    else if ( numstep_ > 500 ) {
      // Output in steps of en percent 
      if ((double)nstep >= percentCounter  ){		
	Info->WriteTimeStep(ptPDE_->GetName(), nstep+stepOffset_, 
			    steptime+timeOffset_);
      percentCounter += timeStepPercent;
      }
    }
    
    ptPDE_->GetSolveStep()->PreStepTrans(nstep, steptime, level, updatesysmat);
    ptPDE_->GetSolveStep()->SolveStepTrans(nstep, steptime, level, updatesysmat);
    ptPDE_->GetSolveStep()->PostStepTrans(nstep,steptime,level);
    
    // writing results in output-file
    if (nstep == stepsave && (nstep <= isaveend_)) { 
      ptPDE_->PostProcess(level);
      ptPDE_->WriteResultsInFile(stepOffset_, timeOffset_);
      stepsave+=isaveincr_;
    }
    
    steptime+=dt;	 
  }
  
}

} // end of namespace
