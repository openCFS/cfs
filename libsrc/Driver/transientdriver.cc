#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <math.h>

#include "transientdriver.hh"
#include "stdSolveStep.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "PDE/StdPDE.hh"
#include "Domain/domain.hh"

namespace CoupledField {


  // ===============
  //   Constructor
  // ===============
  TransientDriver::TransientDriver(Domain * adomain, 
                                   UInt stepOffset,
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
    if(isavebegin_ <= 0)  {
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
  
    Double  steptime  = firstdt_;
    UInt stepsave  = isavebegin_;
  
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
    ptPDE_->GetSolveStep()->SetTimeStep(dt);

    ptPDE_->GetSolveStep()->SetNumTimeSteps(numstep_);
  
    // if multiSequence is performed, the ms-driver
    // writes out the grid one time
    if (! isPartOfSequence_)
      ptdomain_->PrintGrid();
  
    ptPDE_->WriteGeneralPDEdefines();
  
    UInt nstep;
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
        // Output in steps of ten percent 
        if ((double)nstep >= percentCounter  ){           
          Info->WriteTimeStep(ptPDE_->GetName(), nstep+stepOffset_, 
                              steptime+timeOffset_);
          percentCounter += timeStepPercent;
        }
      }

      ptPDE_->GetSolveStep()->SetActTime(steptime);
      ptPDE_->GetSolveStep()->SetActStep(nstep);
      ptPDE_->GetSolveStep()->PreStepTrans(updatesysmat);
      ptPDE_->GetSolveStep()->SolveStepTrans(updatesysmat);
      ptPDE_->GetSolveStep()->PostStepTrans();
    
      // writing results in output-file
      if (nstep == stepsave && (nstep <= isaveend_)) { 
        ptPDE_->PostProcess();
        ptPDE_->WriteResultsInFile(nstep, steptime, stepOffset_, timeOffset_);
        stepsave+=isaveincr_;
      }
    
      steptime+=dt;        

      SETPROFILE("After Transient Step");
    }
  
  }

} // end of namespace
