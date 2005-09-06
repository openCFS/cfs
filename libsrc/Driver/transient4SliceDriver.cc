#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include <math.h>

#include "transient4SliceDriver.hh"
#include "stdSolveStep.hh"
#include "PDE/StdPDE.hh"
#include "Utils/vector.hh"

#include "DataInOut/GMV/outGMV.hh"
#include "Domain/domain.hh"
//#include "CoupledPDE/basecoupledpde.hh"
//#include "CoupledPDE/itercoupledpde.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
//For testing reasons
#include "Domain/domain.hh"

namespace CoupledField {


// ===============
//   Constructor
// ===============
Transient4SliceDriver::Transient4SliceDriver(Domain * adomain, 
								 Integer stepOffset,
								 Double timeOffset,
								 std::string driverTag,
								 Boolean isPartOfSequence)
  : SingleDriver(adomain, stepOffset, timeOffset, 
				 driverTag, isPartOfSequence)
{
  ENTER_FCN( "Transient4SliceDriver::Transient4SliceDriver" );
  //either for testing reasons only
  domain_ =  adomain;  
  // vecotrs for accessing parameters
  StdVector<std::string> keyVec, attrVec, valVec;

  attrVec = "tag";
  valVec = driverTag_;
   

  // Get time stepping information from parameter object
  keyVec = "transient4Slice", "numSteps";
  params->Get(keyVec, attrVec, valVec, numstep_);
  
  keyVec = "transient4Slice", "firstDt";
  params->Get(keyVec, attrVec, valVec, firstdt_);
  
  keyVec = "transient4Slice", "stepSaveBeg";
  params->Get(keyVec, attrVec, valVec, isavebegin_);

  keyVec = "transient4Slice", "stepSaveEnd";
  params->Get(keyVec, attrVec, valVec, isaveend_);

  keyVec = "transient4Slice", "stepSaveInc";
  params->Get(keyVec, attrVec, valVec, isaveincr_);

  //Get variables for calculating time to shift grid
  std::string pdename = "acoustic";

  keyVec  = pdename, "sliceData", "pulseTime";
  params->Get( keyVec, pulseTime_);

  keyVec  = pdename, "sliceData", "pulseOffset";
  params->Get( keyVec, pulseOffset_);

  keyVec  = pdename, "sliceData", "soundSpeed";
  params->Get( keyVec, soundSpeed_);

  keyVec  = pdename, "sliceData", "waveLength";
  params->Get( keyVec, waveLength_);

  keyVec  = pdename, "sliceData", "safetyRegion";
  params->Get( keyVec, safetyRegion_);

  keyVec  = pdename, "sliceData", "addNrOfWaveLength";
  params->Get( keyVec, tol_);

  keyVec  = pdename, "sliceData", "elementsPerWavelength";
  params->Get( keyVec, elementsPerWavelength_);

  keyVec  = pdename, "sliceData", "dimension", "dimZ";
  params->Get( keyVec, dimZ_);

  if ( safetyRegion_ < 0 ) {
    Error("safetyRegion value has to be positive",__FILE__,__LINE__);
  }

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
Transient4SliceDriver::~Transient4SliceDriver()
{
  ENTER_FCN( "Transient4SliceDriver::~Transient4SliceDriver" );
}


// =================
//   Solve Problem
// =================
void Transient4SliceDriver::SolveProblem()
{
  ENTER_FCN( "Transient4SliceDriver::SolveProblem" );

  Double  steptime  = firstdt_;
  Integer stepsave  = isavebegin_;

  Double  timeFirst, timeOffset, timeShift;

  //All the variables I need for saving the nodes
  //****************************************************************************
  UInt nodeShift = 0;//for function call SaveNodes
  UInt shiftFactor=0;// for function call SaveNodes
  UInt numShift=0; // for function call: SaveNodes ... number of shifts made
  UInt elemgrid = 0;
  Double  meshsize = 0;
  //****************************************************************************


  Double  dt = firstdt_;
  Boolean updatesysmat=FALSE;
  Boolean mkdir=FALSE;
  
  // timeFirst is the time, the pulse needs at the first time, he runs in the slice
  timeFirst   = 2*pulseTime_  + (tol_ * waveLength_) / soundSpeed_;
  timeFirst  -= (safetyRegion_ * waveLength_) / soundSpeed_;
  timeFirst  += pulseOffset_; 
  // timeOffset is the time, the pulse needs to run thru a shifted slice
  timeOffset  = pulseTime_  + (tol_ * waveLength_) / soundSpeed_;
  timeOffset -= (safetyRegion_ * waveLength_) / soundSpeed_;

  // if driver is not part of multiSequence Driver, get list
  // of pdes which have to be solved and intialize them

  if (isPartOfSequence_ == FALSE) {
    GetMyPDEs();
    Info->StartProgress ("Starting to solve the problem", FALSE);
  }

  // Set the time step value
  ptPDE_->GetSolveStep()->SetTimeStep(dt);
  

  // if multiSequence is performed, the ms-driver
  // writes out the grid one time
  if (! isPartOfSequence_)
    ptdomain_->PrintGrid();
  
  ptPDE_->WriteGeneralPDEdefines();
  
  //time for transformation reached?
  Boolean doTransform = FALSE;

  timeShift = timeFirst;
  
  Integer nstep;
  for (nstep = 1; nstep <= numstep_; nstep++) {
    
    if ( numstep_ < 50 )
      Info->WriteTimeStep(ptPDE_->GetName(), nstep+stepOffset_,
                          steptime+timeOffset_);
    
    
    if ( doTransform ) {
      std::cout << "Perform Transformation at timestep : " <<nstep << std::endl;

      ptPDE_->GetSolveStep()->TransformSol4Slice(shiftFactor,
		    nodeShift, elemgrid, meshsize,(UInt) 0);

      doTransform = FALSE;
   
      numShift++;
    }
    
    if (numShift == 1) {
      updatesysmat = TRUE;
    }
    else {
      updatesysmat = FALSE;
    }
 
    ptPDE_->GetSolveStep()->SetActTime(steptime);
    ptPDE_->GetSolveStep()->SetActStep(nstep);
    ptPDE_->GetSolveStep()->PreStepTrans(updatesysmat);
    ptPDE_->GetSolveStep()->SolveStepTrans4Slice(updatesysmat);
    ptPDE_->GetSolveStep()->PostStepTrans();

    //write history data
    ptPDE_->WriteHistoryInFile(nstep, steptime, stepOffset_, timeOffset_);

    // writing results in output-file
    if (nstep == stepsave && (nstep <= isaveend_)) {
      ptPDE_->PostProcess();
      ptPDE_->WriteResultsInFile(nstep, steptime, stepOffset_, timeOffset_);
      stepsave+=isaveincr_;
    }
    /*============================================================================================
     *
     * Filtering Section: - Calculating the Global Node Number
     *		      - Pick up the Nod'z of interest
     *                    - Save for each node, the pressure Value for the actual Time-Step
     *===========================================================================================*/
    //SaveNodes is implemented in StdPDE
    if ( nstep == 1 ){
      //ptPDE_->GetSolveStep()->TransformSol4Slice((UInt)nodeShift, (UInt)shiftFactor, (UInt)1);
      ptPDE_->GetSolveStep()->TransformSol4Slice(shiftFactor,
		nodeShift, elemgrid, meshsize, 1);  
    }
     
    //The first time the function SaveNodes has been called, the directory has to be opened  
    if(!mkdir){
      /*
	The First time the function is called, it need,
	numstep_ = referring to the number of steps, that have to be made
	meshsize = referring to the meshsize of the grid
      */
      
      ptPDE_->GetSolveStep()->SaveNodes(numstep_, meshsize, numShift, -1, elemgrid);
      mkdir = TRUE;
    }		
    ptPDE_->GetSolveStep()->SaveNodes(shiftFactor, steptime, numShift, nodeShift, elemgrid);
    

    steptime+=dt;

    std::cout << "acttime=" <<  steptime << " timeshift=" << timeShift << std::endl; 
    
    if ( steptime >= timeShift ) {
        doTransform = TRUE;
	timeShift  += timeOffset;
    }
    
  }

}

} // end of namespace

