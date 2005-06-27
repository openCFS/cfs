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

  keyVec  = pdename, "sliceData", "soundSpeed";
  params->Get( keyVec, soundSpeed_);

  keyVec  = pdename, "sliceData", "waveLength";
  params->Get( keyVec, waveLength_);

  keyVec  = pdename, "sliceData", "nrSafetyElements";
  params->Get( keyVec, sik_);

  keyVec  = pdename, "sliceData", "addNrOfWaveLength";
  params->Get( keyVec, tol_);

  keyVec  = pdename, "sliceData", "elementsPerWavelength";
  params->Get( keyVec, elementsPerWavelength_);

  keyVec  = pdename, "sliceData", "dimension", "dimZ";
  params->Get( keyVec, dimZ_);


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
  Integer elemgrid;
  Double  timegrid, timefirst;
  Integer tol;
  Integer sik=sik_;//10;
  Integer nodeShift;//for function call SaveNodes
  Integer shiftFactor;// for function call SaveNodes
  Integer numShift=0; // for function call: SaveNodes ... number of shifts made
  Integer hel = TRUE;
  Double T,wavesPerPulse,way;
  Double  dt = firstdt_;
  Boolean updatesysmat=FALSE;
  Boolean mkdir=FALSE;
  Double meshsize = waveLength_/elementsPerWavelength_;

  tol = tol_;// 2*elementsPerWavelength_;
  //std::cout << " transient -Tol_ = " << tol << std::endl;
  T = waveLength_/soundSpeed_;
  //std::cout << " transient -T_ = " << T << std::endl;
  wavesPerPulse = pulseTime_/T;
  //std::cout << " transient -wavesperpulse_ = " << wavesPerPulse << std::endl;
  way =wavesPerPulse*elementsPerWavelength_;
  //std::cout << " transient -way_ = " << way<< std::endl;
  elemgrid= 2*(Integer)way+2*sik+ (Integer)tol;

  //std::cout << " transient -elemgrid = " << elemgrid << std::endl;
  //elemgrid =  2*((Integer) ((pulseTime_*soundSpeed_/waveLength_)*elementsPerWavelength_)) + 2*sik +tol;
  //std::cout << "transient- elemgrid = " << elemgrid << std::endl;
  timefirst = ((((elemgrid - sik)*1.0)/elementsPerWavelength_)*waveLength_)/soundSpeed_;
  
  //timegrid = ((((tol*1.0)/2)/elementsPerWavelength_)*waveLength_)/soundSpeed_ + pulseTime_;
  Integer numelemshift= elemgrid*1/2-sik;
  //std::cout << "transient- numelemshift = " << numelemshift << std::endl;
  timegrid = (numelemshift*meshsize)/soundSpeed_;
  //std::cout << "transient- timefirst = " << timefirst << std::endl;
  //std::cout << "transient- timegrid = " << timegrid << std::endl;

  //Integer numelemshift= (Integer) (((tol/2))*elementsPerWavelength_+((pulseTime_*soundSpeed_)/waveLength_)* elementsPerWavelength_);

  Integer testmaxelem = numelemshift;
  //std::cout << "testmaxelem init : " << testmaxelem << std::endl;

  //  Integer maxelem = (Integer) (dimZ_ * elementsPerWavelength_)/waveLength_;

  //std::cout << "dim = " << dimZ_ << std::endl;
  //std::cout << "elemper wave = " << elementsPerWavelength_ << std::endl;
  //std::cout << "wavelength = " << waveLength_ << std::endl;
  //std::cout << "maxelem = " << maxelem << std::endl;
/*
  //Just for testing reasons
  FileType * InFile_ = domain_-> GetInFile(); 
  TimeFunc *  zeitfunc_ = new TimeFunc(InFile_);
  zeitfunc_->SetStartTimeVector(4*4);

  //Testing the Transformation for the 3D-Case
  Integer Nodeshift;
  Grid * ptgrid_ = domain_->GetGrid();
  ptgrid_->TransformGridStruct(Nodeshift); 	
*/
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

  //already reached end of Domain?
  Boolean max = FALSE;
  
  Double timeshift = timegrid-timefirst;
  
  Integer nstep;
  Integer transformStep=0;
  Integer numTrans = 0;
  for (nstep = 1; nstep <= numstep_; nstep++) {
    
    if ( numstep_ < 50 )
      Info->WriteTimeStep(ptPDE_->GetName(), nstep+stepOffset_,
                          steptime+timeOffset_);
    
    
    if ( doTransform ) {
      std::cout << "Perform Transformation at timestep : " <<nstep << std::endl;

      testmaxelem += numelemshift;

      std::cout << "testmaxelem : " << testmaxelem << std::endl;
      std::cout << "maxelem     : " << dimZ_/(waveLength_/elementsPerWavelength_) << std::endl;

      if(testmaxelem > (dimZ_/(waveLength_/elementsPerWavelength_))) {
        max = TRUE;
        std::cout << "end of Domain reached -> No shift" << std::endl;
      } 
      else {
        std::cout << "shift" << std::endl;
      }
       
      ptPDE_->GetSolveStep()->TransformSol4Slice((UInt)nodeShift, (UInt)shiftFactor, (UInt)0);

      numShift++; //the number of transformations is incremented
      timeshift = 0.0;
      
      doTransform = FALSE;
      transformStep = nstep;
      
      numTrans++;
    }
    
    if (numTrans == 1 && transformStep == nstep) {
      updatesysmat = TRUE;
    }
    else {
      updatesysmat = FALSE;
    }

    //    ptPDE_->SetNumTransSlice(numTrans);
 
    ptPDE_->GetSolveStep()->SetActTime(steptime);
    ptPDE_->GetSolveStep()->SetActStep(nstep);
    ptPDE_->GetSolveStep()->PreStepTrans(updatesysmat);
    ptPDE_->GetSolveStep()->SolveStepTrans4Slice(updatesysmat);
    ptPDE_->GetSolveStep()->PostStepTrans();

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
//     if(numShift == 0 && hel){
//       ptPDE_->GetSolveStep()->TransformSol4Slice((UInt)nodeShift, (UInt)shiftFactor, (UInt)1);
//       hel = FALSE;
//     }

//     //The first time the function SaveNodes has been called, the directory has to be opened  
//     if(!mkdir){
//       ptPDE_->GetSolveStep()->SaveNodes(numstep_, meshsize, numShift, -1, elemgrid);
//       mkdir = TRUE;
//     }		
//     ptPDE_->GetSolveStep()->SaveNodes(shiftFactor, steptime, numShift, nodeShift, elemgrid);
    
    steptime+=dt;
    std::cout << "acttime=" <<  steptime << " timeshift=" << timeshift << " timegrid=" << timegrid << std::endl; 
    
    if (max==FALSE) {
      timeshift +=dt;
      if(timeshift>=timegrid) {
        doTransform = TRUE;
      }
    }
  }

}

} // end of namespace

