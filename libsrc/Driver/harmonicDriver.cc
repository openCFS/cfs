#include <fstream>
#include <iostream>
#include <string>

#include "harmonicDriver.hh"

#include "DataInOut/GMV/outGMV.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

#include <PDE/basePDE.hh>


namespace CoupledField
{

HarmonicDriver :: HarmonicDriver(Domain * adomain, 
				 Integer stepOffset,
				 Double timeOffset,
				 std::string driverTag,
				 Boolean isPartOfSequence)
:SingleDriver(adomain, stepOffset, timeOffset, 
	      driverTag, isPartOfSequence)
{
  ENTER_FCN( " HarmonicDriver::HarmonicDriver" );

  // vectors for accessing parameters
  StdVector<std::string> keyVec, attrVec, valVec;
  
  attrVec = "tag";
  valVec = driverTag_;
  
  // Get time stepping information from parameter object
  keyVec = "harmonic", "startFreq";
  params->Get(keyVec, attrVec, valVec, startFreq_);
  
  keyVec = "harmonic", "stopFreq";
  params->Get(keyVec, attrVec, valVec, stopFreq_);
  
  keyVec = "harmonic", "numFreq";
  params->Get(keyVec, attrVec, valVec, numFreq_);
    
}

HarmonicDriver :: ~HarmonicDriver()
{
  ENTER_FCN( " HarmonicDriver::~HarmonicDriver" );

}

void HarmonicDriver :: SolveProblem()
{
  ENTER_FCN( " HarmonicDriver::SolveProblem" );

  Integer level=0;
  Boolean reset = TRUE;
  
  if (! isPartOfSequence_)
  ptdomain_->PrintGrid(level);

  if (PrintGridOnly)
      exit(0);

  // if driver is not part of multiSequence Driver, get list
  // of pdes which have to be solved and intialize them
  if (isPartOfSequence_ == FALSE){     
    GetMyPDEs();
    Info->StartProgress ("Starting to solve problem", FALSE);
  }
  
  pdes_[0]->WriteGeneralPDEdefines();
      
  Integer fstep;
  Double actFreq  = startFreq_;
  Double freqIncr = (stopFreq_ - startFreq_) / numFreq_;
  std::string errMsg;
  
  if (pdes_.GetSize() <= 1) {

    // branch for single PDE
    for (fstep = 1; fstep <= numFreq_; fstep++) {
      Info->WriteHarmonicStep(pdes_[0]->GetName(), fstep, actFreq);
      
      pdes_[0]->GetSolveStep()->PreStepHarmonic(fstep, actFreq, level, reset);
      pdes_[0]->GetSolveStep()->SolveStepHarmonic(fstep, actFreq, level, reset);
      pdes_[0]->GetSolveStep()->PostStepHarmonic(fstep, actFreq, level, reset);
           
      // writing results in output-file
      pdes_[0]->PostProcess(level);
      pdes_[0]->WriteResultsInFile();
           
      actFreq += freqIncr;
    }
  } 
  else {
    errMsg  = "HarmonicDriver::Solve: Harmonic simulation is ";
    errMsg += "not yet implemented for coupled PDEs, sorry!";
    Error(errMsg.c_str(), __FILE__, __LINE__);
  }

}


} // end of namespace
