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
    
  adjustDamping_ = FALSE;
  std::string damp;
  keyVec = "harmonic", "adjustDamping";
  //   adjustDamping_ =  params->IsSet(keyVec, attrVec, valVec, adjustDamping_);
  adjustDamping_ = params->IsSet("adjustDamping",  "harmonic");
  //  if (damp == "yes" ) {
  //    adjustDamping_ = TRUE;
  //  }

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
  
  ptPDE_->WriteGeneralPDEdefines();
      
  Integer fstep;
  Double actFreq  = startFreq_;
  Double freqIncr;
  if (numFreq_ > 1) {
    freqIncr = (stopFreq_ - startFreq_) / (numFreq_-1);
  }
  else {
    freqIncr = stopFreq_ - startFreq_ ;
  }

  std::string errMsg;
  
  if ( adjustDamping_ ) {
    ptPDE_->getPDE_assemble()->SetStartFrequency(startFreq_);
  }

  for (fstep = 1; fstep <= numFreq_; fstep++) {
    Info->WriteHarmonicStep(ptPDE_->GetName(), fstep, actFreq);
    
    ptPDE_->GetSolveStep()->PreStepHarmonic(fstep, actFreq, level, reset);
    ptPDE_->GetSolveStep()->SolveStepHarmonic(fstep, actFreq, level, reset);
    ptPDE_->GetSolveStep()->PostStepHarmonic(fstep, actFreq, level, reset);
    
    // writing results in output-file
    ptPDE_->PostProcess(level);
    ptPDE_->WriteResultsInFile();
    
    actFreq += freqIncr;
  }
  
}


} // end of namespace
