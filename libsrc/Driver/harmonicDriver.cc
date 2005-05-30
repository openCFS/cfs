#include <fstream>
#include <iostream>
#include <string>

#include "harmonicDriver.hh"
#include "stdSolveStep.hh"
#include "assemble.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "PDE/StdPDE.hh"
#include "Domain/domain.hh"


namespace CoupledField
{

  HarmonicDriver :: HarmonicDriver(Domain * adomain, 
                                   UInt stepOffset,
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

    Boolean reset = TRUE;
  
    if (! isPartOfSequence_)
      ptdomain_->PrintGrid();

    // if driver is not part of multiSequence Driver, get list
    // of pdes which have to be solved and intialize them
    if (isPartOfSequence_ == FALSE){     
      GetMyPDEs();
      Info->StartProgress ("Starting to solve problem", FALSE);
    }
  
    ptPDE_->WriteGeneralPDEdefines();
      
    actFreq_  = startFreq_;

    UInt fstep;
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
      Info->WriteHarmonicStep(ptPDE_->GetName(), fstep, actFreq_);    
      ptPDE_->GetSolveStep()->PreStepHarmonic(fstep, actFreq_, reset);
      ptPDE_->GetSolveStep()->SolveStepHarmonic(fstep, actFreq_, reset);
      ptPDE_->GetSolveStep()->PostStepHarmonic(fstep, actFreq_, reset);
    
      // writing results in output-file
      ptPDE_->PostProcess();
      ptPDE_->WriteResultsInFile();
    
      actFreq_ += freqIncr;
    }
  
  }


} // end of namespace
