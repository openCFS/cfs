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

HarmonicDriver :: HarmonicDriver(Domain * adomain)
:BaseDriver(adomain)
{
  ENTER_FCN( " HarmonicDriver::HarmonicDriver" );

#ifndef XMLPARAMS
  // get time steps information from conf-file

  conf->get("startFreq",startFreq_);
  conf->get("stopFreq",stopFreq_);
  conf->get("numFreq", numFreq_);

  saveType_ = 1;
  conf->ifget("saveType",saveType_);
  
#else
  Error("Frequency analysis in xml currently not supported");
#endif

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
  Integer pdenumber  = 0;
  ptdomain_->PrintGrid(level);

  if (PrintGridOnly)
      exit(0);

  BasePDE * actPDE = ptdomain_->GetPDE(pdenumber);
  actPDE->WriteGeneralPDEdefines();
      
  Integer fstep;
  Double actFreq  = startFreq_;
  Double freqIncr = (stopFreq_ - startFreq_) / numFreq_;

  for (fstep = 1; fstep <= numFreq_; fstep++) {
    Info->WriteHarmonicStep(actPDE->GetName(), fstep, actFreq);

    actPDE->PreStepHarmonic(fstep, actFreq, level, reset);
    actPDE->SolveStepHarmonic(fstep, actFreq, level, reset);
    actPDE->PostStepHarmonic(fstep, actFreq, level, reset);

    // writing results in output-file
    actPDE->PostProcess(level);
    actPDE->WriteResultsInFile();
    
    actFreq += freqIncr;
  }


}


} // end of namespace
