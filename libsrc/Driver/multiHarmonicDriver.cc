#include <fstream>
#include <iostream>
#include <string>

#include "multiHarmonicDriver.hh"
#include "stdSolveStep.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Domain/domain.hh"
#include "PDE/StdPDE.hh"


#ifdef __sgi
#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#define POW pow
#else
#include <cstdarg>
#include <cstdio>
#include <cmath>
#define POW std::pow
#endif





namespace CoupledField
{

  // ========================================================================
  // ========================= piezoParamIdent - Part ===========================
  // ========================================================================

  // constructor
  // opens datafiles: measuredData.dat for input, imedCurve.dat and piezoLog.dat for output

  MultiHarmonicDriver::MultiHarmonicDriver(Domain * adomain,
                                           UInt stepOffset,
                                           Double timeOffset,
                                           std::string driverTag,
                                           Boolean isPartOfSequence)
    :SingleDriver(adomain, stepOffset, timeOffset, 
                  driverTag, isPartOfSequence){

    ENTER_FCN( "multiharmonic::multiharmonic" );

    //  ptDomain = adomain;

    StdVector<std::string> keyVec, attrVec, valVec;
    
    attrVec = "tag";
    valVec = driverTag_;
    
    // Get time stepping information from parameter object
    keyVec = "multiHarmonic", "startFreq";
    params->Get(keyVec, attrVec, valVec, startFreq_);
  
    keyVec = "multiHarmonic", "stopFreq";
    params->Get(keyVec, attrVec, valVec, stopFreq_);
    
    keyVec = "multiHarmonic", "numFreq";
    params->Get(keyVec, attrVec, valVec, numFreq_);

    keyVec = "multiHarmonic", "nrMultiHarmonics";
    params->Get(keyVec, attrVec, valVec, nrMultHarms_);

 
    Char* measuredData="measuredData.dat";
    allMeasuredData = new std::ifstream(measuredData, std::basic_ios<char>::in);

    if (!allMeasuredData)
      {
        std::cerr << "\n File measuredData.dat does not exist!" << std::endl;
        exit(1);
      }

    std::cout<<"\n Opens impedCurve.dat and piezoLog.dat ... "<<std::endl;

    std::string filename= "imped.dat";

    impedCurve = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);

    if (!impedCurve)
      {
        std::cerr << "\n ImpedanceCurve.dat could not be initialized" << std::endl;
      }

    std::string filenameLog= "piezoLog.dat";

    piezoLog = new std::ofstream(filenameLog.c_str(),std::basic_ios<char>::out);
    
    if (!piezoLog)
      {
        std::cerr << "\n piezoLog.dat could not be initialized" << std::endl;
      }

    std::string filenameParLog= "parLog.dat";
    parLog = new std::ofstream(filenameParLog.c_str(),std::basic_ios<char>::out);

    if (!parLog)
      {
        std::cerr << "\n piezoLog.dat could not be initialized" << std::endl;
      }

  } // end of constructor

  // destructor
  
  MultiHarmonicDriver::~MultiHarmonicDriver(){
    ENTER_FCN( "MultiHarmonicDriver::~MultiHarmonicDriver" );
    allMeasuredData->close();
    impedCurve->close();
    piezoLog->close();
    parLog->close();
  }

  void MultiHarmonicDriver::SolveProblem() {
    ENTER_FCN( "MultiHarmonicDriver::SolveProblem" );


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
      
    UInt fstep;
    Double actFreq  = startFreq_;
    Double freqIncr = (stopFreq_ - startFreq_) / numFreq_;
    std::string errMsg;
  
    // branch for single PDE
    for (fstep = 1; fstep <= numFreq_; fstep++) {
      Info->WriteHarmonicStep(ptPDE_->GetName(), fstep, actFreq);
    
      ptPDE_->GetSolveStep()->SetActFreq(actFreq);
      ptPDE_->GetSolveStep()->SetActStep(fstep);

      //      std::cout<<"\n multiHarm: 1 " <<std::endl;
      ptPDE_->GetSolveStep()->PreStepHarmonic(reset);
    
      //      std::cout<<"\n multiHarm: 2"  <<std::endl;
      ptPDE_->GetSolveStep()->SolveStepHarmonic(reset);
    
      //      std::cout<<"\n multiHarm: 3 " <<std::endl;
      ptPDE_->GetSolveStep()->PostStepHarmonic(reset);
    
      //      std::cout<<"\n multiHarm: nrMultHarms_ =  " <<nrMultHarms_ <<std::endl;
    
      // writing results in output-file
      ptPDE_->PostProcess();
      ptPDE_->WriteResultsInFile( fstep, actFreq);
    
      actFreq += freqIncr;
    }
  }



  void MultiHarmonicDriver::updateMaterialData(Vector<Double> & parameter, MaterialData * ptMaterial){
    ENTER_FCN( "multiharmonic::updateMaterialData");
  }


}// end of namespace coupled field
