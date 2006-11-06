#include "multiSequenceDriver.hh"

// include all kinds of single drivers
#include "staticdriver.hh"
#include "transientdriver.hh"
#include "harmonicDriver.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/pdememento.hh"
#include "Domain/domain.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  MultiSequenceDriver::MultiSequenceDriver( ) 
    : BaseDriver() {
    ENTER_FCN( "MultiSequenceDriver::MultiSequenceDriver" );

    analysis_ = MULTI_SEQUENCE;
    
    numSteps_ = 0;
    actStep_ = 0;
    actTime_ = 0;
    actDriver_ = 0;
  }


  // **********************
  //   Default destructor
  // **********************
  MultiSequenceDriver::~MultiSequenceDriver() {
    ENTER_FCN( "MultiSequenceDriver::~MultiSequenceDriver" );
  }


  // *****************
  //   Solve problem
  // *****************
  void MultiSequenceDriver::SolveProblem() {
    ENTER_FCN( "MultiSequenceDriver::SolveProblem" );

    // helper function for getting 
    // transient and harmonic parameters
    StdVector<std::string> keyVec, attrVec, valVec;
    
    // current driver
    //SingleDriver * actDriver_;

    // vector containing pointer to current set of PDEs
    StdVector<SinglePDE *> ptPDEs;
    
    // Vector of memento objects, which save the internal state
    // of a PDE
    StdVector<shared_ptr<PDEMemento> > memento;
    UInt nextStep = 0;
    Double nextTime = 0.0;
    UInt actNumSteps;
    Double actDt = 0.0;
  
    // Print out the grid
    domain->PrintGrid();

    Info->StartProgress("Starting to solve problem",false);

    // outer loop over all single sequences
    for (UInt iStep=0; iStep<numSteps_; iStep++) {
      
      // remeber current sequence step
      curSequenceStep_ = iStep;

      Info->WriteMultiSequenceStep(iStep+1, 
                                   analysisPerStep_[iStep][0]);

      // Since per time step only one type of 
      // analysis is allowed, we simple access
      // the first entry fo analysisPerStep_
      if (analysisPerStep_[iStep][0] == STATIC) {
        actDriver_ = new StaticDriver( actStep_, actTime_,
                                       tagsPerStep_[iStep][0], true );
        nextStep = actStep_ + 1;
      }
      else if (analysisPerStep_[iStep][0] == TRANSIENT) {      
        actDriver_ = new TransientDriver( actStep_, actTime_, 
                                          tagsPerStep_[iStep][0],true );

        // the time step and the current time have to adapted 
        // after solving the first multiSequence step
        attrVec = "tag";
        valVec = tagsPerStep_[iStep][0];

        keyVec = "transient", "numSteps"; 
        params->Get(keyVec, attrVec, valVec, actNumSteps);
          
        keyVec = "transient", "firstDt"; 
        params->Get(keyVec, attrVec, valVec, actDt);
        
        nextStep = actStep_ + actNumSteps;
        nextTime = actTime_ + actNumSteps * actDt;
      }
      else if (analysisPerStep_[iStep][0] == HARMONIC) {
        actDriver_ = new HarmonicDriver( actStep_, actTime_,
                                         tagsPerStep_[iStep][0], true );

        nextStep = actStep_ + 1;
      }
      

      // Initialize all PDEs
      domain->CreatePDEs(pdesPerStep_[iStep],iStep+1, tagsPerStep_[iStep]);
      
      domain->SetDriver( actDriver_ );
      actDriver_->SetPDE(domain->GetBasePDE());

      // Give the according pdes to the driver
      ptPDEs.Clear();
      ptPDEs.Resize(pdesPerStep_[iStep].GetSize());
      for (UInt iPde=0; iPde < pdesPerStep_[iStep].GetSize(); iPde++)
        ptPDEs[iPde]=domain->GetSinglePDE( pdesPerStep_[iStep][iPde] );

      // Initialize driver objects
      actDriver_->Init();

      // After the first run, initialize this PDE
      // with the solution of the previous run
      if (iStep > 0) {
        for (UInt iPde = 0; iPde < pdesPerStep_[iStep].GetSize(); iPde++ ) {
          if( memento[iPde] != NULL) {
            if (memento[iPde]->IsSet()) {
              ptPDEs[iPde]->SetMemento( memento[iPde], 
                                        valueUsagePerStep_[iStep][iPde] );
            }
          }
        }
      }

      //! Initialize Pdes, after having set the memento object
      domain->InitPDEs(iStep+1, tagsPerStep_[iStep]);
      
      // Solve Problem
      actDriver_->SolveProblem();

      // Get solution for next step and delete
      // all PDEs
      if (iStep < numSteps_-1) {
        memento.Resize(pdesPerStep_[iStep+1].GetSize());

        // Iterate over all PDEs in the next step
        for (UInt iPde = 0; iPde < pdesPerStep_[iStep+1].GetSize(); iPde++) {
          // Iterate over all PDEs in the current step
          for(UInt kPde = 0; kPde < pdesPerStep_[iStep].GetSize(); kPde++) {
            // If both match, then save the result of this step
            // for the next step
            if (pdesPerStep_[iStep+1][iPde] == pdesPerStep_[iStep][kPde])
              //dynamic_cast<const NodeStoreSol<Double>& >
              ptPDEs[kPde]->GetMemento(memento[iPde]);
          }
        }
        
        // delete PDEs
        domain->ResetPDEs();
      }


      // delete analysistypes
      delete actDriver_;
      
      // increase stepNumber and time
      actStep_ = nextStep;
      actTime_ = nextTime;
    } // iStep
  }


  UInt MultiSequenceDriver::GetActStep( const std::string& pdename ) {
    ENTER_FCN( "MultiSequenceDriver::GetActStep" );
    assert( actDriver_);
    return actDriver_->GetActStep( pdename );
  }

  // *****************
  //   Initializer
  // *****************
  void MultiSequenceDriver::Init() 
  {
    ENTER_FCN( "MultiSequenceDriver::Init" );
    
    StdVector<UInt> steps;
    StdVector<std::string> pdesAux, tagsAux, analysisAux, valueUsageAux;

    // 1.) Read in all steps
    params->GetList("index", steps, "multiSequence", "step");

    // 2.) ensure that all steps occur and none is left out
    UInt currStep = 0;
    for (UInt iStep=1; iStep<=steps.GetSize(); iStep++) {
      for (UInt j=0; j<steps.GetSize(); j++) {
        if (steps[j] == iStep) {
          currStep++;
          break;
        }
      }
    }
    
    std::string errMsg;
    if (currStep != steps.GetSize()) {
      errMsg = "MultiStepSequence: The indices for a multiSequence must start ";
      errMsg += "with 1 and have to be in consecutive order.\n";
      errMsg += "Please check your parameter file!";
      Error(errMsg.c_str(), __FILE__, __LINE__);
    }
    numSteps_ = steps.GetSize();
    
    // 3.) Resize 'outer' vectors
    pdesPerStep_.Resize(numSteps_);
    tagsPerStep_.Resize(numSteps_);
    analysisPerStep_.Resize(numSteps_);
    valueUsagePerStep_.Resize(numSteps_);
    
    AnalysisType analysisType;
    PDEMemento::ValueUsageType usage;

    // 4.) Read in all pdes, tags and simulation types

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    std::string stepString;

    attrVec = "", "index", "";
   

    for (UInt iStep=0; iStep<numSteps_; iStep++) {
      stepString = GenStr(iStep+1);
      //params->GetList("tag", tagsAux, "step", "pde");
      valVec  = "", stepString, "";

      keyVec  = "multiSequence", "step", "pde", "refTag";      
      params->GetList( keyVec, attrVec, valVec, tagsAux );

       keyVec  = "multiSequence", "step", "pde", "type";      
      params->GetList( keyVec, attrVec, valVec, pdesAux );    

      keyVec  = "multiSequence", "step", "pde", "analysis";      
      params->GetList( keyVec, attrVec, valVec, analysisAux );

      keyVec  = "multiSequence", "step", "pde", "usage";      
      params->GetList( keyVec, attrVec, valVec, valueUsageAux );

      //       std::cerr << "The pdes for step " << stepString << " are: " << std::endl;
      //       std::cerr << pdesAux << std::endl << std::endl;

      //       std::cerr << "The analysis for step " << stepString << " are: " << std::endl;
      //       std::cerr << analysisAux << std::endl << std::endl;

      // Add pdes, tags, and analysistypes of the current step
      // in the same order as defined by Domain
      for (UInt iPDE=0; iPDE<pdesAux.GetSize(); iPDE++) {
        pdesPerStep_[iStep].Push_back(pdesAux[iPDE]);
        String2Enum(analysisAux[iPDE], analysisType);
        
        tagsPerStep_[iStep].Push_back(tagsAux[iPDE]);
        analysisPerStep_[iStep].Push_back(analysisType);
        
        usage = PDEMemento::String2Enum( valueUsageAux[iPDE] );
        valueUsagePerStep_[iStep].Push_back( usage );
      }
      
    }
    
      
      
    // 5.) Perform final consistency checks
    
    // iterate over all steps
    bool pdeFound = false;
    bool tagFound = false;
    StdVector<std::string> tagStrings;
    StdVector<std::string> tagsLocal;

    // ********************************
    // Check if specified PDE types are 
    // defined in the 'pdeList' and 
    // occur only one time
    // ********************************

    for (UInt iStep=0; iStep<steps.GetSize(); iStep++) {
      
      // iterate over all pdes in current step
      for (UInt iPDE=0; iPDE<pdesPerStep_[iStep].GetSize(); iPDE++) {

        // ********************************
        // Check if specified PDE types are 
        // defined in the 'pdeList' and 
        // occur only one time
        // ********************************
        StdVector<std::string> pdeNames;
        params->GetPDEList( pdeNames );
        pdeFound = false;
        for (UInt kPDE=0; kPDE<pdeNames.GetSize(); kPDE++)
          if (pdesPerStep_[iStep][iPDE] == pdeNames[kPDE])
            pdeFound = true;

        if (!pdeFound) {
          errMsg = "MultiSequenceDriver::Init(): The PDE '";
          errMsg += pdesPerStep_[iStep][iPDE];
          errMsg += "' is not mentioned in the 'pdeList'!";
          Error(errMsg.c_str(), __FILE__, __LINE__);
        }

        // ************************************************
        // Check if the specified tags per pde and per step
        // occur in the 'bcsAndLoads'-section
        // ************************************************

        //get tags for current pde
        keyVec = "pdeList", pdesPerStep_[iStep][iPDE], "bcsAndLoads", "tag";
        attrVec = "", "", "";
        valVec = "", "", "";
        //      std::cerr << "pdesPerStep_[iStep][iPDE] = " << pdesPerStep_[iStep][iPDE]->GetName() << std::endl;
        params->GetList(keyVec, attrVec, valVec, tagsAux);

	tagFound = false;
        
        // Iterate over all 'bcsAndLoads' occurences for this PDE
        for (UInt iTag=0; iTag<tagsAux.GetSize(); iTag++) {

          // split the tag-string into a vector of 
          // seperate tags      
          SplitStringList(tagsAux[iTag], tagsLocal);
          
          // Iterate over all tags in one special 'bcsAndLoads' occurence
          for (UInt iTagLocal=0; iTagLocal<tagsLocal.GetSize(); iTagLocal++)
            if (tagsPerStep_[iStep][iPDE] == tagsLocal[iTagLocal]
                || tagsLocal[iTagLocal] == "anyTag" ) {
              tagFound = true;
              break;
            }
        }

	// NOTE: It can be okay to have no BcsAndLoads section
	//       e.g. acoustic-mechanic
	// Tag not found
        // if (!tagFound) {
        //   (*error) << "MultiSequenceDriver::Init(): The tag '"
	// 	   << tagsPerStep_[iStep][iPDE]
	// 	   << "' was not found in PDE '"
	// 	   << pdesPerStep_[iStep][iPDE] << "'!";
        //   Error( __FILE__, __LINE__ );
        // }
          
        // ************************************************
        // Check if one PDE occured more than 
        // one time per step and if for each step
        // only one type of analysistype occured
        // ************************************************
        for (UInt kPDE=iPDE+1; kPDE<pdesPerStep_[iStep].GetSize(); kPDE++) {
          if (pdesPerStep_[iStep][iPDE] == pdesPerStep_[iStep][kPDE]) {
            errMsg  = "MultiSequenceDriver::Init(): The PDE '";
            errMsg += pdesPerStep_[iStep][kPDE];
            errMsg += "' occured more than one time in step ";
            errMsg += GenStr(iStep+1);
            Error(errMsg.c_str(), __FILE__, __LINE__); 
          }
          if (analysisPerStep_[iStep][iPDE] != analysisPerStep_[iStep][kPDE]) {
            errMsg  = "MultiSequenceDriver::Init(): There were different "; 
            errMsg += "analysisypes defined for multi-Sequence step ";
            errMsg += GenStr(iStep+1);
            errMsg += ". Please correct parameter file!";
            Error(errMsg.c_str(), __FILE__, __LINE__); 

          }
        }
        
        // ************************************************
        // Check if each tag for analysistype matches an according
        // section (like e.g. <transient tag="myTag">   
        // ************************************************
        std::string analysis;
        
        // for static analysis we have no special section
        // to specify
        if (analysisPerStep_[iStep][iPDE] != STATIC) {

          Enum2String(analysisPerStep_[iStep][iPDE], analysis);

          keyVec = analysis, "tag";
          attrVec = "";
          valVec = "";
          params->GetList(keyVec, attrVec, valVec, tagStrings);

          tagFound = false;       
          for (UInt l=0; l<tagStrings.GetSize(); l++) {
            SplitStringList(tagStrings[l], tagsAux);
            for (UInt iTag=0; iTag<tagsAux.GetSize(); iTag++) 
              //            if (tagsAux[iTag] == "anyTag" ||
              //                tagsAux[iTag] == tagsPerStep_[iStep][iPDE]) {
              if (tagsAux[iTag] == tagsPerStep_[iStep][iPDE] ||
                  tagsAux[iTag] == "anyTag") {
                tagFound = true;
                break;
              }
          }
          
          if (! tagFound) {
            *error << "MultiSequenceDriver::Init(): The section for"
                   << "analysistype '" << analysis <<"' with the tag '"
                   << tagsPerStep_[iStep][iPDE] << "' was not found in step "
                   << GenStr(iStep+1);
            Error(__FILE__, __LINE__ );
          }
          
        }
      } // iPDE
    } // iStep
  } 
  
} // end of namespace
