#include "multiSequenceDriver.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "Utils/tools.hh"
#include "PDE/basePDE.hh"
#include "staticdriver.hh"
#include "transientdriver.hh"
#include "harmonicDriver.hh"
#include "Utils/nodestoresol.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  MultiSequenceDriver::MultiSequenceDriver(Domain * adomain) 
    : BaseDriver(adomain) {
    ENTER_FCN( "MultiSequenceDriver::MultiSequenceDriver" );

    numSteps_ = 0;
    actStep_ = 0;
    actTime_ = 0;
    
    // Not much to do...
    Init();
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
    SingleDriver * actDriver;

    // vector containing pointer to current set of PDEs
    StdVector<StdPDE *> ptPDEs;
    
    // Vector of memento objects, which save the internal state
    // of a PDE
    StdVector<PDEMemento> memento;
    Integer level = 0;
    Integer nextStep = 0;
    Double nextTime = 0.0;
    Integer actNumSteps;
    Double actDt = 0.0;
    Double actFrequency;
  
    // helper variables
    Integer iPDE, kPDE;

    // Print out the grid
    ptdomain_->PrintGrid(level);

    Info->StartProgress("Starting to solve problem",FALSE);

    // outer loop over all single sequences
    for (Integer iStep=0; iStep<numSteps_; iStep++) {
      
      Info->WriteMultiSequenceStep(iStep+1, 
				   analysisPerStep_[iStep][0]);

      // Since per time step only one type of 
      // analysis is allowed, we simple access
      // the first entry fo analysisPerStep_
      if (analysisPerStep_[iStep][0] == STATIC) {
	actDriver = new StaticDriver(ptdomain_, actStep_, actTime_,
				     tagsPerStep_[iStep][0], TRUE);
	nextStep = actStep_ + 1;
      }
      else if (analysisPerStep_[iStep][0] == TRANSIENT) {      
	actDriver = new TransientDriver(ptdomain_, actStep_, actTime_, 
					tagsPerStep_[iStep][0],TRUE);

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
	actDriver = new HarmonicDriver(ptdomain_, actStep_, actTime_,
				       tagsPerStep_[iStep][0], TRUE);
      }
      
      // Initialize all PDEs
      ptdomain_->InitPDEs(pdesPerStep_[iStep],iStep+1, tagsPerStep_[iStep]);
      actDriver->SetPDE(ptdomain_->GetBasePDE());

      // Give the according pdes to the driver
      ptPDEs.Clear();
      ptPDEs.Resize(pdesPerStep_[iStep].GetSize());
      for (iPDE=0; iPDE<pdesPerStep_[iStep].GetSize(); iPDE++)
	ptPDEs[iPDE]=ptdomain_->GetStdPDE(pdesPerStep_[iStep][iPDE]);



      // After the first run, initialize this PDE
      // with the solution of the previous run
      if (iStep > 0) {
	std::string transFromTo = "standard";
	//check, if hamonic analysis is followed by transient one
	if ( analysisPerStep_[iStep-1][0] == HARMONIC && 
	     analysisPerStep_[iStep][0] == TRANSIENT) {
	  transFromTo = "complexToReal";
	}
	for (Integer i=0; i<pdesPerStep_[iStep].GetSize(); i++) {
	  if (memento[i].IsSet()) {
	    ptPDEs[i]->SetMemento( memento[i], transFromTo, actFrequency );
	  }
	}
      }
      
      // Solve Problem
      actDriver->SolveProblem();

      //if harmonic analysis, save the used frequencies
      if ( analysisPerStep_[iStep][0] == HARMONIC ) {
	actFrequency = actDriver->GetActFrequency();
      }
      else {
	actFrequency = 0;
      }


      // Get solution for next step and delete
      // all PDEs
      if (iStep < numSteps_-1) {
	memento.Resize(pdesPerStep_[iStep+1].GetSize());


	// Iterate over all PDEs in the next step
	for (iPDE=0; iPDE<pdesPerStep_[iStep+1].GetSize(); iPDE++) {
	  // Iterate over all PDEs in the current step
	  for(kPDE=0; kPDE<pdesPerStep_[iStep].GetSize(); kPDE++) {
	    // If both match, then save the result of this step
	    // for the next step
	    if (pdesPerStep_[iStep+1][iPDE] == pdesPerStep_[iStep][kPDE])
	      //dynamic_cast<const NodeStoreSol<Double>& >
	      ptPDEs[kPDE]->GetMemento(memento[iPDE]);
	  }
	}
	
	// delete PDEs
	ptdomain_->ResetPDEs();
      }


      // delete analysistypes
      delete actDriver;
      
      // increase stepNumber and time
      actStep_ = nextStep;
      actTime_ = nextTime;
    } // iStep
  }


  // *****************
  //   Initializer
  // *****************
  void MultiSequenceDriver::Init() 
  {
    ENTER_FCN( "MultiSequenceDriver::Init" );
    
    StdVector<Integer> steps;
    StdVector<std::string> pdesAux, tagsAux, analysisAux;

    // 1.) Read in all steps
    params->GetList("index", steps, "multiSequence", "step");

    // 2.) ensure that all steps occur and none is left out
    Integer currStep = 0;
    for (Integer iStep=1; iStep<=steps.GetSize(); iStep++) {
      for (Integer j=0; j<steps.GetSize(); j++) {
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
      errMsg += "Please check your paramter file!";
      Error(errMsg.c_str(), __FILE__, __LINE__);
    }
    numSteps_ = steps.GetSize();
    
    // 3.) Resize 'outer' vectors
    pdesPerStep_.Resize(numSteps_);
    tagsPerStep_.Resize(numSteps_);
    analysisPerStep_.Resize(numSteps_);
    
    AnalysisType analysisType;

    // 4.) Read in all pdes, tags and simulation types

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    std::string stepString;

    attrVec = "", "index", "";
   

    for (Integer iStep=0; iStep<numSteps_; iStep++) {
      stepString = Info->GenStr(iStep+1);
      //params->GetList("tag", tagsAux, "step", "pde");
      valVec  = "", stepString, "";

      keyVec  = "multiSequence", "step", "pde", "refTag";      
      params->GetList( keyVec, attrVec, valVec, tagsAux );

      //       std::cerr << "The tags for step " << stepString << " are: " << std::endl;
      //       std::cerr << tagsAux << std::endl << std::endl;
      

      keyVec  = "multiSequence", "step", "pde", "type";      
      params->GetList( keyVec, attrVec, valVec, pdesAux );    

      keyVec  = "multiSequence", "step", "pde", "analysis";      
      params->GetList( keyVec, attrVec, valVec, analysisAux );

      //       std::cerr << "The pdes for step " << stepString << " are: " << std::endl;
      //       std::cerr << pdesAux << std::endl << std::endl;

      //       std::cerr << "The analysis for step " << stepString << " are: " << std::endl;
      //       std::cerr << analysisAux << std::endl << std::endl;

      // Add pdes, tags, and analysistypes of the current step
      // in the same order as defined by Domain
      for (Integer iPDE=0; iPDE<pdesAux.GetSize(); iPDE++) {
	pdesPerStep_[iStep].Push_back(pdesAux[iPDE]);
	String2Enum(analysisAux[iPDE], analysisType);
	
	tagsPerStep_[iStep].Push_back(tagsAux[iPDE]);
	analysisPerStep_[iStep].Push_back(analysisType);
      }
      
    }
    
      
      
      // 5.) Perform final consistency checks
    
    // iterate over all steps
    Boolean pdeFound = FALSE;
    Boolean tagFound = FALSE;
    StdVector<std::string> tagStrings;
    StdVector<std::string> tagsLocal;

    // ********************************
    // Check if specified PDE types are 
    // defined in the 'pdeList' and 
    // occur only one time
    // ********************************

    for (Integer iStep=0; iStep<steps.GetSize(); iStep++) {
      
      // iterate over all pdes in current step
      for (Integer iPDE=0; iPDE<pdesPerStep_[iStep].GetSize(); iPDE++) {

	// ********************************
	// Check if specified PDE types are 
	// defined in the 'pdeList' and 
	// occur only one time
	// ********************************
	StdVector<std::string> pdeNames;
	params->GetPDEList( pdeNames );
	pdeFound = FALSE;
	for (Integer kPDE=0; kPDE<pdeNames.GetSize(); kPDE++)
	  if (pdesPerStep_[iStep][iPDE] == pdeNames[kPDE])
	    pdeFound = TRUE;

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
// 	std::cerr << "pdesPerStep_[iStep][iPDE] = " << pdesPerStep_[iStep][iPDE]->GetName() << std::endl;
	params->GetList(keyVec, attrVec, valVec, tagsAux);

	tagFound = FALSE;
	
	// Iterate over all 'bcsAndLoads' occurences for this PDE
	for (Integer iTag=0; iTag<tagsAux.GetSize(); iTag++) {

	  // split the tag-string into a vector of 
	  // seperate tags	
	  SplitStringList(tagsAux[iTag], tagsLocal);
	  
	  // Iterate over all tags in one special 'bcsAndLoads' occurence
	  for (Integer iTagLocal=0; iTagLocal<tagsLocal.GetSize(); iTagLocal++)
	    if (tagsPerStep_[iStep][iPDE] == tagsLocal[iTagLocal]
 		|| tagsLocal[iTagLocal] == "anyTag" ) {
	    tagFound = TRUE;
	    break;
	    }
	}
	
	if (!tagFound) {
	  errMsg  = "MultiSequenceDriver::Init(): The tag '";
	  errMsg += tagsPerStep_[iStep][iPDE];
	  errMsg += "' was not found in PDE '";
	  errMsg += pdesPerStep_[iStep][iPDE];
	  errMsg += "'!";
	  Error(errMsg.c_str(), __FILE__, __LINE__);
	}
	  
	// ************************************************
	// Check if one PDE occured more than 
	// one time per step and if for each step
	// only one type of analysistype occured
	// ************************************************
	for (Integer kPDE=iPDE+1; kPDE<pdesPerStep_[iStep].GetSize(); kPDE++) {
	  if (pdesPerStep_[iStep][iPDE] == pdesPerStep_[iStep][kPDE]) {
	    errMsg  = "MultiSequenceDriver::Init(): The PDE '";
	    errMsg += pdesPerStep_[iStep][kPDE];
	    errMsg += "' occured more than one time in step ";
	    errMsg += Info->GenStr(iStep+1);
	    Error(errMsg.c_str(), __FILE__, __LINE__); 
	  }
	  if (analysisPerStep_[iStep][iPDE] != analysisPerStep_[iStep][kPDE]) {
	    errMsg  = "MultiSequenceDriver::Init(): There were different "; 
	    errMsg += "analysisypes defined for multi-Sequence step ";
	    errMsg += Info->GenStr(iStep+1);
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

	  tagFound = FALSE;	  
	  for (Integer l=0; l<tagStrings.GetSize(); l++) {
	    SplitStringList(tagStrings[l], tagsAux);
	    for (Integer iTag=0; iTag<tagsAux.GetSize(); iTag++) 
	      // 	    if (tagsAux[iTag] == "anyTag" ||
	      // 		tagsAux[iTag] == tagsPerStep_[iStep][iPDE]) {
	      if (tagsAux[iTag] == tagsPerStep_[iStep][iPDE] ||
		  tagsAux[iTag] == "anyTag") {
		tagFound = TRUE;
		break;
	      }
	  }
	  
	  if (! tagFound) {
	    errMsg  = "MultiSequenceDriver::Init(): The section for ";
	    errMsg += "analysistype '";
	    errMsg += analysis;
	    errMsg += "' with the tag '";
	    errMsg += tagsPerStep_[iStep][iPDE];
	    errMsg += "' was not found in step ";
	    errMsg += Info->GenStr(iStep+1);
	    Error(errMsg.c_str(), __FILE__, __LINE__);
	  }
	  
	}
      } // iPDE
    } // iStep
  } 
  
} // end of namespace
