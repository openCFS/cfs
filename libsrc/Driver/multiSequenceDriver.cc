#include "multiSequenceDriver.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "Utils/tools.hh"

namespace CoupledField {

  // ***************
  //   Constructor
  // ***************
  MultiSequenceDriver::MultiSequenceDriver(Domain * adomain) 
    : BaseDriver(adomain) {
    ENTER_FCN( "MultiSequenceDriver::MultiSequenceDriver" );

    std::cerr << " ---- Creating 'MultiSequenceDriver' ---- " << std::endl;

    numSteps_ = 0;
    
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
    std::cerr << std::endl;
    std::cerr << "steps = " << std::endl << steps << std::endl;

    // 2.) ensure that all steps occur and none is left out
    Integer currStep = 0;
    for (Integer iStep=1; iStep<=steps.GetSize(); iStep++) {
      for (Integer j=0; j<steps.GetSize(); j++) {
	if (steps[j] == iStep) {
	  //std::cerr << "iStep = steps[" << j << "] = " << iStep << std::endl;
	  currStep++;
	  break;
	}
      }
    }
    std::cerr << "currStep = " << currStep << std::endl;
    
    std::string errMsg;
    if (currStep != steps.GetSize()) {
      errMsg = "MultiStepSequence: The indices for a multiSequence must start ";
      errMsg += "with 1 and have to be in consecutive order.\n";
      errMsg += "Please check your paramter file!";
      Error(errMsg.c_str(), __FILE__, __LINE__);
    }
    numSteps_ = steps.GetSize();
    
    // 3.) Resize 'outer' vectors
    std::cerr << "numSteps_ = " << numSteps_ << std::endl;
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
      Boolean isTrue = (stepString == "1");
      if (isTrue)
	std::cerr << "TRUE" << std::endl;
      std::cerr << "stepString = '" << stepString << "'" << std::endl;
      keyVec  = "multiSequence", "step", "pde", "tag";      
      params->GetList( keyVec, attrVec, valVec, tagsAux );

      std::cerr << "The tags for step " << stepString << " are: " << std::endl;
      std::cerr << tagsAux << std::endl << std::endl;
      

      keyVec  = "multiSequence", "step", "pde", "type";      
      params->GetList( keyVec, attrVec, valVec, pdesAux );    

      keyVec  = "multiSequence", "step", "pde", "analysis";      
      params->GetList( keyVec, attrVec, valVec, analysisAux );

      std::cerr << "The pdes for step " << stepString << " are: " << std::endl;
      std::cerr << pdesAux << std::endl << std::endl;

      std::cerr << "The analysis for step " << stepString << " are: " << std::endl;
      std::cerr << analysisAux << std::endl << std::endl;

      // Add pdes, tags, and analysistypes of current step
      for (Integer iPDE=0; iPDE<pdesAux.GetSize(); iPDE++) {
	
	String2Enum(analysisAux[iPDE], analysisType);
	
	pdesPerStep_[iStep].Push_back(pdesAux[iPDE]);
	tagsPerStep_[iStep].Push_back(tagsAux[iPDE]);
	analysisPerStep_[iStep].Push_back(analysisType);
	
      }
    }

    
    
    // 5.) Perform final consistency checks

    // iterate over all steps
    Boolean pdeFound = FALSE;
    Boolean tagFound = FALSE;
    std::string tagString;
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

	pdeFound = FALSE;
	for (Integer kPDE=0; kPDE<ptdomain_->GetNumPDE(); kPDE++)
	  if (pdesPerStep_[iStep][iPDE] == ptdomain_->GetPDE(kPDE)->GetName())
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
	std::cerr << "pdesPerStep_[iStep][iPDE] = " << pdesPerStep_[iStep][iPDE] << std::endl;
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
	// one time per step
	// ************************************************
	for (Integer kPDE=iPDE+1; kPDE<pdesPerStep_[iStep].GetSize(); kPDE++)
	  if (pdesPerStep_[iStep][iPDE] == pdesPerStep_[iStep][kPDE]) {
	     errMsg  = "MultiSequenceDriver::Init(): The PDE '";
	     errMsg += pdesPerStep_[iStep][kPDE];
	     errMsg += "' occured more than one time in step ";
	     errMsg += Info->GenStr(iStep+1);
	     Error(errMsg.c_str(), __FILE__, __LINE__); 
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
	  params->Get(keyVec, attrVec, valVec, tagString);
	  
	  SplitStringList(tagString, tagsAux);

	  tagFound = FALSE;
	  for (Integer iTag=0; iTag<tagsAux.GetSize(); iTag++) 
	    if (tagsAux[iTag] == "anyTag" ||
		tagsAux[iTag] == tagsPerStep_[iStep][iPDE]) {
	      tagFound = TRUE;
	      break;
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
