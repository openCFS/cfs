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

    // current driver
    SingleDriver * actDriver;

    // Vector of solutions
    StdVector<Vector<Double> > sols;
    Integer level = 0;
    Integer nextStep;
  
    // Print out the grid
    ptdomain_->PrintGrid(level);

    Info->StartProgress("Starting to solve problem",FALSE);
    // outer loop over all single sequences
    for (Integer iStep=0; iStep<numSteps_; iStep++) {
      
      
      Info->WriteMultiSequenceStep(iStep+1, 
				   analysisPerStep_[iStep][0]);
      // The first time the PDEs are already
      // initialized
      if (iStep > 0)
	ptdomain_->ResetPDEs();

      // Since per time step only one type of 
      // analysis is allowed, we simple access
      // the first entry fo analysisPerStep_
      if (analysisPerStep_[iStep][0] == STATIC) {
	actDriver = new StaticDriver(ptdomain_, actStep_, actTime_,
				     tagsPerStep_[iStep][0], TRUE);
	nextStep = actStep_ + 1;;
      }
      else if (analysisPerStep_[iStep][0] == TRANSIENT) {      
	actDriver = new TransientDriver(ptdomain_, actStep_, actTime_, 
					tagsPerStep_[iStep][0],TRUE);
	// here the time step and time 
	// have to be imcremented
	
      }
      else if (analysisPerStep_[iStep][0] == HARMONIC) {
	actDriver = new HarmonicDriver(ptdomain_, actStep_, actTime_,
				       tagsPerStep_[iStep][0], TRUE);
      }
      
      // Initialize all PDEs
      ptdomain_->InitPDEs(actStep_+1, tagsPerStep_[iStep]);

      // Give the according pdes to the driver
      actDriver->SetPDEs(pdesPerStep_[iStep]);

      // After the first run, initialize this PDE
      // with the solution of the previous run
      if (iStep > 0) {
	for (Integer i=0; i<pdesPerStep_[iStep].GetSize(); i++)
	  pdesPerStep_[iStep][i]->SetSolution(sols[i]);
      }
      
      // Solve Problem
      actDriver->SolveProblem();
      
      // Get solution for next step and delete
      // all PDEs
      Integer iPDE, kPDE;
      if (iStep < numSteps_-1) {
	sols.Resize(pdesPerStep_[iStep+1].GetSize());
	
	// Iterate over all PDEs in the next step
	for (iPDE=0; iPDE<pdesPerStep_[iStep+1].GetSize(); iPDE++)
	  // Iterate over all PDEs in the current step
	  for(kPDE=0; kPDE<pdesPerStep_[iStep].GetSize(); kPDE++)
	    // If both match, then save the result of this step
	    // for the next step
	    if (pdesPerStep_[iStep+1][iPDE]->GetName() ==
		pdesPerStep_[iStep][kPDE]->GetName())
	      //dynamic_cast<const NodeStoreSol<Double>& >
	      pdesPerStep_[iStep][kPDE]->GetSolution().GetAlgSysVector(sols[iPDE]);
      }
      
      // delete analysistypes
      delete actDriver;
      
      // increase stepNumber
      actStep_ = nextStep;
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
//     std::cerr << std::endl;
//     std::cerr << "steps = " << std::endl << steps << std::endl;

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
//     std::cerr << "currStep = " << currStep << std::endl;
    
    std::string errMsg;
    if (currStep != steps.GetSize()) {
      errMsg = "MultiStepSequence: The indices for a multiSequence must start ";
      errMsg += "with 1 and have to be in consecutive order.\n";
      errMsg += "Please check your paramter file!";
      Error(errMsg.c_str(), __FILE__, __LINE__);
    }
    numSteps_ = steps.GetSize();
    
    // 3.) Resize 'outer' vectors
//     std::cerr << "numSteps_ = " << numSteps_ << std::endl;
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
//       if (isTrue)
// 	std::cerr << "TRUE" << std::endl;
//       std::cerr << "stepString = '" << stepString << "'" << std::endl;
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
      for (Integer l=0; l<ptdomain_->GetNumPDE(); l++) 
	for (Integer iPDE=0; iPDE<pdesAux.GetSize(); iPDE++) {
	  if (ptdomain_->GetPDE(l)->GetName() == pdesAux[iPDE]) {
	    pdesPerStep_[iStep].Push_back(ptdomain_->GetPDE(l));
	    String2Enum(analysisAux[iPDE], analysisType);
	      
	    tagsPerStep_[iStep].Push_back(tagsAux[iPDE]);
	    analysisPerStep_[iStep].Push_back(analysisType);
// 	    std::cerr << "--> true " << std::endl;
	    break;
	  } 
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
	  if (pdesPerStep_[iStep][iPDE]->GetName() == 
	      ptdomain_->GetPDE(kPDE)->GetName())
	    pdeFound = TRUE;

	if (!pdeFound) {
	  errMsg = "MultiSequenceDriver::Init(): The PDE '";
	  errMsg += pdesPerStep_[iStep][iPDE]->GetName();
	  errMsg += "' is not mentioned in the 'pdeList'!";
	  Error(errMsg.c_str(), __FILE__, __LINE__);
	}

	// ************************************************
	// Check if the specified tags per pde and per step
	// occur in the 'bcsAndLoads'-section
	// ************************************************

	//get tags for current pde
	keyVec = "pdeList", pdesPerStep_[iStep][iPDE]->GetName(), "bcsAndLoads", "tag";
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
	  errMsg += pdesPerStep_[iStep][iPDE]->GetName();
	  errMsg += "'!";
	  Error(errMsg.c_str(), __FILE__, __LINE__);
	}
	  
	// ************************************************
	// Check if one PDE occured more than 
	// one time per step and if for each step
	// only one type of analysistype occured
	// ************************************************
	for (Integer kPDE=iPDE+1; kPDE<pdesPerStep_[iStep].GetSize(); kPDE++) {
	  if (pdesPerStep_[iStep][iPDE]->GetName() == 
	      pdesPerStep_[iStep][kPDE]->GetName()) {
	    errMsg  = "MultiSequenceDriver::Init(): The PDE '";
	    errMsg += pdesPerStep_[iStep][kPDE]->GetName();
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
