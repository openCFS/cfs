// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "SimState.hh"

#include <map>
#include <boost/bind/bind.hpp>

// include hdf5 cpp file
#include "H5Cpp.h"

#include "Utils/mathParser/mathParser.hh"
#include "General/Environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/SimInputHDF5.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "DataInOut/ParamHandling/XMLMaterialHandler.hh"
#include "PDE/SinglePDE.hh"
#include "Domain/Domain.hh"
#include "DataInOut/DefineInOutFiles.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "Driver/SingleDriver.hh"

namespace CoupledField {

// declare class specific logging stream
DEFINE_LOG(simState, "simState")

class MaterialHandler;

  SimState::SimState( bool useAsInput, Domain* parentDomain ) {
    LOG_DBG(simState) << " Created new SimState (address: " << this << ")";
    LOG_DBG(simState) << " \tHasInput: " << (useAsInput ? "yes" : "no");
    ownWriter_ = false;
    isInitialized_ = false;
    hasInput_ = useAsInput;
    sequenceStep_ = 0;
    domain_ = NULL;
    parentDomain_ = parentDomain;
    parentParser_ = NULL;
    parentHandle_= 0;
    interpol_ = NO_INTERPOLATION;
    
  }
  
  SimState::~SimState() {
    LOG_DBG(simState) << " Deleted SimState (address: " << this << ")";
    LOG_DBG(simState) << " \tHasInput: " << (hasInput_ ? "yes" : "no");
    // close in file
    //inFile_.reset();
    outFile_.reset();
    feFcts_.clear();
    
    // de-register callback
    conn_.disconnect();
    if( parentParser_ ) {
      parentParser_->ReleaseHandle(parentHandle_);
      parentParser_ = NULL;
    }
    
  }
  
  void SimState::SetInputHdf5Reader(shared_ptr<SimInputHDF5> reader )  {
    inFile_ = reader;
    LOG_DBG(simState) << " Set input reader with file '" 
                      << reader->GetFileName() << "'";
        
  }

  Domain* SimState::GetDomain(UInt sequenceStep, const GridMap& map ) {

    LOG_TRACE(simState) << " GetDomain for sequenceStep " << sequenceStep;
    if( map.size() > 0 ) {
      LOG_TRACE(simState) << "\t=> Grid map was provided";
    }
    // Ensure, that hdf5 input class is present
    assert(inFile_);

    // Initialize in/out readers
    inFile_->InitModule();
    inFile_->DB_Init();
    
    // Check, that sequence step is set correctly 
    UInt lastMsStep = GetLastMsStepNum();
    if( sequenceStep > lastMsStep) {
      EXCEPTION( "Can not return domain for sequence step " << sequenceStep
                 << " as only " << lastMsStep << " multisequence steps are "
                 << "defined" );
    }

    // Attention: currently we may not re-set the sequence step, if the simState
    // is primary intended for output
    sequenceStep_ = sequenceStep;

    // Get content of param and material file
    std::string paramContent, matContent;
    inFile_->DB_GetParamFileContent( paramContent );
    inFile_->DB_GetMatFileContent( matContent );
    LOG_DBG3(simState) << "Content of Parameter file:\n" << paramContent;
    LOG_DBG3(simState) << "Content of Material file:\n" << matContent;

    // Generate xml parameter reader
    LOG_TRACE(simState) << "Generating parameter node from xml file";
    std::string schema = progOpts->GetSchemaPathStr() + "/CFS-Simulation/CFS.xsd";
    PtrParamNode rootNode = XmlReader::ParseString(paramContent, schema,
                                                   "http://www.cfs++.org/simulation");


    // Generate material reader
    
    shared_ptr<XMLMaterialHandler> matHandler = boost::make_shared<XMLMaterialHandler>();
    matHandler->LoadFromString( matContent );

    // Create dummy info node
    PtrParamNode infoNode = ParamNode::GenerateWriteNode("", ""); // empty filename means we don't write and ignore ParamNode::ToFile()

    // Load grids only if not provided 
    std::map<std::string, shared_ptr<SimInput> > inFiles;
    std::map<std::string, StdVector<shared_ptr<SimInput> > > gridInputs;
    
    // If no external grid map is provided, we have to generate the
    // input readers
    if( map.size() == 0 ) {
      DefineInOutFiles fileHandler;
      fileHandler.CreateSimInputFiles( rootNode, infoNode, inFiles, gridInputs );
    }

    std::map<std::string, StdVector<shared_ptr<SimInput> > >::const_iterator it;
    it = gridInputs.begin();

    ResultHandler * resHandler = parentDomain_->GetResultHandler();
    // Create new domain
    LOG_TRACE(simState) << "Generating domain";
    domain_ = new Domain( gridInputs, resHandler, matHandler, shared_from_this(),
                          rootNode, infoNode, false );
    LOG_DBG(simState) << "\t=> Generated new Domain (address: " << domain_ << ")";
    // Provide external grid map if non-empty
    if( map.size() > 0 ) {
      domain_->SetGridMap(map);
    }
    
    LOG_TRACE(simState) << "Generating grid";
    domain_->CreateGrid();

    // create driver for domain and set it to given sequenceStep
    LOG_TRACE(simState) << "Initializing drivers and PDEs";
    domain_->PostInit(sequenceStep);
    
    LOG_TRACE(simState) << "Finished GetDomain for sequenceStep " << sequenceStep;
    return domain_;
  }
  
  bool SimState::SetInputReaderToSameOutput() {
    // assure, that an output reader is present
    assert( outFile_ );

    // if input reader is already set to outpu writer, just leave
    if(inFile_)
      return true;
    try {
      // get file name of output writer
      fs::path h5FileName = outFile_->GetFileName();

      // Create dummy info node
      PtrParamNode infoNode(new ParamNode(ParamNode::INSERT, ParamNode::ELEMENT ));

      // generate input reader for same file
      PtrParamNode node(new ParamNode());
      boost::shared_ptr<SimInputHDF5> in;
      inFile_.reset(new SimInputHDF5(h5FileName.string(), node, infoNode));

      // Initialize module
      inFile_->InitModule();
      inFile_->DB_Init();
      return true;
    } catch( Exception& e ) {
      return false;
    }
    return false;

    
  }
  
  void SimState::UpdateToStep( UInt sequenceStep, UInt stepNum ) {
    LOG_TRACE(simState) << " Updating input to sequence step "
        << sequenceStep << " with step number " << stepNum;
    
    std::set<shared_ptr<BaseFeFunction> >::iterator it = feFcts_.begin();
    for( ; it != feFcts_.end(); ++it ) {
      // get pdeName
      std::string pdeName = (*it)->GetPDE()->GetName();
      std::string feName = 
          SolutionTypeEnum.ToString((*it)->GetResultInfo()->resultType);
      inFile_->DB_GetFeFctCoefs(sequenceStep, stepNum, pdeName, feName,
                                (*it)->GetSingleVector() );
    }
    
    // In case there is a domain object, we also have to notify the 
    // Driver
    if( domain_ ) {
      // notify also the driver about the current time / frequency step
      SingleDriver * ptDriver = domain_->GetSingleDriver();
      Integer index = stepNums_.Find(stepNum);
      if( index == -1 ) {
        EXCEPTION( "Step number " << stepNum << " not defined" );
      }
      Double stepVal = stepVals_[index];
      ptDriver->SetToStepValue( stepNum, stepVal);
    }
  }
  
  void SimState:: SetInterpolation( InterpolType type, MathParser * parser,
                                    BasePDE::AnalysisType analysis,
                                    Double offset ) {
    
    // ensure, that domain was already obtained
    if( ! domain_ ) {
      EXCEPTION( "Please call GetDomin() before setting interpolation");
    }
    
    // store type of interpolation
    interpol_ = type;
    
    
    // store parser and register callback function for step / freq update
    parentParser_ = parser;
    parentHandle_ = parentParser_->GetNewHandle();
    // register call-back function for update
    if(analysis == BasePDE::TRANSIENT ||
        analysis == BasePDE::STATIC ) { 
      parentParser_->SetExpr(parentHandle_, "t");
    } else {
      parentParser_->SetExpr(parentHandle_, "f");
    }
    
    // only add callback, if interpolation is not CONSTANT
    if( type != CONSTANT ) {
      conn_ = parentParser_->AddExpChangeCallBack(
          boost::bind(&SimState::UpdateTimeFreqStep, this ),
          parentHandle_ );
    }
    
    
    // Loop over all registered FeFunctions
    std::map<std::string, std::set<std::string> > coefFcts;
    std::map<std::string, std::set<std::string> >::const_iterator pdeIt;
    inFile_->DB_GetAvailPdeCoefFcts( sequenceStep_, coefFcts );

    std::map<UInt,Double> stepNumVals;
    // Loop over all available PDEs
    for(pdeIt = coefFcts.begin(); pdeIt != coefFcts.end(); ++pdeIt ) {

      const std::string & pdeName = pdeIt->first;
      const std::set<std::string> & coefs = pdeIt->second;
      std::set<std::string>::const_iterator coefIt = coefs.begin();

      // Loop over all available CoefFcts
      for( ; coefIt != coefs.end(); ++coefIt ) {
        const std::string & coefName = *coefIt; 
        std::map<UInt, Double> stepValues;
        inFile_->DB_GetStepValues( sequenceStep_, pdeName, coefName, stepValues );

        // insert step values in to own structure
        std::map<UInt, Double>::const_iterator stepIt = stepValues.begin();
        for( ; stepIt != stepValues.end(); ++stepIt ) {
          stepNumVals[stepIt->first] = stepIt->second;
        }
      }
    } // loop: pdes
    
    // Now we have all step numbers / vals. Copy them to stepNums_ and 
    // stepVals_
    stepNums_.Resize(stepNumVals.size());
    stepVals_.Resize(stepNumVals.size());
    std::map<UInt, Double>::const_iterator stepIt = stepNumVals.begin();
    UInt pos = 0;
    for( ; stepIt != stepNumVals.end(); ++stepIt, ++pos ) {
      stepNums_[pos] = stepIt->first;
      stepVals_[pos] = stepIt->second;
      LOG_DBG(simState) << "stepNums_["<<pos<<"] = " <<  stepNums_[pos] << "\n";
      LOG_DBG(simState) << "stepVals_["<<pos<<"] = " <<  stepVals_[pos] << "\n";
    }
    
  }
  
  void SimState::SetOutputHdf5Writer( shared_ptr<SimOutputHDF5> writer ) {
    LOG_TRACE(simState) << " Setting HDF5 writer";
    
    if( writer ) {
    
      LOG_TRACE(simState) << "\t=> External writer with name '" 
          << writer->GetFileName() << "' was provided.";
      outFile_ = writer;
      ownWriter_ = false;
    
    } else {
      // create new writer with the name of the simulation
      std::string name = progOpts->GetSimName();
      LOG_TRACE(simState) << "\t=> Creating own writer with base name '" 
                << name << "'.";
      
      // check for restart
      bool restart = progOpts->GetRestart();
      
      PtrParamNode infoNode = ParamNode::GenerateWriteNode("", ""); // empty filename means we don't write and ignore ParamNode::ToFile()
      PtrParamNode h5Node (new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      PtrParamNode eFiles (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      eFiles->SetName("externalFiles");
      eFiles->SetValue( "false" );
      h5Node->AddChildNode(eFiles);
      outFile_.reset(new SimOutputHDF5( name, h5Node, infoNode, restart ));
      
      ownWriter_ = true;
    }
    
  }
  
  void SimState::GetSequenceSteps( std::map<UInt, BasePDE::AnalysisType>& analysis,
                                   std::map<UInt, Double>& accTime,
                                   std::map<UInt, bool>& isFinished ) {
    inFile_->DB_GetNumMultiSequenceSteps( analysis, accTime, isFinished );
  }
  
  bool SimState::IsCompleted( UInt sequenceStep ) {
    std::map<UInt, BasePDE::AnalysisType>  analysis;
    std::map<UInt, Double> accTime;
    std::map<UInt, bool>  isFinished;
    inFile_->DB_GetNumMultiSequenceSteps( analysis, accTime, isFinished );

    if( isFinished.find(sequenceStep) == isFinished.end() ) {
      EXCEPTION( "Sequence step " << sequenceStep << " is not contained in "
          << "simulation file '" << inFile_->GetFileName() << "'.");
    }
    return isFinished[sequenceStep];
  }
  
  UInt SimState::GetLastMsStepNum() {
    
    UInt lastMsStep = 0;
    std::map<UInt, BasePDE::AnalysisType>  analysis;
    std::map<UInt, Double> accTime;
    std::map<UInt, bool>  isFinished;
    inFile_->DB_GetNumMultiSequenceSteps( analysis, accTime, isFinished );
    std::map<UInt, BasePDE::AnalysisType> ::const_iterator msIt = analysis.begin();
    for( ; msIt != analysis.end(); ++msIt) {
      lastMsStep = std::max(lastMsStep, msIt->first);
    }
    return lastMsStep;
  }
  
  void SimState::GetLastStepNum(UInt sequenceStep, UInt& lastStepNum,
                                Double& lastStepVal ) {
   lastStepNum = 0;
   lastStepVal = 0.0;
   std::map<std::string, std::set<std::string> > coefFcts;
   std::map<std::string, std::set<std::string> >::const_iterator pdeIt;
   inFile_->DB_GetAvailPdeCoefFcts( sequenceStep, coefFcts );
   
   // Loop over all available PDEs
   for(pdeIt = coefFcts.begin(); pdeIt != coefFcts.end(); ++pdeIt ) {
   
     const std::string & pdeName = pdeIt->first;
     const std::set<std::string> & coefs = pdeIt->second;
     std::set<std::string>::const_iterator coefIt = coefs.begin();
     
     // Loop over all available CoefFcts
     for( ; coefIt != coefs.end(); ++coefIt ) {
       const std::string & coefName = *coefIt; 
       std::map<UInt, Double> stepValues;
       inFile_->DB_GetStepValues( sequenceStep, pdeName, coefName, stepValues );
       
       // get maximum
       lastStepNum = std::max(lastStepNum, stepValues.rbegin()->first);
       lastStepVal = std::max(lastStepVal, stepValues.rbegin()->second);
     }
   } // loop: pdes
    
  }

  void SimState::SetMatParamFile( const fs::path& paramFile, 
                                  const fs::path& matFile ) {
    paramFile_ = paramFile;
    matFile_ = matFile;
  }


  void SimState::BeginMultiSequenceStep( UInt step, BasePDE::AnalysisType type ) {
    LOG_TRACE(simState) << "Begin new MS step " << step  << " of type "   << BasePDE::analysisType.ToString( type );
    // Ensure initialized object
    Init();
    
    sequenceStep_ = step;
    outFile_->DB_BeginMultiSequenceStep(step, type);

  }

  void SimState::WriteStep( UInt stepNum, Double stepVal ) {
    LOG_TRACE(simState) << "Writing step " << stepNum << " with value "
        << stepVal << " s / Hz";

    // Write all coefficient function of the current step to the HDF5 file
    outFile_->DB_BeginStep(stepNum, stepVal);
    std::set<shared_ptr<BaseFeFunction> >::iterator it = feFcts_.begin();
    for( ; it != feFcts_.end(); ++it ) {
      std::string pdeName = (*it)->GetPDE()->GetName();
      std::string feName = 
          SolutionTypeEnum.ToString((*it)->GetResultInfo()->resultType);
      LOG_DBG(simState) << "\tWriting feFunction " << feName 
                        << " of PDE " << pdeName;
      outFile_->DB_WriteFeFunction( pdeName, feName, (*it)->GetSingleVector() );
    }
  }

  void SimState::FinishMultiSequenceStep( bool completed, Double accTime) {
    outFile_->DB_FinishMultiSequenceStep( completed, accTime );
    
    // delete all registered feFunctions
    feFcts_.clear();
    
  }

  void SimState::RegisterFeFct( shared_ptr<BaseFeFunction> feFct) {
    feFcts_.insert( feFct );
  }

  void SimState::Finalize() {
      feFcts_.clear();
      inFile_.reset();
    }
  
  void SimState::Init() {
    // Leave if already initialized
    if( isInitialized_) 
      return;

    LOG_TRACE(simState) << "Initializing simState for the first time";
    // Initialize external hdf5 file
    outFile_->DB_Init();

    outFile_->DB_WriteXmlFiles( paramFile_, matFile_ );
    
    isInitialized_ = true;
    LOG_TRACE(simState) << "Finished initialization";
  }

  void SimState::UpdateTimeFreqStep() {
    Double actStepVal = parentParser_->Eval( parentHandle_ ); 

    LOG_TRACE(simState)  << "Updating Time / Freq value to " 
        << actStepVal << " s / Hz.";

    // get current time / freq step
    SingleDriver * ptDriver = domain_->GetSingleDriver();

    // ======================================
    // determine "step" in the child driver
    // ======================================
    Double w1, w2;
    Integer index1, index2;

    Double yValue = 0.0;

    // get index of last element
    const UInt kend = stepVals_.GetSize() - 1;
    LOG_TRACE(simState)  << "stepVals_ = " << stepVals_.ToString() << " kend = " << kend << "\n";

    // if coordinate is out of bounds or we have just one entry,
    // return boundary value (i.e.first or last) 
    if ( actStepVal > stepVals_[kend] || kend == 0) {
      yValue = stepNums_[kend];
      index1 = index2 = kend;
      w1 = 1.0;
      w2 = 0.0;
    }
    else if ( actStepVal < stepVals_[0] ) {
      yValue = stepNums_[0];
      index1 = index2 = 0;
      w1 = 1.0;
      w2 = 0.0;
    }
    else {

      UInt k;
      index1=0;
      index2=kend;
      // We will find the right place in the table by means of bisection.
      //  index1 and index2 bracket the input value of xEntry
      while (index2-index1 > 1) {
        k=(index2+index1) >> 1; // binary right shift
        if (stepVals_[k] > actStepVal)
          index2=k;
        else
          index1=k;
      }

      // size of x interval
      Double dxVal=stepVals_[index2] - stepVals_[index1];

      // The x-values must be distinct!
      if (dxVal == 0.0) {
        EXCEPTION("You cannot have two equal x values!" );
      }

      // relative distance of xEntry to x-Value bounds
      w1= ( stepVals_[index2] - actStepVal )/dxVal;
      w2 = ( actStepVal - stepVals_[index1] )/dxVal;


      // Perform interpolation    
      switch ( interpol_ ) {

        case NEAREST_NEIGHBOR:
          if ( w1 <= 0.5 ) {
            yValue = stepNums_[index2];
            w1 = 1;
            w2 = 0;
            index1 = index2;
            index2 = 0;
          } else {
            yValue = stepNums_[index1];
            w1 = 1;
            w2 = 0;
            index2 = 0;
          }
          break;

        case LINEAR:
          yValue = w1 * stepNums_[index1] + w2 * stepNums_[index2];
          break;
        default:
          EXCEPTION( "Interpolation type not known");
          break;
      }
    }

    LOG_DBG(simState) << "\trequested stepNum is " << yValue << " (interpolated)";
    LOG_DBG(simState) << "\tindex1: " << index1+1 << " (weight: " << w1 << ")";
    LOG_DBG(simState) << "\tindex2: " << index2+1 << " (weight: " << w2 << ")";

    // ===================================================
    //  Update all coefficients with correspdoning weight
    // ===================================================
    std::set<shared_ptr<BaseFeFunction> >::iterator it = feFcts_.begin();
    for( ; it != feFcts_.end(); ++it ) {
      std::string pdeName = (*it)->GetPDE()->GetName();
      std::string feName = 
          SolutionTypeEnum.ToString((*it)->GetResultInfo()->resultType);
      LOG_DBG(simState) << "\tUpdating feFunction " << feName 
          << " of PDE " << pdeName;

      if( interpol_ == NEAREST_NEIGHBOR ) {
        LOG_DBG(simState) << "\t=> Nearest neighbor update ";
        // In this case we can directly access the values
        inFile_->DB_GetFeFctCoefs(sequenceStep_, index1+1, pdeName, feName,
                                  (*it)->GetSingleVector() );

      } else if ( interpol_ == LINEAR ) {
        LOG_DBG(simState) << "\t=> Linear interpolation update ";

        // get both involved step indices and calculate the weighted sum
        if( ptDriver->IsComplex()) {
          Vector<Complex> & coefVec = 
              dynamic_cast<Vector<Complex>& >(*(*it)->GetSingleVector());
          Vector<Complex> vals1, vals2;
          vals1.Resize(coefVec.GetSize());
          vals2.Resize(coefVec.GetSize());
          // get values for both steps 
          inFile_->DB_GetFeFctCoefs(sequenceStep_, index1+1, pdeName, feName, &vals1 );
          inFile_->DB_GetFeFctCoefs(sequenceStep_, index2+1, pdeName, feName, &vals2 );

          // assemble weighted sum
          coefVec = vals1 * w1;
          coefVec += (vals2 * w2);
        } else {

          Vector<Double> & coefVec = 
              dynamic_cast<Vector<Double>& >(*(*it)->GetSingleVector());
          Vector<Double> vals1, vals2;
          vals1.Resize(coefVec.GetSize());
          vals2.Resize(coefVec.GetSize());
          // get values for both steps 
          inFile_->DB_GetFeFctCoefs(sequenceStep_, index1+1, pdeName, feName, &vals1 );
          inFile_->DB_GetFeFctCoefs(sequenceStep_, index2+1, pdeName, feName, &vals2 );

          // assemble weighted sum
          coefVec = vals1 * w1;
          coefVec += (vals2 * w2);
        }
      }
    }

    // update driver to time / freq step (initially 
    // we can also just set the value of the related mathParser
    // to the specified value and circumvent the driver itself)
    ptDriver->SetToStepValue(index1+1, actStepVal);

  }

  // ************************************************************************
  // ENUM INITIALIZATION
  // ************************************************************************

  // Definition of finite element space types
   static EnumTuple interpolTypeTuples[] = {
     EnumTuple(SimState::NO_INTERPOLATION, "undef"),
     EnumTuple(SimState::CONSTANT, "constant"), 
     EnumTuple(SimState::NEAREST_NEIGHBOR, "nearestNeighbor"), 
     EnumTuple(SimState::LINEAR,           "linear")
   };
   
   Enum<SimState::InterpolType> SimState::InterpolTypeEnum = \
      Enum<SimState::InterpolType>("Types of interpolation",
          sizeof(interpolTypeTuples) / sizeof(EnumTuple),
          interpolTypeTuples);

} // namespace
