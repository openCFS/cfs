// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "SimState.hh"

#include <map>

// include hdf5 cpp file
#include "cpp/H5Cpp.h"

#include "General/Environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/SimInputHDF5.hh"
#include "DataInOut/ParamHandling/Xerces.hh"
#include "DataInOut/ParamHandling/XMLMaterialHandler.hh"
#include "PDE/SinglePDE.hh"
#include "Domain/Domain.hh"
#include "DataInOut/DefineInOutFiles.hh"
#include "DataInOut/ResultHandler.hh"


namespace CoupledField {

class MaterialHandler;

  SimState::SimState( bool useAsInput ) {
    ownWriter_ = false;
    isInitialized_ = false;
    hasInput_ = useAsInput;
    sequenceStep_ = 0;
    domain_ = NULL;
  }
  
  SimState::~SimState() {
    
  }
  
  void SimState::SetInputHdf5Reader(shared_ptr<SimInputHDF5> reader )  {
    inFile_ = reader;
  }

  Domain * SimState::GetDomain(UInt sequenceStep) {

    // Ensure, that hdf5 input class is present
    assert(inFile_);

    inFile_->InitModule();
    inFile_->DB_Init();

    // Get content of param and material file
    std::string paramContent, matContent;
    inFile_->DB_GetParamFileContent( paramContent );
    std::cerr << "parameter content is " << paramContent << std::endl;
    inFile_->DB_GetMatFileContent( matContent );
    std::cerr << "Material content is " << matContent<< std::endl;

    // Generate Xerces parameter reader
    std::string schema = progOpts->GetSchemaPathStr();
    schema += "/CFS-Simulation/CFS.xsd";

    Xerces * reader = new Xerces(schema);
    reader->SetString( paramContent );
    PtrParamNode rootNode = reader->CreateParamNodeInstance();
    delete reader;

    // Generate material reader
    DefineInOutFiles fileHandler;
    XMLMaterialHandler * matHandler = new XMLMaterialHandler();
    matHandler->LoadFromString( matContent );

    // Create dummy info node
    PtrParamNode infoNode(new ParamNode(ParamNode::INSERT, ParamNode::ELEMENT ));

    // Load grids
    std::map<std::string, shared_ptr<SimInput> > inFiles;
    std::map<std::string, StdVector<shared_ptr<SimInput> > > gridInputs;
    fileHandler.CreateSimInputFiles( rootNode, infoNode, inFiles, gridInputs );

    std::map<std::string, StdVector<shared_ptr<SimInput> > >::const_iterator it;
    it = gridInputs.begin();

    ResultHandler * resHandler = NULL;
    // Create new domain
    domain_ = new Domain( gridInputs, resHandler, matHandler, shared_from_this(),
                          rootNode, infoNode );
    std::cerr << "NEW DOMAIN : " << domain_ << std::endl;
    domain_->CreateGrid();

    // create driver for domain
    domain_->PostInit();
    return domain_;
  }
  
  void SimState::SetInputReaderToSameInput() {
    // assure, that an output reader is present
    assert( outFile_ );
    
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

    
  }
  
  void SimState::UpdateToStep( UInt stepNum ) {

    std::set<shared_ptr<BaseFeFunction> >::iterator it = feFcts_.begin();
    for( ; it != feFcts_.end(); ++it ) {
      // get pdeName
      std::string pdeName = (*it)->GetPDE()->GetName();
      std::string feName = 
          SolutionTypeEnum.ToString((*it)->GetResultInfo()->resultType);
      inFile_->DB_GetFeFctCoefs(1, stepNum, pdeName, feName,
                                (*it)->GetSingleVector() );
    }

  }
  
  void SimState::SetOutputHdf5Writer( shared_ptr<SimOutputHDF5> writer ) {
    if( writer ) {
      outFile_ = writer;
      ownWriter_ = false;
    } else {
      // create new writer with the name of the simulation
      std::string name = progOpts->GetSimName();
      
      // check for restart
      bool restart = progOpts->GetRestart();
      
      PtrParamNode infoNode(new ParamNode(ParamNode::INSERT, ParamNode::ELEMENT ));
      PtrParamNode h5Node (new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      PtrParamNode eFiles (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      eFiles->SetName("externalFiles");
      eFiles->SetValue( "false" );
      h5Node->AddChildNode(eFiles);
      outFile_.reset(new SimOutputHDF5( name, h5Node, infoNode, restart ));
      
      ownWriter_ = true;
    }
    
  }
  
  
  UInt SimState::GetLastStepNum() {
    // assert
   UInt lastStepNum = 0;
   
   std::map<std::string, std::set<std::string> > coefFcts;
   std::map<std::string, std::set<std::string> >::const_iterator pdeIt;
   inFile_->DB_GetAvailPdeCoefFcts( sequenceStep_, coefFcts );
   
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
       
       // get maximum
       lastStepNum = std::max(lastStepNum, stepValues.rbegin()->first);
     }
     
   } // loop: pdes
   return lastStepNum;
    
  }

  void SimState::SetMatParamFile( const fs::path& paramFile, 
                                  const fs::path& matFile ) {
    paramFile_ = paramFile;
    matFile_ = matFile;
  }


  void SimState::BeginMultiSequenceStep( UInt step, BasePDE::AnalysisType type ) {
    // Ensure initialized object
    Init();
    
    sequenceStep_ = step;
    outFile_->DB_BeginMultiSequenceStep(step, type);

  }

  void SimState::WriteStep( UInt stepNum, Double stepVal ) {

    // Write all coefficient function of the current step to the HDF5 file
    outFile_->DB_BeginStep(stepNum, stepVal);
    std::set<shared_ptr<BaseFeFunction> >::iterator it = feFcts_.begin();
    for( ; it != feFcts_.end(); ++it ) {
      std::string pdeName = (*it)->GetPDE()->GetName();
      std::string feName = 
          SolutionTypeEnum.ToString((*it)->GetResultInfo()->resultType);
      outFile_->DB_WriteFeFunction( pdeName, feName, (*it)->GetSingleVector() );
    }
  }

  void SimState::FinishMultiSequenceStep( ) {
    outFile_->FinishMultiSequenceStep();
    
    // delete all registered feFunctions
    feFcts_.clear();
    
  }

  void SimState::RegisterFeFct( shared_ptr<BaseFeFunction> feFct) {
    feFcts_.insert( feFct );
  }

  void SimState::Init() {
    // Leave if already initialized
    if( isInitialized_) 
      return;

    // Initialize external hdf5 file
    outFile_->DB_Init();

    outFile_->DB_WriteXmlFiles( paramFile_, matFile_ );
    
    isInitialized_ = true;
  }

} // namespace
