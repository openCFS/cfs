// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "Domain.hh"

#include <set>
#include <map>
#include <memory>
#include <boost/filesystem.hpp>

#include "def_use_openmp.hh"

#include "General/Environment.hh"
#include "General/Exception.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "Domain/CoordinateSystems/CartesianCoordSystem.hh"
#include "Domain/CoordinateSystems/CylCoordSystem.hh"
#include "Domain/CoordinateSystems/PolarCoordSystem.hh"
#include "Domain/CoordinateSystems/TrivialCartesianCoordSystem.hh"
#include "Domain/CoordinateSystems/DefaultCoordSystem.hh"
#include "Domain/Results/BaseResults.hh"

#ifdef USE_OPENMP
#include "Utils/mathParser/mathParserOMP.hh"
#endif
#include "Utils/mathParser/mathParser.hh"

#include "DataInOut/SimInput.hh"
#include "DataInOut/ParamHandling/MaterialHandler.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/SimInput.hh"
#include "General/Exception.hh"

#include "Optimization/Design/DensityFile.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Optimization.hh"


#include "DataInOut/ResultHandler.hh"

// Single Field PDEs
#include "PDE/AcousticPDE.hh"
#include "PDE/AcousticPDE_BEM.hh"
#include "PDE/AcousticMixedPDE.hh"
#include "PDE/AcousticSplitPDE.hh"
#include "PDE/ElecPDE.hh"
#include "PDE/PerturbedFlowPDE.hh"
#include "PDE/FlowPDE.hh"
#include "PDE/LinFlowPDE.hh"
#include "PDE/HeatPDE.hh"
#include "PDE/MagneticPDE.hh"
#include "PDE/MagneticScalarPotentialPDE.hh"
#include "PDE/MagEdgePDE.hh"
#include "PDE/MagEdgeMixedAVPDE.hh"
#include "PDE/MagEdgeSpecialAVPDE.hh"
#include "PDE/DarwinPDE.hh"
#include "PDE/MechPDE.hh"
#include "PDE/SmoothPDE.hh"
#include "PDE/TestPDE.hh"
#include "PDE/ElecCurrentPDE.hh"
#include "PDE/ElecQuasiStaticPDE.hh"
#include "PDE/WaterWavePDE.hh"
#include "PDE/LatticeBoltzmannPDE.hh"

// Coupling of Single PDEs
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "CoupledPDE/IterCoupledPDE.hh"
#include "CoupledPDE/PiezoCoupling.hh"
#include "CoupledPDE/MagnetoStrictCoupling.hh"
#include "CoupledPDE/AcouMechCoupling.hh"
#include "CoupledPDE/FluidMechCoupling.hh"
#include "CoupledPDE/WaterWaveAcousticsCoupling.hh"
#include "CoupledPDE/WaterWaveMechCoupling.hh"
#include "CoupledPDE/LinFlowHeatCoupling.hh"
#include "CoupledPDE/LinFlowAcouCoupling.hh"
#include "CoupledPDE/LinFlowMechCoupling.hh"
// Include driver
#include "Driver/BaseDriver.hh"
#include "Driver/SingleDriver.hh"
#include "Driver/MultiSequenceDriver.hh"

namespace CoupledField
{

// **************
//   Construtor
// **************

Domain::Domain(
    std::map<std::string, StdVector<shared_ptr<SimInput> > >& gridInputs,
    ResultHandler * handler, MaterialHandler * ptMat, 
    shared_ptr<SimState> simState,
    PtrParamNode rootNode, PtrParamNode infoNode,
    bool output)
{
  // initialize data
  numSinglePde_ = 0;
  numDirectCoupledPde_ = 0;
  numIterCoupledStdPde_ = 0;
  useExternalGridMap_ = false;
  isParentDomain_ = output;

  // Create new MathParser
#ifdef USE_OPENMP
  mathParser_ = new MathParserOMP();
#else
  mathParser_ = new MathParser();
#endif

  // assign pointers
  gridInputs_ = gridInputs;
  simState_ = simState;
  resultHandler_ = handler;
  ptIterCoupledPde_ = NULL;
  ptSingleDriver_ = NULL;
  multiSequenceDriver_ = NULL;
  param_ = rootNode;
  info_ = infoNode;
  
  ptMatHandler_ = ptMat;
  ptMatHandler_->SetDomain( this );
  
  // register variables defined in "variableList" element
  RegisterVariables();
}

void Domain::CreateGrid()
{
  boost::shared_ptr<Timer> timer = info_->Get(ParamNode::HEADER)->Get("domain/grids/timer")->AsTimer();
  timer->Start(); // called in simulation once, or with simstate agains

  std::string probGeo;
  Double depth2d = 1.0;
  
  param_->Get("domain")->GetValue("geometryType", probGeo);
  if (probGeo == "3d")
  {
    dim_ = 3;
  }
  else if (probGeo == "axi" || probGeo == "plane")
  {
    dim_ = 2;
    
    if( (probGeo == "plane") && (param_->Get("domain")->Has("depth2dPlane")) ){      
      depth2d = param_->Get("domain")->Get("depth2dPlane")->As<double>();
    }
  }
  else
  {
    EXCEPTION( "Wrong Dimension parameter in xml-File" );
  }

  
  // Check, if there is already a pre-defined set of grids
  if (gridMap_.size() > 0 && useExternalGridMap_) {
    // -------------------
    //  External Grid Map
    // -------------------
    // grids are already present, only coordinate systems are missing
    CreateCoordinateSystems();
    
  } else {

    // -------------------
    //  Own Grid Map
    // -------------------

    // make sure that there is a grid with ID "default", which is read first
    if ( gridInputs_.size() == 1 ) {
      if (gridInputs_.find("default") == gridInputs_.end())
      {
        WARN("Grid '" << gridInputs_.begin()->first
            << "' was renamed to 'default', as it is the only one.");
      }
      
      ReadGrid("default", gridInputs_.begin()->second, probGeo == "axi", depth2d);
      gridInputs_.clear();
    } else {
      if (gridInputs_.find("default") == gridInputs_.end()) {
        EXCEPTION("There is no grid with ID 'default'.");
      }
      
      ReadGrid("default", gridInputs_["default"], probGeo == "axi", depth2d);
    }

    // Create grids using input readers.
    std::map<std::string, StdVector<shared_ptr<SimInput> > >::const_iterator
        gridIt = gridInputs_.begin(),
        endIt = gridInputs_.end();

    // iterate over all other grid IDs
    for ( ; gridIt != endIt; ++gridIt )
    {
      if ( gridIt->first != "default") {
        ReadGrid( gridIt->first, gridIt->second, probGeo == "axi", depth2d );
      }
    } // loop: input readers
  } // if: use of external grids


  if (!progOpts->GetPrintGrid() == true)
  {
    // Initialize resultHandler
    if( resultHandler_ && isParentDomain_) {
      resultHandler_->Init( gridMap_, false);
    }

    // print grid information to result file if requested
    if(param_->Get("domain")->Get("printGridInfo")->As<bool>() ) {
      // only print grid information, if we are the master domain (i.e. not
      // within a child domain when reading from an external simulation) and if
      // we are not in a restarted state (i.e. we assume that the grid information was
      // printed already in the first attempt)
      if( resultHandler_ && 
          isParentDomain_&&
          !progOpts->GetRestart()) {
        gridMap_["default"]->CreateGridInformation(resultHandler_, coordSys_);
      }
    }
  }

  timer->Stop();
}

void Domain::ReadGrid(const std::string & gridId,
                      const StdVector< shared_ptr<SimInput> > & inputs,
                      bool isAxi, Double depth2d)
{
  // Type of mesh library. May either be cfsGrid or adaptGrid
  // (if AdaptGrid should be ever revived.)
  std::string libmesh = "cfsGrid";

  if( isParentDomain_) {
    std::cout << "++ Reading mesh ... " << std::flush; // endl after check for regularity
  }

  // iterate over all inputs for the current grid and init reader
  for (UInt iFile=0, numFiles=inputs.GetSize(); iFile < numFiles; ++iFile)
  {
    shared_ptr<SimInput> actInFile = inputs[iFile];
    actInFile->InitModule();
  }

  // create new grid
  Grid * actGrid = NULL;
  if (libmesh == "cfsGrid")
  {
    if (gridId == "default")
    {
      actGrid = new GridCFS(dim_, param_, info_, gridId);
    }
    else
    {
      actGrid = new GridCFS( inputs[0]->GetDim(), param_, info_, gridId);
    }
  }
  else
  {
    EXCEPTION( "Type of mesh_library should be one of "
        << "'cfsgrid' or 'adaptgrid', but is '" << libmesh << "'" );
  }

  // set flag about axisymmetry
  actGrid->SetAxi( isAxi );
  
  // set flag about depth in 2dplane
  if( (!isAxi) && (dim_ == 2) ){
    actGrid->SetDepth2dPlane(depth2d);
  }

  // add grid to internal map
  gridMap_[gridId] = actGrid;

  // Read in coordinate systems
  // This has to be done, before the grid gets finalized,
  // as some methods in the grid rely on the existence of coordinate
  // systems
  if (gridId == "default")
    CreateCoordinateSystems();

  // iterate over all inputs for the current grid and read mesh
  for (UInt iFile=0, numFiles=inputs.GetSize(); iFile < numFiles; ++iFile)
  {
    shared_ptr<SimInput> actInFile = inputs[iFile];

    if( isParentDomain_) { 
      boost::filesystem::path p(actInFile->GetFileName());
      // BOOST PROBLEM!!!
      // if p.leaf() makes problem, p.filename() would maybe work.
      // p.filename() does not compile for me!!
      // What should work:
      // boost::filesystem::base(p) << "." << boost::filesystem::extension(p)
      std::cout << p.leaf() << " " << std::flush;
    }

    actInFile->ReadMesh(actGrid);
  }

  actGrid->FinishInit();

  // Initialize non-conforming interfaces
  if (gridId == "default") {
    // initialize PML layer generation, return if not specified
    actGrid->TriggerAutoLayerGeneration();
    // Initialize non-conforming interfaces, return if not specified
    actGrid->InitNcInterfacesFromXML();
  }


  // check the grid on regularity
  bool regular = false;
  bool non_regular = false;
  for (UInt r=0, numRegs=actGrid->regionData.GetSize(); r < numRegs; ++r)
  {
    if (actGrid->regionData[r].regular)
      regular = true;
    else
      non_regular = true;
  }

  // finish output for grid reading
  if( isParentDomain_) {
    std::cout << "-> ";
    if (regular && non_regular)
      std::cout << "partially ";
    if (!regular && non_regular)
      std::cout << "not ";
    std::cout << "regular" << std::endl; // also regular && !non_regular
  }
  
}

void Domain::PostInit(UInt sequenceStep)
{
  // set up the driver first
  // SetDriver extracts the SingleDriver which is what CreateInstance returns
  // but in the case of a MultiSequenceDriver the SingeDriver is NULL up to init (which is done later!).

  // we do not have to delete driver as it is due to SetDriver() deleted
  // either via ptSingleDriver_ or multiSequenceDriver_ in the destructor
  BaseDriver* driver = BaseDriver::CreateInstance( simState_, this, param_, info_ );
  SetDriver(driver);
  bool restart = isParentDomain_ ? progOpts->GetRestart() : false;

  std::string name = BasePDE::analysisType.ToString(driver->GetAnalysisType());

  shared_ptr<Timer> timer = driver->GetInfo()->Get("init_" + name + "/timer")->AsTimer();

  // the single driver initilized the pde which cannot be done prior initialization of optimization,
  // the multisequence intitialization does not set up the pdes itself and is necessary to init optimization
  if(domain->GetMultiSequenceDriver() != NULL) {
    timer->Start();
    driver->Init(restart);
    timer->Stop();
  }

  // check if we have to do optimization. Do it before driver->Init() to construct the CoefFunctionOpt material
  if(GetParamRoot()->Has("optimization"))
  {
    Optimization::CreateInstance(); // has an SetOptimization() included
  }
  else
  {
    // check if we simulate with ersatz material - after driver and only if not used with optimization
    // if used with optimization loadErsatzMaterial specifies the starting point for optimization
    // and is loaded from Optimization::PostInit because scaling (and EvalObjectiveGradient) is already done before we reach here
    if(DensityFile::NeedLoadErsatzMaterial())
      designSpace_ = DensityFile::ReadErsatzMaterial();
  }

  // For optimization the design needs to be already set to initialize the proper material coefficients
  // in the multisequence case init does something else and was already called above
  // note that the multi sequence driver does not initilize the single pdes yet within Domain::PostInit()
  if(domain->GetMultiSequenceDriver() == NULL) {
    timer->Start();
    driver->Init(restart);
    timer->Stop();
  }

  // we need driver->Init() first
  if(optimization_ != NULL)
  {
    // second initialization phase, constructs material
    optimization_->PostInit();
    // third initialization phase, constructs optimizer
    optimization_->PostInitSecond();
  }

  // note that in the common case isParentDomain_ is true
  if(multiSequenceDriver_ && !isParentDomain_)
    multiSequenceDriver_->SetSequenceStep(sequenceStep);
}

// **************
//   Destructor
// **************
Domain::~Domain()
{
  // delete all grids only, if no external grid were used
  if( !useExternalGridMap_) {
    std::map<std::string, Grid*>::iterator gridIt;
    for (gridIt = gridMap_.begin(); gridIt != gridMap_.end(); gridIt++)
    {
      delete gridIt->second;
      gridIt->second = NULL;
    }
  }

  delete ptIterCoupledPde_;
  ptIterCoupledPde_ = NULL;

  // When the StdVector ptpde_ is destroyed, only the pointers to the PDEs,
  // but not the PDEs themselves will be destroyed
  for (UInt i = 0; i < ptSinglePde_.GetSize(); i++)
  {
    delete (ptSinglePde_[i]);
  }
  ptSinglePde_.Clear();

  // Delete all coordinate systems
  std::map<std::string, CoordSystem*>::iterator it;
  for (it = coordSys_.begin(); it != coordSys_.end(); it++)
  {
    delete (*it).second;
  }
  coordSys_.clear();

  // Delete all direct coupled PDEs
  for (UInt i = 0; i < ptDirectCoupledPde_.GetSize(); i++)
  {
    delete (ptDirectCoupledPde_[i]);
  }
  ptDirectCoupledPde_.Clear();

  // Destructor of IterCoupledPDE deletes couplings!

  // stuff set up by PostInit at the end

  // If our driver is a MultiSequenceDriver we don't delete ptSingleDriver_!
  if (multiSequenceDriver_ != NULL)
  {
    delete multiSequenceDriver_;
    multiSequenceDriver_ = NULL;
  }
  else if (ptSingleDriver_ != NULL)
  {
    delete ptSingleDriver_;
    ptSingleDriver_ = NULL;
  }
  
  if(mathParser_) {
    delete mathParser_;
    mathParser_ = NULL;
  }

  // the optimization is optional. Important, before ersatzMaterial!
  if (optimization_ != NULL) {
    delete optimization_;
    optimization_ = NULL;
  }

  // ersatzMaterial is either set by PostInit()->ReadErsatzMaterial or Optimization
  if(designSpace_ != NULL) {
    delete designSpace_;
    designSpace_ = NULL;
  }

}

void Domain::SolveProblem()
{
  BaseDriver* driver = multiSequenceDriver_;
  if(driver == NULL)
    driver = ptSingleDriver_;

  // PostInit needs to be called in advance!
  if(optimization_ != NULL)
    optimization_->SolveProblem(); // will call multiple driver->SolveProblem
  else
    driver->SolveProblem();
}

  // **********************
  //   Getter function
  // **********************
  StdPDE* Domain::GetStdPDE(const std::string pdeName)
  {
    bool pdeFound = false;
    UInt i;
    std::string errMsg;

    // search the direct coupled pdes
    for (i=0; i<ptDirectCoupledPde_.GetSize(); i++) {
      if ( ptDirectCoupledPde_[i]->GetName() == pdeName ) {
        pdeFound = true;
        break;
      } 
    }

    if ( pdeFound == true)
      return ptDirectCoupledPde_[i];
    else {
      // search the single pdes
      std::map<SinglePDE*,bool>::iterator it;
      
      for (i=0; i<ptSinglePde_.GetSize(); i++) {
        
        if (ptSinglePde_[i]->GetName() == pdeName) {
          
          // check if SinglePDE is not coupled directly
          it = isDirectCoupled_.find( ptSinglePde_[i] );
          
          if ( it == isDirectCoupled_.end() ) {
            EXCEPTION("It was impossible to determine, if the PDE "
                      << "with the name '" << pdeName << "' is directly "
                      << "coupled or not" );
          }
          
          //         if ( (*it).second == false ) {
          //           pdeFound = true;
          //           break;
          //         }
          pdeFound = true;
          break;
        }
      }
   
      
      if (pdeFound == true)
        return ptSinglePde_[i];
      else {
        EXCEPTION( "Domain:GetPDE: PDE with name '" << pdeName
                   << "' was not found/created!." );
        return NULL;
      }
    }
  }




SinglePDE* Domain::GetSinglePDE(const std::string pdeName,  bool throw_exception)
{
  // check for the pde an return
  for (UInt i = 0, s = ptSinglePde_.GetSize(); i < s; ++i)
  {
    if (ptSinglePde_[i]->GetName() == pdeName)
      return ptSinglePde_[i];
  }

  // nothing found
  if (throw_exception)
  {
    EXCEPTION("PDE with name '" << pdeName << "' not found");
    return NULL;
  }
  else
    return NULL;
}

BasePDE* Domain::GetBasePDE()
{

  // if only one SinglePDE exists
  if (numSinglePde_ == 1)
    return ptSinglePde_[0];

  // if only one DirectCouplePDE exists
  else if (numDirectCoupledPde_ == 1 && ptIterCoupledPde_ == NULL)
    return ptDirectCoupledPde_[0];

  // one or more Single/DirectCoupledPDEs
  // are coupled in an iterative way
  else
    return ptIterCoupledPde_;

}

Grid* Domain::GetGrid(const std::string& id)
{
  if (gridMap_.find(id) == gridMap_.end())
    EXCEPTION( "Grid with id '" << id << "' is not defined" );
  return gridMap_[id];
}


CoordSystem* Domain::GetCoordSystem(const std::string& name)
{
  auto it =coordSys_.find(name);

  if (it == coordSys_.end())
    EXCEPTION("A coordinate system with name '" << name << "' is not registered!" );

  return (*it).second;
}

StdVector<std::string> Domain::GetCoordSystems() const
{
  StdVector<std::string> res;

  for(std::map<std::string, CoordSystem*>::const_iterator it = coordSys_.begin(); it != coordSys_.end(); ++it)
    res.Push_back(it->first);

  return res;
}

// **************************
//   Initialization of PDEs
// **************************
void Domain::CreatePDEs(UInt sequenceStep, PtrParamNode infoNode)
{
  // create single pde(s)
  CreateSinglePDEs(sequenceStep, infoNode);

  // create direct coupled pde(s)
  CreateDirectCoupledPDEs(sequenceStep, infoNode);

  // Create iterative coupled pde
  CreateIterCoupledPDE(sequenceStep, infoNode);
}

void Domain::RestorePDEs(StdVector<SinglePDE*>& single)
{
  ptSinglePde_ = single;
  assert(ptSinglePde_.GetSize() == single.GetSize());
  numSinglePde_ = single.GetSize();

  // restoring direct coupled and iterative coupled not yet implemented
  isDirectCoupled_.clear();
  for(unsigned int i = 0; i < single.GetSize(); i++)
    isDirectCoupled_[single[i]] = false;
}


void Domain::InitPDEs(UInt sequenceStep)
{
  // in case we have an iterative coupled PDE,
  // we take its info pointer and use it
  // as base for the coupled ones
  PtrParamNode base;
  if (ptIterCoupledPde_) {
    base = ptIterCoupledPde_->GetInfoNode();
  }

  // Initialize those PDEs which are not directly coupled
  std::map<SinglePDE*, bool>::iterator it;

  for( UInt iStage = 0; iStage < 1; ++iStage ) {
    for (UInt i = 0; i < numSinglePde_; i++) {
      it = isDirectCoupled_.find(ptSinglePde_[i]);
      if ((*it).second == false) {
        switch(iStage) {
          case 0:
            ptSinglePde_[i]->Init_Stage1(sequenceStep,base);
            break;
          case 1:
            ptSinglePde_[i]->Init_Stage2();
            break;
          case 2:
            ptSinglePde_[i]->Init_Stage3();
            break;
          default:
            EXCEPTION( "Only 3 stages of initialization known");
            break;
        }
      }
    }
  }

  // initialize direct coupled pde(s)
  // -> this triggers also the initialization of
  // those single PDEs which are directly coupled
  for (UInt i = 0; i < numDirectCoupledPde_; i++)
  {
    if( isParentDomain_) {
	std::cout << "++ Initializing direct coupling" << std::endl;
	}
	//std::cout << "Domain.cc - preInit: pde->Name()? " << ptDirectCoupledPde_[i]->GetName() << std::endl;
	//std::cout << "Domain.cc - preInit: pde->IsNonLin()? " << ptDirectCoupledPde_[i]->IsNonLin() << std::endl;
	
    ptDirectCoupledPde_[i]->Init(sequenceStep);
    ptDirectCoupledPde_[i]->DefineAlgSys();
    
    //std::cout << "Domain.cc - postInit: pde->IsNonLin()? " << ptDirectCoupledPde_[i]->IsNonLin() << std::endl;
  }

  
  // Initialize those PDEs which are not directly coupled
  for( UInt iStage = 1; iStage < 3; ++iStage ) {
    for (UInt i = 0; i < numSinglePde_; i++) {
      it = isDirectCoupled_.find(ptSinglePde_[i]);
      if ((*it).second == false) {
        switch(iStage) {
          case 0:
            ptSinglePde_[i]->Init_Stage1(sequenceStep,base);
            break;
          case 1:
            ptSinglePde_[i]->Init_Stage2();
            break;
          case 2:
            ptSinglePde_[i]->Init_Stage3();
            break;
          default:
            EXCEPTION( "Only 3 stages of initialization known");
            break;
        }
      }
    }
  }


  // Initialize algebraic system of each SinglePDE
  // Note: DefineAlgSys() triggers only the initialization
  // of those SinglePDEs, which are not directly coupled
  for (UInt i = 0; i < numSinglePde_; i++)
  {
    it = isDirectCoupled_.find(ptSinglePde_[i]);
    if ((*it).second == false)
    {
      ptSinglePde_[i]->DefineAlgSys();
    }
  }
}

// **************************
//   Creation of PDEs
// **************************
void Domain::CreateSinglePDEs(UInt sequenceStep, PtrParamNode infoNode)
{
  // default grid
  Grid* defaultGrid = gridMap_["default"];

  ParamNodeList pdeNodes = param_->GetByVal("sequenceStep", std::string("index"), sequenceStep)->Get("pdeList")->GetChildren();

  ptSinglePde_.Resize(pdeNodes.GetSize());
  ptSinglePde_.Init();
  numSinglePde_ = pdeNodes.GetSize();
  
  for (UInt i = 0; i < pdeNodes.GetSize(); i++)
  {

    std::string actPdeName = pdeNodes[i]->GetName();
    PtrParamNode actPdeNode = pdeNodes[i];
    if( isParentDomain_) 
      std::cout << "++ Creating PDE '" + actPdeName + "' for analysis '"
                << BasePDE::analysisType.ToString(domain->GetSingleDriver()->GetAnalysisType()) << "'" << std::endl;

    if (actPdeName == "electrostatic")
      ptSinglePde_[i] = new ElecPDE(defaultGrid, actPdeNode, infoNode, simState_, this);

    else if (actPdeName == "mechanic")
      ptSinglePde_[i] = new MechPDE(defaultGrid, actPdeNode, infoNode, simState_, this);

    else if (actPdeName == "acoustic") {
        ptSinglePde_[i] = new AcousticPDE(defaultGrid, actPdeNode, infoNode, simState_, this);
    }
    // else if (actPdeName == "acoustic_BEM") {
    //     ptSinglePde_[i] = new AcousticPDE_BEM(/* defaultGrid, actPdeNode, infoNode, simState_, this */);
    // }
    else if (actPdeName == "split") {
      ptSinglePde_[i] = new AcousticSplitPDE(defaultGrid, actPdeNode, infoNode,
          simState_, this );
    }
    else if (actPdeName == "acousticMixed")
      ptSinglePde_[i] = new AcousticMixedPDE(defaultGrid, actPdeNode, infoNode, simState_, this);

    else if (actPdeName == "magnetic"){
      std::string formulation = actPdeNode->Get("formulation")->As<std::string>();
      if (formulation == "A") {
        ptSinglePde_[i] = new MagneticPDE(defaultGrid, actPdeNode, infoNode, simState_, this);
      }else if(formulation == "Psi") {
        ptSinglePde_[i] = new MagneticScalarPotentialPDE(defaultGrid, actPdeNode, infoNode, simState_, this);
      }else{
        EXCEPTION("Formulation of MagEdgePDE not known!");
      }
    }

    else if (actPdeName == "magneticEdge"){
      std::string formulation = actPdeNode->Get("formulation")->As<std::string>();
      if (formulation == "A") {
        ptSinglePde_[i] = new MagEdgePDE(defaultGrid, actPdeNode, infoNode, simState_, this);
      }else{
        if(formulation == "A-V"){
          ptSinglePde_[i] = new MagEdgeMixedAVPDE(defaultGrid, actPdeNode, infoNode, simState_, this);
        }else if(formulation == "specialA-V"){
          ptSinglePde_[i] = new MagEdgeSpecialAVPDE(defaultGrid, actPdeNode, infoNode, simState_, this);
        }else if(formulation == "Darwin" || "Darwin_doubleLagrange"){
          ptSinglePde_[i] = new DarwinPDE(defaultGrid, actPdeNode, infoNode, simState_, this);
        }else{
          EXCEPTION("Formulation of MagEdgePDE not known!");
        }
      }
    }

    else if (actPdeName == "heatConduction")
      ptSinglePde_[i] = new HeatPDE(defaultGrid, actPdeNode, infoNode, simState_, this);

    else if (actPdeName == "fluidMechLin")
      ptSinglePde_[i] = new LinFlowPDE(defaultGrid, actPdeNode, infoNode, simState_, this);

    else if (actPdeName == "fluidMech") {
      std::string formulation = actPdeNode->Get("formulation")->As<std::string>();

      if (formulation == "perturbed") {
        ptSinglePde_[i] = new PerturbedFlowPDE(defaultGrid, actPdeNode, infoNode, simState_, this);
      }
      else {
        ptSinglePde_[i] = new FlowPDE(defaultGrid, actPdeNode, infoNode, simState_, this);
      }
    }
    else if (actPdeName == "testPDE") {
        ptSinglePde_[i] = new TestPDE(defaultGrid, actPdeNode, infoNode, simState_, this);
    }
    else if (actPdeName == "elecConduction") {
        ptSinglePde_[i] = new ElecCurrentPDE(defaultGrid, actPdeNode, infoNode, simState_, this);
    }
    else if (actPdeName == "elecQuasistatic") {
        ptSinglePde_[i] = new ElecQuasistaticPDE(defaultGrid, actPdeNode, infoNode, simState_, this);
    }
    else if (actPdeName == "waterWave") {
      ptSinglePde_[i] = new WaterWavePDE(defaultGrid, actPdeNode, infoNode, simState_, this);
    }
    else if (actPdeName == "LatticeBoltzmann") {
        ptSinglePde_[i] = new LatticeBoltzmannPDE(defaultGrid, actPdeNode, infoNode, simState_, this);
    }
    else if (actPdeName == "smooth") {
        ptSinglePde_[i] = new SmoothPDE(defaultGrid, actPdeNode, infoNode, simState_, this);
    }
    else
    {
      EXCEPTION( actPdeName << " - this type of pdes is unknown" );
    }

    // by default, not single pde is directly coupled
    isDirectCoupled_[ptSinglePde_[i]] = false;
    
    // Ensure, that at least one PDE is present
    if( numSinglePde_ == 0 ) {
      WARN( "No PDEs were created. Please check your parameter file.")
    }
    // Initialize current PDE
    // -> This step has now moved to method InitPDEs
    //ptSinglePde_[i]->Init();
  }
}

void Domain::CreateIterCoupledPDE(UInt sequenceStep, PtrParamNode infoNode)
{

  std::string errMsg;

  // check if more than one PDEs are defined
  if (numIterCoupledStdPde_ <= 1) {
    ptIterCoupledPde_ = NULL;
    return;
  }

  if( isParentDomain_) 
    std::cout << "++ Creating coupling" << std::endl;

  // ================================
  //   Check for iterative coupling
  // ================================

  // check for presence of "couplingList" and "iterative" element
  PtrParamNode iterNode;
  PtrParamNode couplingNode =
      param_->GetByVal("sequenceStep", std::string("index"), sequenceStep) 
      ->Get("couplingList", ParamNode::PASS);
  if ( couplingNode) {
    iterNode = couplingNode->Get("iterative", ParamNode::PASS);
  }

  // Create IterCoupledPDE and pass all StdPDEs to it
  ptIterCoupledPde_ = new IterCoupledPDE( ptSinglePde_,
                                          ptDirectCoupledPde_,
                                          iterNode, infoNode, simState_, this  );
  
  // Loop over all SinglePDEs and pass pointer to iterative coupled PDE
  for( UInt i = 0; i < ptSinglePde_.GetSize(); ++i ) {
    // std::cout << "PDE: " << ptSinglePde_[i]->GetName() << std::endl;
    ptSinglePde_[i]->SetIterCoupledPDE( ptIterCoupledPde_ );
  }
  
}

// ***************************
//   CreateDirectCoupledPDEs
// ***************************
void Domain::CreateDirectCoupledPDEs(UInt sequenceStep, PtrParamNode infoNode)
{

  numIterCoupledStdPde_ = numSinglePde_;

  // if only one pde exists, there is nothing to do
  if (numSinglePde_ <= 1)
  {
    return;
  }

  // get "couplingList" node (must exist)
  PtrParamNode couplingNode =  param_->GetByVal("sequenceStep", std::string("index"), sequenceStep)->Get("couplingList");
  PtrParamNode directNode = couplingNode->Get("direct", ParamNode::PASS);
  if (!directNode)
    return;

  // get nodes of pairwise direct couplings
  ParamNodeList pairNodes = directNode->GetChildren();

  SinglePDE *pde1 = NULL;
  SinglePDE *pde2 = NULL;
  BasePairCoupling *coupling = NULL;

  // HARD CODED: At the moment we allow only one direct coupled pde
  // with only on pairwise coupling
  StdVector<SinglePDE*> singlePdes;
  StdVector<BasePairCoupling*> DirectCouplingPairs;
  std::set<std::string> setSinglePDEs;
  std::string couplingName;

  for (UInt i = 0; i < pairNodes.GetSize(); i++)
  {

    // get couplingName
    couplingName = pairNodes[i]->GetName();

    // *** PIEZO Coupling ***
    if (couplingName == "piezoDirect")
    {

      pde1 = GetSinglePDE("mechanic");
      pde2 = GetSinglePDE("electrostatic");

      // in the case of piezo coupling, the electrostatic
      // entries have to be multiplied by -1
      dynamic_cast<ElecPDE*> (pde2)->SetPiezoCoupling();

      coupling = new PiezoCoupling(pde1, pde2, pairNodes[i], info_,
                                   simState_, this );

    }
//    // *** magnetostriction Coupling *** 
//    else if ( couplingName == "magnetoStrictionDirect" ) {
//
//      pde1 = GetSinglePDE( "mechanic" );
//      pde2 = GetSinglePDE( "magneticScalar" );
//
//      // in the case of  Mag-Mech coupling, the magnetic field
//      // entries have to be multiplied by -1
//      dynamic_cast<MagScalarPDE*>(pde2)->SetMagStrictCoupling();
//      coupling = new MagStrictCoupling( pde1, pde2, pairNodes[i] );
//
//    }
//    
    // *** ACOU-MECH Coupling ***
    else if (couplingName == "acouMechDirect")
    {

      pde1 = GetSinglePDE("mechanic");
      pde2 = GetSinglePDE("acoustic");

      // in the case of acou-Mech coupling, the acoustic
      // entries have to be multiplied by -1
      dynamic_cast<AcousticPDE*> (pde2)->SetMechanicCoupling();

      coupling = new AcouMechCoupling(pde1, pde2, pairNodes[i], info_,
                                      simState_, this );
    }
    
    else if (couplingName == "magnetoStrictDirect")
    {

      pde1 = GetSinglePDE("mechanic");
      pde2 = GetSinglePDE("magnetic");

	//pass mechanic pde to magnetic for the case of nonlinear magnetostriction
      dynamic_cast<MagneticPDE*> (pde2)->SetMagnetoStrictCoupling(pde1);

      coupling = new MagnetoStrictCoupling(pde1, pde2, pairNodes[i], info_,
                                      simState_, this );
    }
    
    // *** FLUID-MECH Coupling ***
    else if (couplingName == "fluidMechDirect")
    {

      pde1 = GetSinglePDE("mechanic");
      pde2 = GetSinglePDE("fluidMechPerturbed");

      // in the case of acou-Mech coupling, the acoustic
      // entries have to be multiplied by -1
//      dynamic_cast<PerturbedFlowPDE*> (pde2)->SetMechanicCoupling();

      coupling = new FluidMechCoupling(pde1, pde2, pairNodes[i], info_,
                                       simState_, this );
    }
    // *** WATER WAVE-ACOUSTIC Coupling ***
    else if (couplingName == "waterWaveAcouDirect")
    {
      pde1 = GetSinglePDE("waterWave");
      pde2 = GetSinglePDE("acoustic");

      coupling = new WaterWaveAcousticCoupling(pde1, pde2, pairNodes[i], info_,
                                               simState_, this );
    }
    // *** Water Wave-MECH Coupling ***
    else if (couplingName == "waterWaveMechDirect")
    {

      pde1 = GetSinglePDE("mechanic");
      pde2 = GetSinglePDE("waterWave");

      coupling = new WaterWaveMechCoupling(pde1, pde2, pairNodes[i], info_,
                                           simState_, this );
    }
    // *** Linear flow coupled with heat ***
    else if (couplingName == "linFlowHeatDirect")
    {

      pde1 = GetSinglePDE("fluidMechLin");
      pde2 = GetSinglePDE("heatConduction");

      coupling = new LinFlowHeatCoupling(pde1, pde2, pairNodes[i], info_,
    		                             simState_, this );
    }
    // *** Linear flow coupled with acoustic wave equation***
    else if (couplingName == "linFlowAcouDirect")
    {

      pde1 = GetSinglePDE("fluidMechLin");
      pde2 = GetSinglePDE("acoustic");

      coupling = new LinFlowAcouCoupling(pde1, pde2, pairNodes[i], info_,
    		                             simState_, this );
    }
    else if (couplingName == "linFlowMechDirect")
    {

      pde1 = GetSinglePDE("fluidMechLin");
      pde2 = GetSinglePDE("mechanic");

      coupling = new LinFlowMechCoupling(pde1, pde2, pairNodes[i], info_,
                                     simState_, this );
    }
//    // ------------------------------------------------------------------------
//    // *** THERMO-MECH Coupling ***
//    else if (couplingName == "thermoMechDirect")
//    {
//
//      pde1 = GetSinglePDE("mechanic");
//      pde2 = GetSinglePDE("heatConduction");
//
//      //turn the coupling on in heat-conduction pde
//      //in order to create the heat-conduction damp matrix
//      dynamic_cast<HeatCondPDE*> (pde2)->SetMechCoupling();
//
//      coupling = new ThermoMechCoupling(pde1, pde2, pairNodes[i]);
//    }
//    // *** THERMO-ELECTRIC Coupling ***
//    else if (couplingName == "thermoElectricDirect")
//    {
//
//      pde1 = GetSinglePDE("electrostatic");
//      pde2 = GetSinglePDE("heatConduction");
//
//      //turn the coupling on in heat-conduction pde
//      //in order to create the heat-conduction damp matrix
//      dynamic_cast<HeatCondPDE*> (pde2)->SetElectroCoupling();
//
//      // in the case of thermo coupling, the electrotstatic
//      // entries have to be multiplied by -1
//      dynamic_cast<ElecPDE*> (pde1)->SetThermoCoupling();
//
//      coupling = new ThermoElectricCoupling(pde1, pde2, pairNodes[i]);
//    }
//    //------------------------------------------------------------------------
//
    else
    {
      EXCEPTION( "The direct coupling '" << couplingName
          << "' is not implemented!" << std::endl );
    }

    // set flag for direct coupling
    isDirectCoupled_[pde1] = true;
    isDirectCoupled_[pde2] = true;

    // add single PDEs and couplings into collections
    setSinglePDEs.insert(pde1->GetName());
    setSinglePDEs.insert(pde2->GetName());
    DirectCouplingPairs.Push_back(coupling);
  }

  // check if any pair coupling was found
  if (coupling == NULL)
    return;

  // Transform set of PDEs into a vector
  std::set<std::string>::iterator itSet;

  for (itSet = setSinglePDEs.begin(); itSet != setSinglePDEs.end(); itSet++)
  {
    singlePdes.Push_back(GetSinglePDE(*itSet));
  }

  ptDirectCoupledPde_.Push_back(new DirectCoupledPDE(gridMap_["default"], 
                                                     PtrParamNode(), info_,simState_,
                                                     this ));
  ptDirectCoupledPde_[0]->SetSinglePDEs(singlePdes);
  ptDirectCoupledPde_[0]->SetCouplings(DirectCouplingPairs);

  // At the moment we allow only one direct coupled pde, so we set the
  // number of direct coupledPDEs to one;
  numDirectCoupledPde_ = 1;

  // now determine, how many SinglePDEs are coupling directly
  std::map<SinglePDE*, bool>::iterator it = isDirectCoupled_.begin();

  numIterCoupledStdPde_ = numDirectCoupledPde_;

  while (it != isDirectCoupled_.end())
  {
    if ((*it).second == false)
      numIterCoupledStdPde_++;
    it++;
  }

}

void Domain::CreateCoordinateSystems()
{
  PtrParamNode in = info_->Get(ParamNode::HEADER)->Get("domain");
  in = in->Get("coordinateSystems", ParamNode::APPEND);
  
  // first create the "global" standard cartesian
  // coordinate system
  std::string defaultname = "default";
  coordSys_[defaultname] = new DefaultCoordSystem(gridMap_["default"] );

  // get nodes with coordinate systems
  PtrParamNode coosyNode = param_->Get("domain")->Get("coordSysList", ParamNode::PASS);
  if (!coosyNode)
    return;

  ParamNodeList coordNodes = coosyNode->GetChildren();

  // iterate over all coordinate system nodes
  for (UInt i = 0; i < coordNodes.GetSize(); i++)
  {

    // fetch coosy name and type
    std::string type = coordNodes[i]->GetName();
    std::string name = coordNodes[i]->Get("id")->As<std::string>();

    CoordSystem * actCoord = NULL;
    if (type == "polar")
    {
      actCoord = new PolarCoordSystem( name, gridMap_["default"], 
                                       coordNodes[i] );
    }
    else if (type == "cylindric")
    {
      actCoord = new CylCoordSystem( name, gridMap_["default"], 
                                     coordNodes[i] );
    }
    else if (type == "cartesian")
    {
      actCoord = new CartesianCoordSystem(name, gridMap_["default"],
          coordNodes[i] );
    }
    else if (type == "trivialCartesian")
    {
      actCoord = new TrivialCartesianCoordSystem(name, gridMap_["default"],
                                                 coordNodes[i] );
    }
    else
    {
      EXCEPTION( "Coordinate system with type '" << type
          << "' not known!" );
    }
    coordSys_[name] = actCoord;
    actCoord->ToInfo(in);
  }
}

void Domain::RegisterVariables() 
{
  PtrParamNode varListNode = param_->Get("domain")->Get("variableList", ParamNode::PASS);
  if( varListNode ) {
   ParamNodeList & varNodes = varListNode->GetChildren();
   ParamNodeList::iterator it = varNodes.Begin();
   std::string varName, valString;
   unsigned int handle;
   handle = mathParser_->GetNewHandle();
   Double value = 0.0;
   for(; it != varNodes.End(); it++ ) {
     (*it)->GetValue("name", varName);
     (*it)->GetValue("value", valString);
     // check for reserved variable names
     if ( (varName == "t") || (varName == "dt") || (varName == "f") || (varName == "step") )
     {
       EXCEPTION("The variable '" << varName
                 << "' is reserved, its value will be set automatically. "
                 << "Please choose a different name.");
     }
     mathParser_->SetExpr(handle, valString);
     value = mathParser_->Eval(handle);
     mathParser_->SetValue(MathParser::GLOB_HANDLER, varName, value);
   }
   mathParser_->ReleaseHandle(handle);
  }
}
// *************
//   SetDriver
// *************
void Domain::SetDriver(BaseDriver * driver)
{
  if (driver->GetAnalysisType() == BasePDE::MULTI_SEQUENCE)
  {
    multiSequenceDriver_ = dynamic_cast<MultiSequenceDriver*>(driver);
    ptSingleDriver_ = multiSequenceDriver_->GetSingleDriver(); // NULL before Init()!!
  }
  else
  {
    ptSingleDriver_ = dynamic_cast<SingleDriver*>(driver);
  }
}

BaseDriver* Domain::GetDriver()
{
  return ptSingleDriver_;
}

// *************
//   ResetPDEs
// *************
void Domain::ResetPDEs(bool keep)
{
  // keep is for optimization with multiple sequences.
  // Then the MultiSequenceDriver keeps drivers and pdes instead of deleting then with a new sequence step
  // Delete single pde(s)
  for(UInt iPDE = 0; iPDE < numSinglePde_; iPDE++)
    if(!keep)
      delete ptSinglePde_[iPDE];

  ptSinglePde_.Clear(); // for keep the drivers are alread in MultiSequenceDriver::keptPDEs_

  // delete direct coupled pde(s)
  for (UInt iPDE = 0; iPDE < numDirectCoupledPde_; iPDE++) {
    assert(!keep); // not yet implemented
    delete ptDirectCoupledPde_[iPDE];
  }
  ptDirectCoupledPde_.Clear();

  // delete iterative coupled pde
  if (ptIterCoupledPde_ != NULL) {
    assert(!keep); // not yet implemented
    delete ptIterCoupledPde_;
  }
  ptIterCoupledPde_ = NULL;

  // Also reset all state variables
  isDirectCoupled_.clear();
  numSinglePde_ = 0;
  numDirectCoupledPde_ = 0;
  numIterCoupledStdPde_ = 0;
}

// *************
//   SetGridMap
// *************
void Domain::SetGridMap( const std::map<std::string, Grid* >& gridMap ) {
  gridMap_ = gridMap;
  useExternalGridMap_ = true;
}


// *************
//   PrintGrid
// *************
void Domain::PrintGrid()
{

  resultHandler_->Init(gridMap_, true);
  resultHandler_->Finalize();
}

void Domain::Dump()
{
  for (UInt i = 0; i < ptSinglePde_.GetSize(); i++)
  {
    std::cout << "  single pde: " << ptSinglePde_[i]->GetName() << std::endl;
  }

  for (UInt i = 0; i < ptDirectCoupledPde_.GetSize(); i++)
  {
    std::cout << "  direct coupled pde: " << ptDirectCoupledPde_[i]->GetName()
        << std::endl;
  }

}


void Domain::ToInfo(PtrParamNode in, bool force_default)
{
  PtrParamNode in_ = in->Get("coordinateSystems");
  for(std::map<std::string, CoordSystem*>::iterator it = coordSys_.begin(); it != coordSys_.end(); ++it)
  {
    PtrParamNode s = in_->Get("system", ParamNode::APPEND);
    s->Get("name")->SetValue(it->first);
    it->second->ToInfo(s);
  }

  // CalcGridBoundingBox() takes some time and most people are not interested in the information
  if(force_default || progOpts->DoDetailedInfo())
  {
    // Get bounding box of the grid
    Matrix<double> m = this->GetGrid()->CalcGridBoundingBox();
    PtrParamNode s = in_->Get("domain", ParamNode::APPEND);
    s->Get("min_x")->SetValue(m[0][0]);
    s->Get("max_x")->SetValue(m[0][1]);
    s->Get("min_y")->SetValue(m[1][0]);
    s->Get("max_y")->SetValue(m[1][1]);
    if (m.GetNumRows() > 2) {
      s->Get("min_z")->SetValue(m[2][0]);
      s->Get("max_z")->SetValue(m[2][1]);
    }
  }

}

bool Domain::HasPerdiodicBC() const
{
  for(unsigned int i = 0; i < ptSinglePde_.GetSize(); i++)
  {
    std::map<SolutionType, shared_ptr<BaseFeFunction> > fes = ptSinglePde_[i]->GetFeFunctions();
    for(std::map<SolutionType, shared_ptr<BaseFeFunction> >::const_iterator it = fes.begin(); it != fes.end(); it++)
      if(it->second->HasPeriodicBC())
        return true;
  }
  return false;
}

}
