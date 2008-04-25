// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "domain.hh"

#include <set>
#include <memory>

#include "General/environment.hh"
#include "General/exception.hh"
#include "Domain/grid.hh"
#include "Domain/GridCFS/grid_cfs.hh"
#include "DataInOut/simInput.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/MaterialHandler.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/Xerces.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/simInput.hh"
#include "General/exception.hh"
#include "Utils/coordSystem.hh"
#include "Utils/cylCoordSys.hh"
#include "Utils/polCoordSys.hh"
#include "Utils/cartesianCoordSys.hh"
#include "Utils/trivialCartesianCoordSys.hh"
#include "Utils/defaultCoordSys.hh"
#include "DataInOut/Scripting/cfsmessenger.hh"
#include "DataInOut/resultHandler.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/DesignSpace.hh"
#include "PDE/pdes_header.hh"
#include "PDE/basePDE.hh"
#include "Utils/result.hh"

// Coupling of single PDEs
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "CoupledPDE/itercoupledpde.hh"
#include "CoupledPDE/coupledpdedef.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "CoupledPDE/PiezoCoupling.hh"
#include "CoupledPDE/AcouMechCoupling.hh"
#include "CoupledPDE/ThermoMechCoupling.hh"
#include "CoupledPDE/ThermoElectricCoupling.hh"

// Include driver
#include "Driver/basedriver.hh"
#include "Driver/singleDriver.hh"
#include "Driver/multiSequenceDriver.hh"
#ifdef ADAPTGRID
#include "domain/AdaptGrid/interface_adgrid.hh"
#endif



namespace CoupledField {

  
  // **************
  //   Construtor
  // **************

  Domain::Domain( std::map<std::string, StdVector<shared_ptr<SimInput> > >& gridInputs,
                  ResultHandler * handler,
                  MaterialHandler * ptMat ) {
    
    
    
    Info->PrintF("","==================\n");
    Info->PrintF("","   DOMAIN SETUP   \n");
    Info->PrintF("","==================\n\n");
    
    // initialize data
    numSinglePde_ = 0;
    numDirectCoupledPde_ = 0;
    numIterCoupledStdPde_ = 0;
    
    // assign pointers
    gridInputs_ = gridInputs;
    ptMatHandler_ = ptMat;
    resultHandler_  = handler;
    ptIterCoupledPde_ = NULL;
    ptSingleDriver_ = NULL;
    multiSequenceDriver_ = NULL;
    optimization_ = NULL;
    ersatzMaterial = NULL;
  }
  
  void Domain::CreateGrid() {
   
    // read type of mesh library
    std::string libmesh = "cfsGrid";
    ParamNode * inputNode = param->Get("input", false );
    if( inputNode )
      inputNode->Get("meshLibrary", libmesh );
    
    std::string probGeo;
    param->Get("domain")->Get( "geometryType", probGeo );
    if ( probGeo == "3d") {
      dim_ = 3;
    } 
    else if ( probGeo == "axi" || probGeo == "plane" ) {
      dim_ = 2;
    }
    else {
      EXCEPTION( "Wrong Dimension parameter in xml-File" );
    }
    
    Info->StartProgress("Reading in the mesh");
    
    std::map<std::string, StdVector<shared_ptr<SimInput> > >
      ::const_iterator gridIt;

    // iterate over all grid ids
    for( gridIt = gridInputs_.begin(); gridIt != gridInputs_.end(); gridIt++ ) {

      std::string gridId = gridIt->first;
      StdVector<shared_ptr<SimInput> > const & inputs = gridIt->second;

      // iterate over all inputs for the current grid and init reader
      for( UInt iFile = 0; iFile < inputs.GetSize(); iFile++ ) {
        shared_ptr<SimInput> actInFile = inputs[iFile];
        actInFile->InitModule();
      }

      // create new grid
      Grid * actGrid = NULL;
      if (libmesh =="cfsGrid") {
        if(gridId == "default")
          actGrid = new GridCFS(dim_);
        else
        {
          actGrid = new GridCFS(inputs[0]->GetDim());
        }
      } else {
        EXCEPTION( "Type of mesh_library should be one of "
                   << "'cfsgrid' or 'adaptgrid', but is '" << libmesh << "'" );
      }
      
      // add grid to internal map
      gridMap_[gridId] = actGrid;

      // Read in coordinate systems
      // This has to be done, before the grid gets finalized,
      // as some methods in the grid rely on the existence of coordinate
      // systems
      if(gridId == "default")
        CreateCoordinateSystems();
      
      // iterate over all inputs for the current grid and read mesh
      for( UInt iFile = 0; iFile < inputs.GetSize(); iFile++ ) {
        shared_ptr<SimInput> actInFile = inputs[iFile];
        actInFile->ReadMesh(actGrid);
      }

      actGrid->FinishInit();
    }
      
    SETPROFILE("After Grid Creation");

    Info->FinishProgress();

    Info->StartProgress("Initializing non-matching interfaces");
    
    // Call the nonmatching grid intersection calculation
    gridMap_["default"]->InitNonmatchingInterfaces();

    Info->PrintF("","\n=========================\n");
    Info->PrintF("","   END OF DOMAIN SETUP   \n");
    Info->PrintF("","=========================\n\n");
    Info->FinishProgress();

    if ( !commandLine->GetPrintGrid() == true ) {
      // Initialize resultHandler
      resultHandler_->Init( gridMap_["default"], false );
    }   

  }

  void Domain::PostInit() 
  {
      // set up the driver first
      // SetDriver extracts the SingleDriver which is what CreateInstance returns
      // but in the case of an MultiSequenceDriver the SingeDriver is NULL up to init.
      
      // we do not have to delete driver as it is due to SetDriver() deleted
      // either via ptSingleDriver_ or multiSequenceDriver_ in the destructor
      BaseDriver* driver = BaseDriver::CreateInstance();  
      SetDriver(driver); // see above!
      Info->FinishProgress();

      // initialize the driver
      driver->Init();

      // check if we simulate with ersatz material - after driver!! 
      //- not when we do optimization and vice versa
      bool esm = false;
      if(param->Has("loadErsatzMaterial"))
      {
        ReadErsatzMaterial(param->Get("loadErsatzMaterial"));
        esm = true;
      }
      // check if we have to do optimization
      if(param->Has("optimization"))
      {
        if(esm) throw Exception("you cannot load ErsatzMaterial and to optimization");
        SetOptimization(Optimization::CreateInstance());
      }
      else SetOptimization(NULL);
  }


  // **************
  //   Destructor
  // **************
  Domain::~Domain() {


    // delete all grid
    std::map<std::string, Grid* >::iterator gridIt;
    for( gridIt = gridMap_.begin(); gridIt != gridMap_.end(); gridIt++ ) {
      delete gridIt->second;
      gridIt->second = NULL;
    }

    delete ptIterCoupledPde_;
    ptIterCoupledPde_ = NULL;

    // When the StdVector ptpde_ is destroyed, only the pointers to the PDEs,
    // but not the PDEs themselves will be destroyed
    for ( UInt i = 0; i < ptSinglePde_.GetSize(); i++ ) {
      delete (ptSinglePde_[i]);
    }
    ptSinglePde_.Clear();

    // Delete all coordinate systems
    std::map<std::string, CoordSystem*>::iterator it;
    for ( it=coordSys_.begin(); it!=coordSys_.end(); it++ ) {
      delete (*it).second;
    }
    coordSys_.clear();

    // Delete all direct coupled PDEs
    for ( UInt i = 0; i < ptDirectCoupledPde_.GetSize(); i++ ) {
      delete (ptDirectCoupledPde_[i]);
    }
    ptDirectCoupledPde_.Clear();

    // Destructor of IterCoupledPDE deletes couplings!
    
    // stuff set up by PostInit at the end
    
    // If our driver is a MultiSequenceDriver we don't delete ptSingleDriver_!
    if(multiSequenceDriver_ != NULL) { delete multiSequenceDriver_; multiSequenceDriver_ = NULL; }
    else if(ptSingleDriver_ != NULL) { delete ptSingleDriver_; ptSingleDriver_ = NULL; }
    
    // the optimization is optional. Important, before ersatzMaterial!
    if(optimization_ != NULL) { delete optimization_; optimization_ = NULL; }
    
    // ersatzMaterial is either set by PostInit()->ReadErsatzMaterial or Optimization 
    if(ersatzMaterial != NULL) { delete ersatzMaterial; ersatzMaterial = NULL; }    
  }


  void Domain::SolveProblem() 
  {
    BaseDriver* driver = multiSequenceDriver_;
    if(driver == NULL) driver = ptSingleDriver_;
    
    // PostInit needs to be called in advance!
    if(GetOptimization() != NULL) 
       GetOptimization()->SolveProblem(); // will call multiple driver-SolveProblem
    else
       driver->SolveProblem();
  }


  // **********************
  //   Getter function
  // **********************
  StdPDE * Domain::GetStdPDE(const std::string pdeName)
  {
    bool pdeFound = false;
    UInt i;
    std::string errMsg;

    // search the direct coupled pdes

    // *** IMPLEMENT ***

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
    }
  }
  
  SinglePDE * Domain::GetSinglePDE(const std::string pdeName)
  {
    bool pdeFound = false;
    UInt i;
    std::string errMsg;

    for (i=0; i<ptSinglePde_.GetSize(); i++) {
      if (ptSinglePde_[i]->GetName() == pdeName) {
        pdeFound = true;
        break;
      }
    }
    
    if (pdeFound == true)
      return ptSinglePde_[i];
    else {
      EXCEPTION( "Domain:GetSinglePDE: PDE with name '" << pdeName
                 << "' was not found/created!." );
    }
  }

  BasePDE* Domain::GetBasePDE()
  {
    
    // if only one SinglePDE exists
    if ( numSinglePde_ == 1)
      return ptSinglePde_[0];

    // if only one DirectCouplePDE exists
    else if ( numDirectCoupledPde_ == 1 &&
              ptIterCoupledPde_ == NULL )
      return ptDirectCoupledPde_[0];

    // one or more Single/DirectCoupledPDEs
    // are coupled in an iterative way
    else 
      return ptIterCoupledPde_;
    
  }

  Grid * Domain::GetGrid( const std::string& id ) {
    if( gridMap_.find( id ) == gridMap_.end() ) {
      EXCEPTION( "Grid with id '" << id 
                 << "' is not defined" );
    }
    return gridMap_[id];    
  }

  CoordSystem * Domain::GetCoordSystem( const std::string & name ) {
    
    
    std::map<std::string, CoordSystem*>::iterator it;

    it = coordSys_.find(name);

    if ( it == coordSys_.end() ) {
      EXCEPTION( "Domain::GetCoordSystem: A coordinate system with name '"
                 << name << "' is not registered!" );
    }

    return (*it).second;
  
  }


  // **************************
  //   Initialization of PDEs
  // **************************
  void Domain::CreatePDEs( UInt sequenceStep ) {


    // create single pde(s)
    CreateSinglePDEs( sequenceStep );    
    
    // create direct coupled pde(s)
    CreateDirectCoupledPDEs( sequenceStep );

    // Create iterative coupled pde
    CreateIterCoupledPDE( sequenceStep );
  }

  void Domain::InitPDEs( UInt sequenceStep ) {
    // intialize single pde(s)

    // Initialize those PDEs which are not
    // directly coupled
    std::map<SinglePDE*,bool>::iterator it;

    for (UInt i=0; i<numSinglePde_; i++) {
      it = isDirectCoupled_.find( ptSinglePde_[i] );
      if ( (*it).second == false) {
        ptSinglePde_[i]->Init( sequenceStep );
      }
    }

    // initialize direct coupled pde(s)
    // -> this triggers also the initialization of
    // those single PDEs which are directly coupled
    for (UInt i=0; i<numDirectCoupledPde_; i++) {
      Info->StartProgress("Initializing direct coupling");
      ptDirectCoupledPde_[i]->Init( sequenceStep );
      ptDirectCoupledPde_[i]->DefineAlgSys();
      Info->FinishProgress();
    }
    
    // Initialize coupledPDE
    if (ptIterCoupledPde_ != NULL) {
      Info->StartProgress("Initializing iterative coupling");
      ptIterCoupledPde_->InitCoupling( );
      Info->FinishProgress();
    }

    // Initialize algebraic system of each SinglePDE
    // Note: DefineAlgSys() triggers only the initialization 
    // of those SinglePDEs, which are not directly coupled
    for (UInt i = 0; i < numSinglePde_; i++ ) {
      it = isDirectCoupled_.find( ptSinglePde_[i] );
      if ( (*it).second == false) {
        ptSinglePde_[i]->DefineAlgSys();
      }
    }

  }

  // **************************
  //   Creation of PDEs
  // **************************
  void Domain::CreateSinglePDEs( UInt sequenceStep ) {


    // default grid
    Grid * defaultGrid = gridMap_["default"];

    StdVector<ParamNode*> pdeNodes;
    pdeNodes =  param->Get("sequenceStep", "index", GenStr(sequenceStep) )
      ->Get("pdeList")->GetChildren();

    ptSinglePde_.Resize( pdeNodes.GetSize() );
    ptSinglePde_.Init();
    numSinglePde_ = pdeNodes.GetSize();

    for (UInt i=0;i< pdeNodes.GetSize();i++) {

      std::string actPdeName = pdeNodes[i]->GetName();
      ParamNode * actPdeNode = pdeNodes[i];
      Info->StartProgress("Creating PDE '" + actPdeName + "'");
      
      if (actPdeName == "electrostatic") 
        ptSinglePde_[i]=new ElecPDE(defaultGrid, actPdeNode );

      else if (actPdeName == "mechanic")
        ptSinglePde_[i]=new MechPDE(defaultGrid, actPdeNode );

      else if (actPdeName == "acoustic") {
        
        std::string acouSubType = actPdeNode->Get("subType")->AsString();
        if (acouSubType == "flowNoise")
          ptSinglePde_[i]=new AcouFlowNoise(defaultGrid, actPdeNode);
        else if  (acouSubType == "combustionNoise")
          ptSinglePde_[i]=new AcouCombustionNoise(defaultGrid, actPdeNode);
	else
          ptSinglePde_[i]=new AcousticPDE(defaultGrid, actPdeNode );
      }
      
      else if (actPdeName == "acousticXYZ")
        ptSinglePde_[i]=new AcousticXYZPDE(defaultGrid, actPdeNode );
      
      else if (actPdeName == "acousticMixed")
        ptSinglePde_[i]=new AcousticMixedPDE(defaultGrid, actPdeNode );

      else if (actPdeName == "smooth")
        ptSinglePde_[i]=new SmoothPDE(defaultGrid, actPdeNode );

      else if (actPdeName == "magnetic") {
        //        if (dim_ == 2)
          ptSinglePde_[i]=new MagPDE(defaultGrid, actPdeNode );
        // else if (dim_ == 3)
        // ptSinglePde_[i]=new MagEdgePDE(ptgrid, actPdeNode); 
          //        else
          //          EXCEPTION( "Magnetic field calculation currently only possible in 2D!" );
      }
      
      else if (actPdeName == "mpcci")
        ptSinglePde_[i]=new MpcciPDE(defaultGrid, actPdeNode );

      else if (actPdeName == "heatConduction")
        ptSinglePde_[i]=new HeatCondPDE(defaultGrid, actPdeNode );

      else if (actPdeName == "stokesFluid")
        ptSinglePde_[i]=new StokesFluidPDE(defaultGrid, actPdeNode );

      else if (actPdeName == "bubble")
        ptSinglePde_[i]=new BubblePDE(defaultGrid, actPdeNode );

      else {
        EXCEPTION( actPdeName << " - this type of pdes is unknown" );
      }     

      // by default, not single pde is directly coupled
      isDirectCoupled_[ptSinglePde_[i]] = false;
      
      // Initialize current PDE
      // -> This step has now moved to method InitPDEs
      //ptSinglePde_[i]->Init();

      Info->FinishProgress();
    }

  } // end of InitPDE()


  void Domain::CreateIterCoupledPDE( UInt sequenceStep ) {
  
    std::string errMsg;

    // check if more than one PDEs are defined
    if (numIterCoupledStdPde_ <= 1) {
      ptIterCoupledPde_ = NULL;
      return;
    }
    
    Info->StartProgress("Creating coupling");
    
    // ================================
    //   Check for iterative coupling
    // ================================
    
    // check for presence of "couplingList" and "iterative" element
    ParamNode * couplingNode = 
      param->Get( "sequenceStep", "index", GenStr(sequenceStep) )
      ->Get("couplingList", false );
    if( !couplingNode) return;
    ParamNode * iterNode = couplingNode->Get("iterative", false );
    if (!iterNode) return;
    StdVector<ParamNode*> iterCplNodes = iterNode->GetChildren();

    
    // iterate over all pairwise iterative couplings
    StdVector<StdPDE*> iterCoupledPDEs;
    for( UInt i = 0; i < iterCplNodes.GetSize(); i++ ) {

      // HACK: Since the child nodes of the "iterative" node
      // do no only contain the two pdes, but also
      // the nonLinear node, we have to catch the latter case,
      // as we are only interested in the names of the pdes here
      if( iterCplNodes[i]->GetName() == "nonLinear") 
        continue;


      // fetch names of related pdes
      StdVector<ParamNode*> pdeNodes = iterCplNodes[i]->GetChildren();
      
      for(UInt iPde = 0; iPde < pdeNodes.GetSize(); iPde++ ) {
        std::string name = pdeNodes[iPde]->GetName();
        if( name!= "method") {
          SinglePDE * pde = GetSinglePDE( name );
          if( iterCoupledPDEs.Find( pde ) < 0 ) {
            iterCoupledPDEs.Push_back( pde );
          } 
        }
      } // pde
    } // couplings
    
    CoupledPDEDef * CouplingDef = new CoupledPDEDef(gridMap_["default"]);
    
    // create coupling objects 
    StdVector<StdPDE*> orderedPdes;
    CouplingDef->CreateCoupling(orderedPdes, couplings_, iterCoupledPDEs, iterNode );
    
    // Sort the different orderedPDEs into the singlePDEs
    StdVector<SinglePDE*> iterSinglePDEs;
    for (UInt i = 0; i< orderedPdes.GetSize(); i++) {
      iterSinglePDEs.Push_back( (SinglePDE*) orderedPdes[i]);
    }

    
    // Delete all of the singlePDEs, which are DirectCoupled and 
    // replace them by the direct-coupled ones. This is necessary,
    // since the iterCoupledPDE has to solve StdPDEs, whereas the
    // pairwise iterative couplings are defined only for 
    // SinglePDEs

    Integer index = 0; for (UInt i = 0; i<ptSinglePde_.GetSize(); i++ ) {
      if (isDirectCoupled_[ptSinglePde_[i]] == true ) {
        index = orderedPdes.Find(ptSinglePde_[i]);
        
        if ( index != -1 ) {
          orderedPdes[index] = ptDirectCoupledPde_[0];
          ptDirectCoupledPde_[0]->InitCoupling(NULL);
        }
      }
    }

    
    // create new iterative coupeld PDE
    ptIterCoupledPde_ = new IterCoupledPDE(orderedPdes, iterSinglePDEs,
                                           couplings_, iterNode);
    
    delete CouplingDef; 
    
    Info->FinishProgress();
  }
    

  // ***************************
  //   CreateDirectCoupledPDEs
  // ***************************
  void Domain::CreateDirectCoupledPDEs( UInt sequenceStep ) {


    numIterCoupledStdPde_ = numSinglePde_;

    // if only one pde exists, there is nothing to do
    if ( numSinglePde_ <= 1 ) {
      return;
    }
    
    // get "couplingList" node (must exist)
    ParamNode * couplingNode = 
      param->Get( "sequenceStep", "index", GenStr(sequenceStep) )
      ->Get("couplingList");
    ParamNode * directNode = couplingNode->Get("direct", false );
    if( !directNode ) return;

    // get nodes of pairwise direct couplings
    StdVector<ParamNode*> pairNodes = directNode->GetChildren();

    SinglePDE *pde1 = NULL;
    SinglePDE *pde2 =  NULL;
    BasePairCoupling *coupling = NULL;
      

    // HARD CODED: At the moment we allow only one direct coupled pde
    // with only on pairwise coupling
    StdVector<SinglePDE*> singlePdes;
    StdVector<BasePairCoupling*> DirectCouplingPairs;
    std::set<std::string> setSinglePDEs;
    std::string couplingName;
    
    for (UInt i=0; i<pairNodes.GetSize(); i++) {

      // get couplingName
      couplingName = pairNodes[i]->GetName();

      // *** PIEZO Coupling ***
      if ( couplingName == "piezoDirect" ) {

        pde1 = GetSinglePDE( "mechanic" );
        pde2 = GetSinglePDE( "electrostatic" );

        // in the case of piezo coupling, the electrotstatic
        // entries have to be multiplied by -1
        dynamic_cast<ElecPDE*>(pde2)->SetPiezoCoupling();
        
        coupling = new PiezoCoupling( pde1, pde2, pairNodes[i] );
        
      } 
      // *** ACOU-MECH Coupling ***
      else if ( couplingName == "acouMechDirect" ) {

        pde1 = GetSinglePDE( "mechanic" );
        pde2 = GetSinglePDE( "acoustic" );

        // in the case of acou-Mech coupling, the acoustic
        // entries have to be multiplied by -1
        dynamic_cast<AcousticPDE*>(pde2)->SetMechanicCoupling();

        coupling = new AcouMechCoupling( pde1, pde2, pairNodes[i] );
      }

// ------------------------------------------------------------------------
      // *** THERMO-MECH Coupling ***
      else if ( couplingName == "thermoMechDirect" ) {

         pde1 = GetSinglePDE( "mechanic" );
         pde2 = GetSinglePDE( "heatConduction" );


 		//turn the coupling on in heat-conduction pde
 		//in order to create the heat-conduction damp matrix
 		dynamic_cast<HeatCondPDE*>(pde2)->SetMechCoupling();


        coupling = new ThermoMechCoupling( pde1, pde2, pairNodes[i] );
      }
      // *** THERMO-ELECTRIC Coupling ***
      else if ( couplingName == "thermoElectricDirect" ) {

        pde1 = GetSinglePDE( "electrostatic" );
        pde2 = GetSinglePDE( "heatConduction" );
		
		//turn the coupling on in heat-conduction pde
		//in order to create the heat-conduction damp matrix
		dynamic_cast<HeatCondPDE*>(pde2)->SetElectroCoupling();
		
        // in the case of thermo coupling, the electrotstatic
        // entries have to be multiplied by -1
        dynamic_cast<ElecPDE*>(pde1)->SetThermoCoupling();


        coupling = new ThermoElectricCoupling( pde1, pde2, pairNodes[i] );
      }
//------------------------------------------------------------------------

      else {
        EXCEPTION( "The direct coupling '" << couplingName
                   << "' is not implemented!" << std::endl );
      }
      
      // set flag for direct coupling
      isDirectCoupled_[pde1] = true;
      isDirectCoupled_[pde2] = true;

      // add single PDEs and couplings into collections
      setSinglePDEs.insert( pde1->GetName() );
      setSinglePDEs.insert( pde2->GetName() );
      DirectCouplingPairs.Push_back(coupling);
    }

    // check if any pair coupling was found
    if (coupling == NULL)
      return;

    // Transform set of PDEs into a vector
    std::set<std::string>::iterator itSet;
    
    for (itSet = setSinglePDEs.begin();  itSet != setSinglePDEs.end(); 
         itSet++ ) {
      singlePdes.Push_back( GetSinglePDE(*itSet) );
    }
    
    ptDirectCoupledPde_.Push_back(new DirectCoupledPDE(gridMap_["default"], NULL) );
    ptDirectCoupledPde_[0]->SetSinglePDEs( singlePdes );
    ptDirectCoupledPde_[0]->SetCouplings( DirectCouplingPairs );

    // At the moment we allow only one direct coupled pde, so we set the
    // number of direct coupledPDEs to one;
    numDirectCoupledPde_ = 1;

    // now determine, how many SinglePDEs are coupling directly
    std::map<SinglePDE*,bool>::iterator it = isDirectCoupled_.begin();

    numIterCoupledStdPde_ = numDirectCoupledPde_;

    while (it != isDirectCoupled_.end() ) {
      if ( (*it).second == false )
        numIterCoupledStdPde_++;
      it++;
    }
    
    
  }

  void Domain::CreateCoordinateSystems() {

    // first create the "global" standard cartesian
    // coordinate system
    std::string defaultname = "default";
    coordSys_[defaultname] = new DefaultCoordSystem(gridMap_["default"]);

    // get nodes with coordinate systems
    ParamNode * coosyNode =
      param->Get("domain")->Get("coordSysList", false);
    if (!coosyNode) return;
    
    StdVector<ParamNode*> coordNodes = coosyNode->GetChildren();
    
    // iterate over all coordinate system nodes
    for( UInt i = 0; i < coordNodes.GetSize(); i++ ) {
      
      // fetch coosy name and type
      std::string type = coordNodes[i]->GetName();
      std::string name = coordNodes[i]->Get( "id" )->AsString();
      
      CoordSystem * actCoord = NULL;
      if( type == "polar") {
        actCoord =  new PolarCoordSystem( name, gridMap_["default"], 
                                          coordNodes[i] );
      } 
      else if ( type == "cylindric" ) {
        actCoord = new CylCoordSystem( name, gridMap_["default"], 
                                       coordNodes[i] );
      }
      else if ( type == "cartesian" ) {
        actCoord = new CartesianCoordSystem( name, gridMap_["default"], 
                                             coordNodes[i] );
      }
      else if ( type == "trivialCartesian" ) {
        actCoord = new TrivialCartesianCoordSystem( name, gridMap_["default"], 
                                                    coordNodes[i] );
      } else {
        EXCEPTION( "Coordinate system with type '" << type 
                   << "' not known!" );
      }
      coordSys_[name] = actCoord;
    }
  }

  

  // *************
  //   SetDriver
  // *************
  void Domain::SetDriver( BaseDriver * driver ) 
  {
    if( driver->GetAnalysisType() == BasePDE::MULTI_SEQUENCE ) {
      multiSequenceDriver_ = dynamic_cast<MultiSequenceDriver*>(driver);
      ptSingleDriver_ = multiSequenceDriver_->GetSingleDriver(); // NULL before Init()!!
    } else {
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
  void Domain::ResetPDEs() {


    // Delete single pde(s)
    for (UInt iPDE=0; iPDE<numSinglePde_; iPDE++) {
      delete ptSinglePde_[iPDE];
    }
    ptSinglePde_.Clear();

    // delete direct coupled pde(s)
    for (UInt iPDE=0; iPDE<numDirectCoupledPde_; iPDE++) {
      delete ptDirectCoupledPde_[iPDE]; 
    }
    ptDirectCoupledPde_.Clear();

    // delete iterative coupled pde
    if (ptIterCoupledPde_ != NULL) {
      delete ptIterCoupledPde_;
    }
    ptIterCoupledPde_ = NULL;

    //  // delete all couplings
    //     for (UInt iCoupl=0; iCoupl<couplings_.GetSize(); iCoupl++) {
    //       delete couplings_[iCoupl];
    //     }
    //     couplings_.Clear();

    // Also reset all state variables
    isDirectCoupled_.clear();
    numSinglePde_ = 0;
    numDirectCoupledPde_ = 0;
    numIterCoupledStdPde_ = 0;
  }

  void Domain::ReadErsatzMaterial(ParamNode* pn)
  {
    DesignElement::SetEnums();
    Optimization::SetEnums();
        
    // we read something like <loadErsatzMaterial region="piezo" file="piezo_density.xml" set="last"/>
    // Initialize our xerces dom parser to handle the external xml file
    Xerces* xerces = new Xerces(pn->Get("file")->AsString());
    // set the global ParamNode tree pointer
    ParamNode* xml = xerces->CreateParamNodeInstance();
    // release the xerces ressources, param is not affected
    delete xerces;
    // check this file
    if(xml->Count("set") == 0) 
      throw Exception("There are no design sets in the ersatz material file");

    // create a own structure for the ersatzMaterial
    std::string reg = pn->Get("region")->AsString();        
    if(!GetGrid()->HasRegion(reg)) 
       throw Exception("region given in loadErsatzMaterial is invalid");      
    RegionIdType regionId = GetGrid()->RegionNameToId(reg);
    
    // find the proper design set. This is either 'first', 'last' or the * in <set id="*"> ...
    ParamNode* set = NULL;
    std::string key = pn->Get("set")->AsString();
    if(key == "first") set = xml->GetList("set")[0];
    if(key == "last")  set = xml->GetList("set").Last(); 
    if(set == NULL)    set = xml->Get("set", "id", key);

    // the header is like
    // <header>
    //   <design name="density" initial="1.0"/>
    //   ...
    //   <transferFunction type="simp" application="mech" design="density" param="1.0"/>
    //   ..      
    // </header>
 
    // the design set consists of entries like
    // <element nr="401" type="density" design="0.886466" gradient="-7.56246e-09" filt_grad="-7.56246e-09"/>
    // only the combination nr and type is unique. E.g. in piezo we have types density and polarization
    StdVector<ParamNode*> des = xml->Get("header")->GetList("design");
    StdVector<ParamNode*> tfs = xml->Get("header")->GetList("transferFunction");      
    StdVector<ParamNode*> res(0); // empty

    // create the design space -> data has initial values!    
    ersatzMaterial = new DesignSpace(regionId, des, tfs, res);

    // read the set and replace the initial values
    StdVector<ParamNode*> elems = set->GetList("element");

    // check the the dimensions! the number of design variables comes from the regions and desings
    if(ersatzMaterial->data.GetSize() != elems.GetSize())
      EXCEPTION("ErsatzMaterialFile has " << elems.GetSize() << " entries, the model has "
                << ersatzMaterial->data.GetSize() << " entries");

    for(unsigned int e = 0; e < elems.GetSize(); e++)
    {
      unsigned int nr = elems[e]->Get("nr")->AsInt();
      DesignElement::Type dt = (DesignElement::Type) DesignElement::type.Parse(elems[e]->Get("type")->AsString());
      double val = elems[e]->Get("design")->AsDouble();
      
      // replace the value of the DesignElement
      DesignElement* de = ersatzMaterial->Find(nr, dt);
      de->SetDesign(val);
    }
  }

  bool Domain::GetErsatzMaterial(const Elem* elem, const BaseForm* form, double& result)
  {
    // is the stuff active at all? 
    if(ersatzMaterial == NULL) return false;

    // we cannot check for the region here, if form is a linear form (e.g.
    // pressure) but the design variable comes from elemens one dimension higher.
    int idx = ersatzMaterial->Find(elem, false);
    if(idx == -1) return false;
    
    // The desing space does the magic stuff. 
    // In the SIMP case we get density of element power param
    // all identified by the form and in piezo coupling case it
    // might even be the product of the transfer funcitons of
    // density and polarization
    result = ersatzMaterial->GetErsatzMaterialFactor(idx, form);
    return true;
  }

  DesignSpace* Domain::GetErsatzMaterial(bool throw_excpetion)
  {
    if(ersatzMaterial == NULL && throw_excpetion) 
      EXCEPTION("No ersatz material defined either via 'loadErsatzMaterial'"
                << " or an appropriate optimization");

    return ersatzMaterial;
  }


  // *************
  //   PrintGrid
  // *************
  void Domain::PrintGrid() {
    

    resultHandler_->Init( gridMap_["default"], true  );
    resultHandler_->Finalize();
  }
  
  void Domain::Dump()
  {
      for(UInt i = 0; i < ptSinglePde_.GetSize(); i++) {
         std::cout << "  single pde: " << ptSinglePde_[i]->GetName() << std::endl;
      }      

      for(UInt i = 0; i < couplings_.GetSize(); i++) {
         std::cout << "  coupled pde: " << couplings_[i]->GetPDE()->GetName() << std::endl;
      }
      
      for(UInt i = 0; i < ptDirectCoupledPde_.GetSize(); i++) {
         std::cout << "  direct coupled pde: " << ptDirectCoupledPde_[i]->GetName() << std::endl;
      }
      
  }

}
