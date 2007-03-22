// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "domain.hh"

#include <set>
#include <memory>

#include "General/environment.hh"
#include "Domain/grid.hh"
#include "Domain/GridCFS/grid_cfs.hh"
#include "DataInOut/simInput.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/MaterialHandler.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/simInput.hh"
#include "General/exception.hh"
#include "Utils/coordSystem.hh"
#include "Utils/cylCoordSys.hh"
#include "Utils/polCoordSys.hh"
#include "Utils/defaultCoordSys.hh"
#include "DataInOut/Scripting/cfsmessenger.hh"
#include "DataInOut/resultHandler.hh"

#include "PDE/pdes_header.hh"
#include "PDE/basePDE.hh"

// Coupling of single PDEs
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "CoupledPDE/itercoupledpde.hh"
#include "CoupledPDE/coupledpdedef.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "CoupledPDE/PiezoCoupling.hh"
#include "CoupledPDE/AcouMechCoupling.hh"

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

  Domain::Domain(SimInput * const aptSimInput, ResultHandler * handler,
                 MaterialHandler * ptMat ) {
    
    
    ENTER_FCN( "Domain::Domain" );
    
    Info->PrintF("","==================\n");
    Info->PrintF("","   DOMAIN SETUP   \n");
    Info->PrintF("","==================\n\n");
    
    // initialize data
    numSinglePde_ = 0;
    numDirectCoupledPde_ = 0;
    numIterCoupledStdPde_ = 0;
    
    // assign pointers
    simInput_ = aptSimInput; 
    ptMatHandler_ = ptMat;
    resultHandler_  = handler;
    ptIterCoupledPde_ = NULL;
    ptSingleDriver_ = NULL;
  }
  
  void Domain::CreateGrid() {
    ENTER_FCN( "Domain::CreateGrid" );
    
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
    
    SETPROFILE("Before Grid-Creation");
    
    // initialize pointer to grid 
    if (libmesh =="cfsGrid") {
      ptgrid_=new GridCFS();
    }
    
#ifdef ADAPTGRID
    else if (libmesh== "adaptGrid") {
      ptgrid_=new InterfaceAdaptGrid<2>(InFile_);
    }
#endif
    
    else {
      EXCEPTION( "Type of mesh_library should be one of "
                 << "'cfsgrid' or 'adaptgrid', but is '" << libmesh << "'" );
    }
    
    Info->StartProgress("Reading in the mesh");
    
    //     try 
    //     {
    //         EXCEPTION("Test WARNING in Domain! " << 2);
    //     }
    //     catch (Exception& ex)
    //     {
    //         RETHROW_EXCEPTION(ex, "Ich habe eine Exception gefangen. " << 3);
    //     }
    
    simInput_->InitModule(ptgrid_);
    simInput_->ReadMesh();
    
    // Read in coordinate systems
    // Note: It is important that is 
    CreateCoordinateSystems();
    ptgrid_->FinishInit();
    
    SETPROFILE("After Grid Creation");

    

    Info->FinishProgress();

    Info->StartProgress("Initializing non-matching interfaces");

    // Call the nonmatching grid intersection calculation
    ptgrid_->InitNonmatchingInterfaces();

    // Initialize resultHandler
    resultHandler_->Init( ptgrid_ );

    Info->PrintF("","\n=========================\n");
    Info->PrintF("","   END OF DOMAIN SETUP   \n");
    Info->PrintF("","=========================\n\n");
    Info->FinishProgress();
  }


  // **************
  //   Destructor
  // **************
  Domain::~Domain() {

    ENTER_FCN( "Domain::~Domain" );

    delete ptgrid_;
    ptgrid_ = NULL;
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

  }


  // **********************
  //   Getter function
  // **********************
  StdPDE * Domain::GetStdPDE(const std::string pdeName)
  {
    ENTER_IFCN( "Domain::GetStdDPE" );
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
    ENTER_IFCN( "Domain::GetSingleDPE" );
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
    ENTER_IFCN( "Domain::GetDPE" );
    
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

  CoordSystem * Domain::GetCoordSystem( const std::string & name ) {
    ENTER_FCN( "Domain::GetCoordSystem" );
    
    
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

    ENTER_FCN( "Domain::CreatePDEs" );

    // create single pde(s)
    CreateSinglePDEs( sequenceStep );    
    
    // create direct coupled pde(s)
    CreateDirectCoupledPDEs( sequenceStep );

    // Create iterative coupled pde
    CreateIterCoupledPDE( sequenceStep );
  }

  void Domain::InitPDEs( UInt sequenceStep ) {

    ENTER_FCN( "Domain::InitPDEs" );

#ifdef USE_SCRIPTING

    // Call initialization procedure
    StdVector<std::string> context;
    context.Push_back( commandLine->GetSimName() );
    messenger->TriggerEvent(CFSMessenger::CFS_Init, context);
#endif

    // intialize single pde(s)

    // Initialize those PDEs which are not
    // directly coupled
    std::map<SinglePDE*,bool>::iterator it;

    for (UInt i=0; i<numSinglePde_; i++) {
      it = isDirectCoupled_.find( ptSinglePde_[i] );
      if ( (*it).second == false) {
        //std::cerr << "Domain: Init() of " 
        //<< ptSinglePde_[i]->GetName() << std::endl;
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

    ENTER_FCN( "Domain::CreatePDEs" );

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
        ptSinglePde_[i]=new ElecPDE(ptgrid_, actPdeNode );

      else if (actPdeName == "mechanic")
        ptSinglePde_[i]=new MechPDE(ptgrid_, actPdeNode );

      else if (actPdeName == "acoustic") {
        
        std::string acouSubType = actPdeNode->Get("subType")->AsString();
        if (acouSubType == "flowNoise")
          ptSinglePde_[i]=new AcouFlowNoise(ptgrid_, actPdeNode);
        else if  (acouSubType == "combustionNoise")
          ptSinglePde_[i]=new AcouCombustionNoise(ptgrid_, actPdeNode);
	else
          ptSinglePde_[i]=new AcousticPDE(ptgrid_, actPdeNode );
      }
      
      else if (actPdeName == "acousticXYZ")
        ptSinglePde_[i]=new AcousticXYZPDE(ptgrid_, actPdeNode );

      else if (actPdeName == "smooth")
        ptSinglePde_[i]=new SmoothPDE(ptgrid_, actPdeNode );

      else if (actPdeName == "magnetic") {
        if (dim_ == 2)
          ptSinglePde_[i]=new MagPDE(ptgrid_, actPdeNode );
        // else if (dim_ == 3)
        // ptSinglePde_[i]=new MagEdgePDE(ptgrid, actPdeNode); 
        else
          EXCEPTION( "Magnetic field calculation currently only possible in 2D!" );
      }
      
      else if (actPdeName == "mpcci")
        ptSinglePde_[i]=new MpcciPDE(ptgrid_, actPdeNode );

      else if (actPdeName == "heatConduction")
        ptSinglePde_[i]=new HeatCondPDE(ptgrid_, actPdeNode );

      else if (actPdeName == "stokesFluid")
        ptSinglePde_[i]=new StokesFluidPDE(ptgrid_, actPdeNode );

      else if (actPdeName == "bubble")
        ptSinglePde_[i]=new BubblePDE(ptgrid_, actPdeNode );

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
    ENTER_FCN( "Domain::CreateIterCoupledPDE" );
  
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
    
    CoupledPDEDef * CouplingDef = new CoupledPDEDef(ptgrid_);
    
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

    ENTER_FCN( "Domain::CreateDirectCoupledPDEs" );

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
    
    ptDirectCoupledPde_.Push_back(new DirectCoupledPDE(ptgrid_, NULL) );
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
    ENTER_FCN( "Domain::CreateCoordinateSystems");

    // first create the "global" standard cartesian
    // coordinate system
    std::string defaultname = "default";
    coordSys_[defaultname] = new DefaultCoordSystem(ptgrid_);

    // get nodes with coordinate systems
    ParamNode * coosyNode =
      param->Get("domain")->Get("coordSys", false);
    if (!coosyNode) return;
    
    StdVector<ParamNode*> coordNodes = coosyNode->GetChildren();
    
    // iterate over all coordinate system nodes
    for( UInt i = 0; i < coordNodes.GetSize(); i++ ) {
      
      // fetch coosy name and type
      std::string type = coordNodes[i]->GetName();
      std::string name = coordNodes[i]->Get( "name" )->AsString();
      
      CoordSystem * actCoord = NULL;
      if( type == "polar") {
        actCoord =  new PolarCoordSystem( name, ptgrid_, coordNodes[i] );
      } 
      else if ( type == "cylindric" ) {
        actCoord = new CylCoordSystem( name, ptgrid_, coordNodes[i] );
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
  void Domain::SetDriver( BaseDriver * driver ) {

    if( driver->GetAnalysisType() == MULTI_SEQUENCE ) {
      ptSingleDriver_ = dynamic_cast<MultiSequenceDriver*>(driver)
        ->GetSingleDriver();
    } else {
      ptSingleDriver_ = dynamic_cast<SingleDriver*>(driver);
    }
  }
  
  // *************
  //   ResetPDEs
  // *************
  void Domain::ResetPDEs() {

    ENTER_FCN( "Domain::ResetDEs" );

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


  // *************
  //   PrintGrid
  // *************
  void Domain::PrintGrid() {
    
    ENTER_FCN( "Domain::PrintGrid" );
    
    resultHandler_->WriteGrid();
  }
  
}
