#include "domain.hh"

#include <set>

#include "General/environment.hh"
#include "Domain/grid.hh"
#include "Domain/GridCFS/interface_gridcfs.hh"
#include "Domain/GridStruct/interface_gridstruct.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/MaterialHandler.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/writeresults.hh"
#include "Utils/coordSystem.hh"
#include "Utils/cylCoordSys.hh"
#include "Utils/polCoordSys.hh"
#include "Utils/defaultCoordSys.hh"
#include "DataInOut/Scripting/cfsmessenger.hh"

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

#ifdef NETGEN
#include "interface_netgen.hh"
#endif

#ifdef GRIDLIB
#include "interface_gridlib.hh"
#endif

#ifdef ADAPTGRID
#include "domain/AdaptGrid/interface_adgrid.hh"
#endif



namespace CoupledField {

  // **************
  //   Construtor
  // **************
  Domain:: Domain(FileType * const aptFileType, WriteResults * ptOut,
                  TimeFunc * aptTimeFunc, MaterialHandler * ptMat) {

    ENTER_FCN( "Domain::Domain" );

    Info->PrintF("","==================\n");
    Info->PrintF("","   DOMAIN SETUP   \n");
    Info->PrintF("","==================\n\n");

    // initialize data
    numSinglePde_ = 0;
    numDirectCoupledPde_ = 0;
    numIterCoupledStdPde_ = 0;

    // assign pointers
    InFile_ = aptFileType; 
    OutFile_ = ptOut;
    ptMatHandler_ = ptMat;
    ptTimeFunc_ = aptTimeFunc;
    ptIterCoupledPde_ = NULL;
    ptSingleDriver_ = NULL;

    // read type of mesh library
    std::string libmesh;
    params->Get( "meshLibrary", libmesh, "input" );

    std::string probGeo;
    params->Get( "type", probGeo, "geometry" );

    if ( libmesh == "structGrid" ) {
      if ( probGeo == "3d") {
        dim_ = 3;
      } 
      else if ( probGeo == "axi" || probGeo == "plane" ) {
        dim_ = 2;
      }
      else {
        Info->Error( "Wrong Dimension parameter in xml-File", 
                     __FILE__, __LINE__ );
      }
    }
    else {
      dim_ = InFile_->GetDim();

      // Check consistency of mesh and geometry type specified
      // in the param file
      params->Get( "type", probGeo, "geometry" );
      if ( !(probGeo == "3d"    && dim_ == 3 ||
             probGeo == "axi"   && dim_ == 2 ||
             probGeo == "plane" && dim_ == 2  ) ) {
        (*error) << "Dimensions in parameter file and geometry file "
                 << "do not match!";
        Error( __FILE__, __LINE__ );
      }
    }

    Info->PrintF( "", "Type of grid: %s\nDimension: %i\n\n",
                  libmesh.c_str(), dim_ );
    SETPROFILE("Before Grid-Creation");

    // initialize pointer to grid 
    if (dim_==2) {

      if (libmesh =="cfsGrid") {
        ptgrid_=new GridInterfaceCFS<2>(InFile_);
      }
      else if (libmesh == "structGrid" ) {
        ptgrid_=new GridInterfaceStruct<2>(dim_);
      }

#ifdef ADAPTGRID
      else if (libmesh== "adaptGrid") {
        ptgrid_=new InterfaceAdaptGrid<2>(InFile_);
      }
#endif

      else {
        (*error) << "Type of mesh_library should be one of "
                 << "'cfsgrid' or 'adaptgrid', but is '" << libmesh << "'";
        Error( __FILE__, __LINE__ );
      }
    }

    if (dim_==3) {
      if (libmesh =="cfsGrid") {
        ptgrid_=new GridInterfaceCFS<3>(InFile_);
      }
      else if (libmesh == "structGrid" ) {
        ptgrid_=new GridInterfaceStruct<3>(dim_);
      }


#ifdef ADAPTGRID
      else if (libmesh== "adaptGrid") {
        ptgrid_=new InterfaceAdaptGrid<3>(InFile_);
      }
#endif
      else {
        Error("Unknown type of mesh_library in conf-file",__FILE__,__LINE__);
      }
    }
    
    Info->StartProgress("Reading in the mesh");

    //read in the mesh information
    ptgrid_->Read();
    SETPROFILE("After Grid Creation");

    // Read in coordinate systems
    CreateCoordinateSystems();

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
          (*error) << "It was impossible to determine, if the PDE "
                   << "with the name '" << pdeName << "' is directly "
                   << "coupled or not";
          Error( __FILE__, __LINE__ );
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
      errMsg = "Domain:GetPDE: PDE with name '";
      errMsg += pdeName;
      errMsg += "' was not found/created!.";
      Error(errMsg.c_str(), __FILE__, __LINE__);
      return NULL;
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
      errMsg = "Domain:GetSinglePDE: PDE with name '";
      errMsg += pdeName;
      errMsg += "' was not found/created!.";
      Error(errMsg.c_str(), __FILE__, __LINE__);
      return NULL;
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
      (*error) << "Domain::GetCoordSystem: A coordinate system with name '"
               << name << "' is not registered!";
      Error( __FILE__, __LINE__ );
    }

    return (*it).second;
  
  }


  // **************************
  //   Initialization of PDEs
  // **************************
  void Domain::CreatePDEs( StdVector<std::string> & pdeNames,
                         UInt sequenceStep,
                         StdVector<std::string> tags ) {

    ENTER_FCN( "Domain::CreatePDEs" );

    // create single pde(s)
    CreateSinglePDEs(pdeNames);    
    
    // create direct coupled pde(s)
    CreateDirectCoupledPDEs(tags);

    // Create iterative coupled pde
    CreateIterCoupledPDE(tags);
  }

  void Domain::InitPDEs(  UInt sequenceStep,
                          StdVector<std::string> tags ) {

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
        ptSinglePde_[i]->Init(sequenceStep,tags[i]);
      }
    }

    // initialize direct coupled pde(s)
    // -> this triggers also the initialization of
    // those single PDEs which are directly coupled
    for (UInt i=0; i<numDirectCoupledPde_; i++) {
      Info->StartProgress("Initializing direct coupling");
      ptDirectCoupledPde_[i]->Init(sequenceStep,tags[i]);
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
  void Domain::CreateSinglePDEs(StdVector<std::string> & pdeNames) {

    ENTER_FCN( "Domain::CreatePDEs" );

    numSinglePde_ = pdeNames.GetSize();

    ptSinglePde_.Resize(numSinglePde_);
    ptSinglePde_.Init();

    for (UInt i=0;i< pdeNames.GetSize();i++) {
      Info->StartProgress("Creating PDE '" + pdeNames[i] + "'");

      if (pdeNames[i] == "electrostatic") 
        ptSinglePde_[i]=new ElecPDE(ptgrid_,ptTimeFunc_,OutFile_);

      else if (pdeNames[i] == "mechanic")
        ptSinglePde_[i]=new MechPDE(ptgrid_,ptTimeFunc_,OutFile_);

      else if (pdeNames[i] == "acoustic") {
        StdVector<std::string> acouSubType;
        params->GetList( "subType", acouSubType,"pdeList", "acoustic");
        if (acouSubType[0] == "flowNoise")
          ptSinglePde_[i]=new AcouFlowNoise(ptgrid_,ptTimeFunc_,OutFile_);
        else if  (acouSubType[0] == "combustionNoise")
          ptSinglePde_[i]=new AcouCombustionNoise(ptgrid_,ptTimeFunc_,OutFile_);
	else
          ptSinglePde_[i]=new AcousticPDE(ptgrid_,ptTimeFunc_,OutFile_);
      }
      
      else if (pdeNames[i] == "acousticXYZ")
        ptSinglePde_[i]=new AcousticXYZPDE(ptgrid_,ptTimeFunc_,OutFile_);

      else if (pdeNames[i] == "nrbc")
        {
          StdVector<SolutionType> type(1);
          type[0] = NRBC_PHI;
          ptSinglePde_[i]=new nrbcPDE(ptgrid_,ptTimeFunc_,OutFile_, "nrbc", type);
        }
      
      else if (pdeNames[i] == "smooth")
        ptSinglePde_[i]=new SmoothPDE(ptgrid_,ptTimeFunc_,OutFile_);

      else if (pdeNames[i] == "magnetic") {
        if (dim_ == 2)
          ptSinglePde_[i]=new MagPDE(ptgrid_,ptTimeFunc_,OutFile_);
        // else if (dim_ == 3)
        // ptSinglePde_[i]=new MagEdgePDE(ptgrid_,ptTimeFunc_,OutFile_); 
        else
          Error( "Magnetic field calculation currently only possible in 2D!",
                 __FILE__, __LINE__);
      }
      
      else if (pdeNames[i] == "mpcci")
        ptSinglePde_[i]=new MpcciPDE(ptgrid_,ptTimeFunc_,OutFile_);

      else if (pdeNames[i] == "heatConduction")
        ptSinglePde_[i]=new HeatCondPDE(ptgrid_,ptTimeFunc_,OutFile_);

      else if (pdeNames[i] == "stokesFluid")
        ptSinglePde_[i]=new StokesFluidPDE(ptgrid_,ptTimeFunc_,OutFile_);

      else if (pdeNames[i] == "bubble")
        ptSinglePde_[i]=new BubblePDE(ptgrid_,ptTimeFunc_,OutFile_);

      else {
        std::string msg=pdeNames[i]+" - this type of pdes is unknown";
        Error(msg.c_str(),__FILE__,__LINE__);
      }     

      // by default, not single pde is directly coupled
      isDirectCoupled_[ptSinglePde_[i]] = false;
      
      // Initialize current PDE
      // -> This step has now moved to method InitPDEs
      //ptSinglePde_[i]->Init();

      Info->FinishProgress();
    }

  } // end of InitPDE()


  void Domain::CreateIterCoupledPDE(StdVector<std::string> & sequenceTags) {
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
    
    StdVector<StdPDE*> iterCoupledPDEs;
    StdVector<SinglePDE*> iterSinglePDEs;
    StdVector<std::string> iterCoupledPDENames;
    StdVector<std::string> methods;
    
    // Check if all sequenceTags are the same.
    // Currently the tag for the current coupling in a multisequence
    // analysis is determined by simply taking the first tag of
    // all tags for the PDEs, which have to be the same.
    std::string firstTag;
    if (sequenceTags.GetSize() > 1) {
      firstTag = sequenceTags[0];
      for ( UInt i = 1; i < sequenceTags.GetSize(); i++ ) {
        if ( sequenceTags[i] != firstTag ) {
          (*error) << "CreateIterCoupledPDE: The tags in the "
                   << "<multiSequence> section "
                   << "are not all the same in each step.\n "
                   << "Coupling is only possible if in each step all the "
                   << "<refTag> tags are the same!";
          Error( __FILE__, __LINE__ );
        }
      }
    }

    params->GetIterCoupledPDEList(iterCoupledPDENames, sequenceTags[0]);
    

    // we have all the names of the PDEs which couple iteratively.
    // Now we have to get the according pointers to the PDEs
    for (UInt i=0; i<iterCoupledPDENames.GetSize(); i++) {
      for (UInt j=0; j<ptSinglePde_.GetSize(); j++) {
        if (iterCoupledPDENames[i] == ptSinglePde_[j]->GetName()) {
          iterCoupledPDEs.Push_back(ptSinglePde_[j]);
        }
      }
    }
  

    // Determine method of coupling (currently only RHS is
    // given)
    params->GetList( "method", methods);
    for (UInt i=0; i<methods.GetSize(); i++)
      if (methods[i] != "RHS")
        {
          errMsg  = "Domain::InitCoupledPDE: Methode '";
          errMsg += methods[i];
          errMsg += "' not implemented for iterative coupling";
          Error(errMsg.c_str(), __FILE__, __LINE__);
        } 
    
    
    CoupledPDEDef * CouplingDef = new CoupledPDEDef(ptgrid_);
    
    // create coupling objects 
    StdVector<StdPDE*> orderedPdes;
    CouplingDef->CreateCoupling(orderedPdes, couplings_, iterCoupledPDEs);
    
    // Sort the different orderedPDEs into the singlePDEs
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
                                           couplings_, sequenceTags[0]);
    
    delete CouplingDef; 
    
    Info->FinishProgress();
  }
    

  // ***************************
  //   CreateDirectCoupledPDEs
  // ***************************
  void Domain::CreateDirectCoupledPDEs(StdVector<std::string> &sequenceTags) {

    ENTER_FCN( "Domain::CreateDirectCoupledPDEs" );

    numIterCoupledStdPde_ = numSinglePde_;

    // if only one pde exists, there is nothing to do
    if ( numSinglePde_ <= 1 ) {
      return;
    }

    // otherwise check if some direct coupling exists
    StdVector<std::string> couplingNames;
    params->GetDirectCouplingList(couplingNames, sequenceTags[0]);


    SinglePDE *pde1 = NULL;
    SinglePDE *pde2 =  NULL;
    BasePairCoupling *coupling = NULL;
      

    // HARD CODED: At the moment we allow only one direct coupled pde
    // with only on pairwise coupling
    StdVector<SinglePDE*> singlePdes;
    StdVector<BasePairCoupling*> DirectCouplingPairs;
    std::set<std::string> setSinglePDEs;

    for (UInt i=0; i<couplingNames.GetSize(); i++) {

      // *** PIEZO Coupling ***
      if ( couplingNames[i] == "piezoDirect" ) {

        pde1 = GetSinglePDE( "mechanic" );
        pde2 = GetSinglePDE( "electrostatic" );

        // in the case of piezo coupling, the electrotstatic
        // entries have to be multiplied by -1
        dynamic_cast<ElecPDE*>(pde2)->SetPiezoCoupling();
        
        coupling = new PiezoCoupling( pde1, pde2 );
        
      } 
      // *** ACOU-MECH Coupling ***
      else if ( couplingNames[i] == "acouMechDirect" ) {

        pde1 = GetSinglePDE( "mechanic" );
        pde2 = GetSinglePDE( "acoustic" );

        // in the case of acou-Mech coupling, the acoustic
        // entries have to be multiplied by -1
        dynamic_cast<AcousticPDE*>(pde2)->SetMechanicCoupling();

        coupling = new AcouMechCoupling( pde1, pde2 );
      }
      else {
        (*error) << "The direct coupling '" << couplingNames[i]
                 << "' is not implemented!" << std::endl;
        Error( __FILE__, __LINE__ );
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
    
    ptDirectCoupledPde_.Push_back(new DirectCoupledPDE(ptgrid_, OutFile_,
                                                       ptTimeFunc_));
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


    // then read in the additional coordinate systems
    StdVector<std::string> keyVec, attrVec, valVec, names;
    attrVec = "", "", "";
    valVec = "", "", "";
    
    // a) cylindrical coordinate systems
    keyVec = "domain", "coordSys", "cylindric", "name";
    params->GetList(keyVec, attrVec, valVec, names);

    for (UInt i = 0; i < names.GetSize(); i++ ) {
      CoordSystem * actCoord = new CylCoordSystem(names[i],ptgrid_);
      coordSys_[names[i]] = actCoord;
    }

    // b) polar coordinate systems
    names.Clear();
    keyVec = "domain", "coordSys", "polar", "name";
    params->GetList(keyVec, attrVec, valVec, names);

    for (UInt i = 0; i < names.GetSize(); i++ ) {
      CoordSystem * actCoord = 
        new PolarCoordSystem(names[i],ptgrid_);
      coordSys_[names[i]] = actCoord;
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

    if ( ptgrid_ == NULL ) {
      Error( "Domain: ptgrid == NULL!", __FILE__, __LINE__ );
    }
    OutFile_->Init( ptgrid_ );
    OutFile_->WriteGrid();
  }

}
