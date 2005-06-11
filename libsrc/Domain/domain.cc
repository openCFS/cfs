#include "domain.hh"

#include <set>

#include "General/environment.hh"
#include "Domain/grid.hh"
#include "Domain/GridCFS/interface_gridcfs.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/writeresults.hh"

#include "PDE/pdes_header.hh"
#include "PDE/basePDE.hh"

// Coupling of single PDEs
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "CoupledPDE/itercoupledpde.hh"
#include "CoupledPDE/coupledpdedef.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "CoupledPDE/PiezoCoupling.hh"
#include "CoupledPDE/AcouMechCoupling.hh"


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
                  TimeFunc * aptTimeFunc) {

    ENTER_FCN( "Domain::Domain" );
  

    // initialize data
    numSinglePde_ = 0;
    numDirectCoupledPde_ = 0;
    numIterCoupledStdPde_ = 0;

    // assign pointers
    InFile_ = aptFileType; 
    OutFile_ = ptOut;
    ptTimeFunc_ = aptTimeFunc;
    ptIterCoupledPde_ = NULL;

    // read type of output results from conf-file
    std::string libmesh;
    params->Get( "meshLibrary", libmesh, "input" );

    UInt dim=InFile_->GetDim();

    std::string probGeo;

    // Check consistency of mesh and geometry type specified
    // in the param file
    params->Get( "type", probGeo, "geometry" );
    if ( !(probGeo == "3d"    && dim == 3 ||
           probGeo == "axi"   && dim == 2 ||
           probGeo == "plane" && dim == 2 ) ) {
      Info->Error( "Dimensions in parameter file and geometry file do not fit",
                   __FILE__, __LINE__ );
    }

    // initialize pointer to grid 
    if (dim==2) {

      if (libmesh =="cfsGrid") {
        ptgrid_=new GridInterfaceCFS<2>(InFile_);
      }

#ifdef ADAPTGRID
      else if (libmesh== "adaptGrid") {
        ptgrid_=new InterfaceAdaptGrid<2>(InFile_);
      }
#endif

      else {
        std::string errmsg = "Type of mesh_library should be one of ";
        errmsg += "'cfsgrid' or 'adaptgrid', but is '" + libmesh + "'";
        Info->Error( errmsg, __FILE__, __LINE__ );
      }
    }

    if (dim==3) {
      if (libmesh =="cfsGrid") {
        ptgrid_=new GridInterfaceCFS<3>(InFile_);
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

    Info->FinishProgress();
    
 
  }


  // **********************
  //   Default destructor
  // **********************
  Domain::~Domain() {
    ENTER_FCN( "Domain::~Domain" );

    delete ptgrid_;
    delete ptIterCoupledPde_;

    // When the StdVector ptpde_ is destroyed, only the pointers to the PDEs,
    // but not the PDEs themselves will be destroyed
    for ( UInt i = 0; i < ptSinglePde_.GetSize(); i++ ) {
      delete (ptSinglePde_[i]);
    }

    //already deleted in destructor of IterCoupledPDE!!!
    //     for ( UInt i = 0; i < couplings_.GetSize(); i++ ) {
    //       delete (couplings_[i]);
    //     }

  }


  // **********************
  //   Getter function
  // **********************
  StdPDE * Domain::GetStdPDE(const std::string pdeName)
  {
    ENTER_IFCN( "Domain::GetStdDPE" );
    Boolean pdeFound = FALSE;
    UInt i;
    std::string errMsg;

    // search the direct coupled pdes

    // *** IMPLEMENT ***

    // search the single pdes
    std::map<SinglePDE*,Boolean>::iterator it;
      
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
        
//         if ( (*it).second == FALSE ) {
//           pdeFound = TRUE;
//           break;
//         }
        pdeFound = TRUE;
        break;
      }
    }

    
    if (pdeFound == TRUE)
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
    Boolean pdeFound = FALSE;
    UInt i;
    std::string errMsg;

    for (i=0; i<ptSinglePde_.GetSize(); i++) {
      if (ptSinglePde_[i]->GetName() == pdeName) {
        pdeFound = TRUE;
        break;
      }
    }
    
    if (pdeFound == TRUE)
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

  // **************************
  //   Initialization of PDEs
  // **************************
  void Domain::InitPDEs( StdVector<std::string> & pdeNames,
                         UInt sequenceStep,
                         StdVector<std::string> tags ) {

    ENTER_FCN( "Domain::InitPDEs" );


    // create single pde(s)
    CreateSinglePDEs(pdeNames);    
    
    // create direct coupled pde(s)
    CreateDirectCoupledPDEs(tags);

    // intialize single pde(s)

    // Initialize those PDEs which are not
    // directly coupled
    std::map<SinglePDE*,Boolean>::iterator it;

    for (UInt i=0; i<numSinglePde_; i++) {
      it = isDirectCoupled_.find( ptSinglePde_[i] );
      if ( (*it).second == FALSE) {
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
    
    // Create iterative coupled pde
    CreateIterCoupledPDE(tags);

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
      ptSinglePde_[i]->DefineAlgSys();
    }

  }

  // **************************
  //   Creation of PDEs
  // **************************
  void Domain::CreateSinglePDEs(StdVector<std::string> & pdeNames) {

    ENTER_FCN( "Domain::CreatePDEs" );

    numSinglePde_ = pdeNames.GetSize();

    ptSinglePde_.Clear();
    ptSinglePde_.Resize(numSinglePde_);



    // Read dimension from mesh file and perform a consistency check
    UInt dim = InFile_->GetDim();

    for (UInt i=0;i< pdeNames.GetSize();i++) {
      Info->StartProgress("Creating PDE '" + pdeNames[i] + "'");

      if (pdeNames[i] == "electrostatic") 
        ptSinglePde_[i]=new ElecPDE(ptgrid_,ptTimeFunc_,OutFile_);

      else if (pdeNames[i] == "mechanic")
        ptSinglePde_[i]=new MechPDE(ptgrid_,ptTimeFunc_,OutFile_);

      else if (pdeNames[i] == "acoustic")
        {
          StdVector<std::string> acouSubType;
          params->GetList( "subType", acouSubType,"pdeList", "acoustic");
          if (acouSubType.GetSize())
            ptSinglePde_[i]=new AcouFlowNoise(ptgrid_,ptTimeFunc_,OutFile_);
          else
            ptSinglePde_[i]=new AcousticPDE(ptgrid_,ptTimeFunc_,OutFile_);
        }
      

      else if (pdeNames[i] == "smooth")
        ptSinglePde_[i]=new SmoothPDE(ptgrid_,ptTimeFunc_,OutFile_);

      else if (pdeNames[i] == "magnetic") 
        {
          if (dim == 2)
            ptSinglePde_[i]=new MagPDE(ptgrid_,ptTimeFunc_,OutFile_);
          else
            Error( "Magnetic field calculation currently only possible in 2D!",
                   __FILE__, __LINE__);
        }

      else if (pdeNames[i] == "piezo")
        ptSinglePde_[i]=new PiezoPDE(ptgrid_,ptTimeFunc_,OutFile_);

      else if (pdeNames[i] == "mpcci")
        ptSinglePde_[i]=new MpcciPDE(ptgrid_,ptTimeFunc_,OutFile_);

      //       else if (pdeNames[i] == "acouflownoise")
      //              ptSinglePde_[i]=new AcouFlowNoise(ptgrid_,ptTimeFunc_,OutFile_);

      //      else if (pdeNames[i] == "smoothlaplace") 
      //        ptSinglePde_[i]=new SmoothLaPlacePDE(ptgrid_,ptTimeFunc_,OutFile_); 


      //      else if (pdeNames[i] == "magnetic") 
      //        if (dim == 3
      //        ptSinglePde_[i]=new MagEdgePDE(ptgrid_,ptTimeFunc_,OutFile_); 
      else
        {
          std::string msg=pdeNames[i]+" - this type of pdes is unknown";
          Error(msg.c_str(),__FILE__,__LINE__);
        }     
      
      // by default, not single pde is directly coupled
      isDirectCoupled_[ptSinglePde_[i]] = FALSE;
      
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
      for (UInt i=1; i<sequenceTags.GetSize(); i++)
        if (sequenceTags[i] != firstTag) {
          errMsg = "CreateIterCoupledPDE: The tags in the <multiSequence> section ";
          errMsg += "are not all the same in each step.\n Coupling is only ";
          errMsg += "possible if in each step all the <refTag> are the same!";
          Error(errMsg.c_str(), __FILE__, __LINE__ );
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
    Integer index = 0;
    for (UInt i = 0; i<ptSinglePde_.GetSize(); i++ ) {
      if (isDirectCoupled_[ptSinglePde_[i]] == TRUE ) {
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
    

  void Domain::CreateDirectCoupledPDEs(StdVector<std::string> & sequenceTags) {
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
    StdVector<BasePairCoupling*> couplings;
    std::set<SinglePDE*> setSinglePDEs;

    for (UInt i=0; i<couplingNames.GetSize(); i++) {
      //std::cerr << "Coupling " << i+1 << " = " << couplingNames[i] << std::endl;

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
    isDirectCoupled_[pde1] = TRUE;
    isDirectCoupled_[pde2] = TRUE;

    // add single PDEs and couplings into collections
    setSinglePDEs.insert(pde1);
    setSinglePDEs.insert(pde2);
    couplings.Push_back(coupling);

    }


    // check if any pair coupling was found
    if (coupling == NULL)
      return;

    // Transform set of PDEs into a vector
    std::set<SinglePDE*>::iterator itSet;
    
    for (itSet = setSinglePDEs.begin();  itSet != setSinglePDEs.end(); 
         itSet++ ) {
      singlePdes.Push_back(*itSet);
    }
    

   
    
   

    ptDirectCoupledPde_.Push_back(new DirectCoupledPDE(ptgrid_, OutFile_,
                                                       ptTimeFunc_));
    ptDirectCoupledPde_[0]->SetSinglePDEs( singlePdes );
    ptDirectCoupledPde_[0]->SetCouplings( couplings );


    // At the moment we allow only one direct coupled pde, so we set the number
    // of direct coupledPDEs to one;
    numDirectCoupledPde_ = 1;

    // now determine, how many SinglePDEs are coupling directly
    std::map<SinglePDE*,Boolean>::iterator it = isDirectCoupled_.begin();

    numIterCoupledStdPde_ = numDirectCoupledPde_;

    while (it != isDirectCoupled_.end() ) {
      if ( (*it).second == FALSE )
        numIterCoupledStdPde_++;
      it++;
    }
    
    
  }

  void Domain::ResetPDEs()
  {
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


  void Domain::PrintGrid( ) {
    ENTER_FCN( "Domain::PrintGrid" );
    if ( ptgrid_ == NULL )
      Error("Domain: ptgrid == NULL!");
    OutFile_->Init(ptgrid_);
    OutFile_->WriteGrid( );
  }


  void Domain::Update( ) {
    ENTER_FCN( "Domain::Update" );
 
    // Init AlgSystem
    UpdateAlgSys();
  }


  void Domain::UpdateAlgSys( ) {
    ENTER_FCN( "Domain::UpdateAlgSys"  );

    //set the algebraic systems and read material data
    for (UInt i=0;i< numSinglePde_;i++)
      {
        ptSinglePde_[i]->DeleteAlgSys();
        
        // the 'Reset' function was never implemented
        // for a PDE, so what sould the function call
        // trigger?
        //ptSinglePde_[i]->Reset();
        ptSinglePde_[i]->DefineAlgSys(); 
      }
  }

}
