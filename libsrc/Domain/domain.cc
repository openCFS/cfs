#include <iostream>
#include <fstream>
#include <string>

#include "domain.hh"
#include "Domain/grid.hh"
#include "Domain/bcs.hh"
#include "Domain/GridCFS/interface_gridcfs.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "CoupledPDE/coupledpdedef.hh"
#include "CoupledPDE/pdecoupling.hh"

#ifdef NETGEN
#include "interface_netgen.hh"
#endif

#ifdef GRIDLIB
#include "interface_gridlib.hh"
#endif

#ifdef ADAPTGRID
#include "domain/AdaptGrid/interface_adgrid.hh"
#endif

#include "PDE/pdes_header.hh"
#include "CoupledPDE/coupled_pdes_header.hh"

#include "Utils/vector.hh"
#include "PDE/basePDE.hh"

namespace CoupledField {

  // **************
  //   Construtor
  // **************
  Domain:: Domain(FileType * const aptFileType, WriteResults * ptOut,
		  TimeFunc * aptTimeFunc) {

    ENTER_FCN( "Domain::Domain" );
  
    numlevel_ = 0;

    InFile_     = aptFileType; 
    OutFile_    = ptOut;
    ptTimeFunc_ = aptTimeFunc;
    ptcoupledpde_ = NULL;

    // read type of output results from conf-file
    std::string libmesh;
    params->Get( "meshLibrary", libmesh, "input" );

    Integer dim=InFile_->ReadDim();

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
    
 
    if (PrintGridOnly) {
      PrintGrid(0);
      exit(0);
    }

    // allocate an object with an information about boundary condition
    ptBCs_=new BCs(InFile_);

    //read restraints information
    ptBCs_->ReadBCs();

  }


  // **********************
  //   Default destructor
  // **********************
  Domain::~Domain() {
    ENTER_FCN( "Domain::~Domain" );
    if (ptgrid_) delete ptgrid_;
    if (ptBCs_) delete ptBCs_;
  }

  // **********************
  //   Getter function
  // **********************
  BasePDE * Domain::GetPDE(const std::string pdeName)
  {
    ENTER_IFCN( "Domain::GetDPE" );

    Boolean pdeFound = FALSE;
    Integer i;
    std::string errMsg;

    for (i=0; i<numpde_; i++)
      if (ptpde_[i]->GetName() == pdeName) {
	pdeFound = TRUE;
	break;
      }
    if (pdeFound == TRUE)
      return ptpde_[i];
    else {
      errMsg = "Domain:GetPDE: PDE with name '";
      errMsg += pdeName;
      errMsg += "' was not found/created!.";
      Error(errMsg.c_str(), __FILE__, __LINE__);
      return NULL;
    }
  }

  // **************************
  //   Initialization of PDEs
  // **************************
  void Domain::InitPDEs( StdVector<std::string> & pdeNames,
			 Integer sequenceStep,
			 StdVector<std::string> tags ) {

    ENTER_FCN( "Domain::InitPDEs" );

    CreatePDEs(pdeNames);    

    for (Integer i=0; i<numpde_; i++) {
      ptpde_[i]->Init(sequenceStep,tags[i]);
     }

    CreateCoupledPDE(tags);

    // Initialize coupledPDE
    if (numpde_ > 1) {
      Info->StartProgress("Initializing Coupling");
      ptcoupledpde_->InitCoupling(numlevel_);
      Info->FinishProgress();
    }

    // Initialize algebraic system of each PDE
    for (int i = 0; i < numpde_; i++ ) {
      ptpde_[i]->SetAlgSys();
    }

  }

  // **************************
  //   Creation of PDEs
  // **************************
  void Domain::CreatePDEs(StdVector<std::string> & pdeNames) {

    ENTER_FCN( "Domain::CreatePDEs" );

    numpde_ = pdeNames.GetSize();

    ptpde_.Clear();
    ptpde_.Resize(numpde_);



    // Read dimension from mesh file and perform a consistency check
    Integer dim = InFile_->ReadDim();

    for (Integer i=0;i< pdeNames.GetSize();i++) {
      Info->StartProgress("Creating PDE '" + pdeNames[i] + "'");
      if (pdeNames[i] == "electrostatic") 
	ptpde_[i]=new ElecPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

      else if (pdeNames[i] == "mechanic")
	ptpde_[i]=new MechPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

      else if (pdeNames[i] == "acoustic")
	{
	  StdVector<std::string> acouSubType;
	  params->GetList( "subType", acouSubType,"pdeList", "acoustic");
	  if (acouSubType.GetSize())
	    ptpde_[i]=new AcouFlowNoise(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);
	  else
	    ptpde_[i]=new AcousticPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);
	}
      

      else if (pdeNames[i] == "smooth")
	ptpde_[i]=new SmoothPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

      else if (pdeNames[i] == "magnetic") 
	{
	  if (dim == 2)
	    ptpde_[i]=new MagPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);
	  else
	    Error( "Magnetic field calculation currently only possible in 2D!",
		   __FILE__, __LINE__);
	}

      else if (pdeNames[i] == "piezo")
	ptpde_[i]=new PiezoPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

      else if (pdeNames[i] == "mpcci")
	ptpde_[i]=new MpcciPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

//       else if (pdeNames[i] == "acouflownoise")
//       	ptpde_[i]=new AcouFlowNoise(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

      //      else if (pdeNames[i] == "smoothlaplace") 
      //	ptpde_[i]=new SmoothLaPlacePDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_); 


      //      else if (pdeNames[i] == "magnetic") 
      //	if (dim == 3
      //	ptpde_[i]=new MagEdgePDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_); 

      else
	{
	  std::string msg=pdeNames[i]+" - this type of pdes is unknown";
	  Error(msg.c_str(),__FILE__,__LINE__);
	}     

      // Initialize current PDE
      // -> This step has now moved to method InitPDEs
      //ptpde_[i]->Init();

      Info->FinishProgress();

    }


} // end of InitPDE()


  void Domain::CreateCoupledPDE(StdVector<std::string> & sequenceTags) {
    ENTER_FCN( "Domain::InitCoupledPDE" );
  
    std::string errMsg;

    // check if more than one PDEs are defined
    if (numpde_ <= 1) {
      ptcoupledpde_ = NULL;
      return;
    }

    Info->StartProgress("Creating coupling");

    // ================================
    //   Check for iterative coupling
    // ================================
    
    StdVector<BasePDE*> iterCoupledPDEs;
    StdVector<std::string> iterCoupledPDENames;
    StdVector<std::string> methods;

    // Check if all sequenceTags are the same.
    // Currently the tag for the current coupling in a multisequence
    // analysis is determined by simply taking the first tag of
    // all tags for the PDEs, which have to be the same.
    std::string firstTag;
    if (sequenceTags.GetSize() > 1) {
      firstTag = sequenceTags[0];
      for (Integer i=1; i<sequenceTags.GetSize(); i++)
	if (sequenceTags[i] != firstTag) {
	  errMsg = "CreateCoupledPDE: The tags in the <multiSequence> section ";
	  errMsg += "are not all the same in each step.\n Coupling is only ";
	  errMsg += "possible if in each step all the <refTag> are the same!";
	  Error(errMsg.c_str(), __FILE__, __LINE__ );
	}
    }

    params->GetIterCoupledPDEList(iterCoupledPDENames, sequenceTags[0]);
    
    iterCoupledPDEs.Clear();

    // we have all the names of the PDEs which couple iteratively.
    // Now we have to get the according pointers to the PDEs
    for (Integer i=0; i<iterCoupledPDENames.GetSize(); i++)
      for (Integer j=0; j<ptpde_.GetSize(); j++) {
	  if (iterCoupledPDENames[i] == ptpde_[j]->GetName())
	      iterCoupledPDEs.Push_back(ptpde_[j]);
      }
    
    params->GetList( "method", methods);
    for (Integer i=0; i<methods.GetSize(); i++)
      if (methods[i] != "RHS")
	{
	  errMsg  = "Domain::InitCoupledPDE: Methode '";
	  errMsg += methods[i];
	  errMsg += "' not implemented for iterative coupling";
	  Error(errMsg.c_str(), __FILE__, __LINE__);
	} 
    
    
    CoupledPDEDef * CouplingDef = new CoupledPDEDef(ptgrid_, ptBCs_);

    // create coupling objects 
    CouplingDef->CreateCoupling(orderedpdes_,couplings_, iterCoupledPDEs);
    
    // create new iterative coupeld PDE
    ptcoupledpde_ = new IterCoupledPDE(orderedpdes_,couplings_,
				       ptgrid_, ptBCs_,InFile_,
				       OutFile_, sequenceTags[0]);
    
    delete CouplingDef;	


    // ================================
    //   Check for direct coupling
    // ================================

    // -- Not implemented yet --

    Info->FinishProgress();
  }

  void Domain::ResetPDEs()
  {
    ENTER_FCN( "Domain::ResetDEs" );


    // Delete single pdes_
    for (Integer iPDE=0; iPDE<numpde_; iPDE++)
      if (ptpde_[iPDE])
	delete ptpde_[iPDE];


    if (ptcoupledpde_ != NULL)
      delete ptcoupledpde_;
  }


  void Domain::PrintGrid(const Integer level) {
    ENTER_FCN( "Domain::PrintGrid" );
    OutFile_->Init(ptgrid_);
    OutFile_->WriteGrid(level);
  }


  void Domain::SetSubdomains() {
    ENTER_FCN( "Domain::SetSubdomains" );
    Error( "Domain:SetSubdomains: Not implemented!", 
	   __FILE__, __LINE__);
  }


  void Domain::Update(const Integer level) {
    ENTER_FCN( "Domain::Update" );
 
    ptBCs_->Update(ptgrid_);

    // Init AlgSystem
    UpdateAlgSys(level);
  }


  void Domain::UpdateAlgSys(const Integer level) {
    ENTER_FCN( "Domain::UpdateAlgSys"  );

    numlevel_ ++;
 
    //set the algebraic systems and read material data
    for (int i=0;i< numpde_;i++)
      {
	ptpde_[i]->DeleteAlgSys(i);
	ptpde_[i]->Reset();
	ptpde_[i]->SetAlgSys(); 
      }
  }

}
