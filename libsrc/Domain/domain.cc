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

    // Integer i;
 
    // read type of output results from conf-file
#ifndef XMLPARAMS
    std::string libmesh="cfsgrid";
    conf->ifget("mesh_library",libmesh);
#else
    std::string libmesh;
    params->Get( "mesh_library", libmesh );
#endif

    Integer dim=InFile_->ReadDim();

    // initialize pointer to grid 
    if (dim==2) {

      if (libmesh =="cfsgrid") {
	ptgrid_=new GridInterfaceCFS<2>(InFile_);
      }

#ifdef ADAPTGRID
      else if (libmesh== "adaptgrid") {
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
      if (libmesh =="cfsgrid") {
	ptgrid_=new GridInterfaceCFS<3>(InFile_);
      }

#ifdef ADAPTGRID
      else if (libmesh== "adaptgrid") {
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

    // Create PDEs
    CreatePDEs();

    // Create Coupled PDE
    CreateCoupledPDE();
  }


  // **********************
  //   Default destructor
  // **********************
  Domain::~Domain() {
    ENTER_FCN( "Domain::~Domain" );
    if (ptgrid_) delete ptgrid_;
    if (ptBCs_) delete ptBCs_;
  }

  // **************************
  //   Initialization of PDEs
  // **************************

  void Domain::InitPDEs(Integer sequenceStep,
			StdVector<std::string> tags)
  {
    ENTER_FCN( "Domain::InitPDEs") ;

    

    // Initialize single PDE
    for (Integer i=0; i<numpde_; i++) {
      Info->StartProgress( "Initializing PDE " + ptpde_[i]->GetName());
      ptpde_[i]->Init(sequenceStep,tags[i]);
      Info->FinishProgress();
    }

    // initialize coupledPDE
    if (numpde_ > 1)
      {
	Info->StartProgress("Initializing Coupling");
	ptcoupledpde_->InitCoupling(numlevel_);
	Info->FinishProgress();
      }
   
 
    // Initialize algebraic system of 
    // each PDE
    // for (int i=0;i< numpde_;i++)
//       ptpde_[i]->SetAlgSys();

  }

  // **************************
  //   Creation of PDEs
  // **************************
  void Domain::CreatePDEs() {

    ENTER_FCN( "Domain::CreatePDEs" );

    // get numbers of PDEs in domain
    StdVector<std::string> pdes;
#ifndef XMLPARAMS
    conf->getliststr("list_pdes",pdes);
#else
    params->GetPDEList( pdes );
#endif

    numpde_ = pdes.GetSize();
    ptpde_.Resize(numpde_);


    // Read dimension from mesh file and perform a consistency check
    Integer dim = InFile_->ReadDim();
#ifdef XMLPARAMS
    std::string probGeo;
    params->Get( "type", probGeo, "geometry" );
    if ( !(probGeo == "3d"    && dim == 3 ||
	   probGeo == "axi"   && dim == 2 ||
	   probGeo == "plane" && dim == 2 ) ) {
      Info->Error( "Dimensions in parameter file and geometry file do not fit",
		   __FILE__, __LINE__ );
    }
#endif
    
    // Allocate all specific PDEs
    // if (!ptalgsys_) Error("You try to allocate object BasePDE with null
    // pointer to AlgSys");


    for (int i=0;i< pdes.GetSize();i++) {
      Info->StartProgress("Creating PDE '" + pdes[i] + "'");
      if (pdes[i] == "electrostatic") 
	ptpde_[i]=new ElecPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

      else if (pdes[i] == "mechanic")
	ptpde_[i]=new MechPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

      else if (pdes[i] == "acoustic")
  	ptpde_[i]=new AcousticPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

      else if (pdes[i] == "smooth")
	ptpde_[i]=new SmoothPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

      else if (pdes[i] == "magnetic") 
	{
	  if (dim == 2)
	    ptpde_[i]=new MagPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);
	  else
	    Error( "Magnetic field calculation currently only possible in 2D!",
		   __FILE__, __LINE__);
	}

      else if (pdes[i] == "piezo")
	ptpde_[i]=new PiezoPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

      else if (pdes[i] == "acouflownoise")
      	ptpde_[i]=new AcouFlowNoise(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

      //      else if (pdes[i] == "smoothlaplace") 
      //	ptpde_[i]=new SmoothLaPlacePDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_); 


      //      else if (pdes[i] == "magnetic") 
      //	if (dim == 3
      //	ptpde_[i]=new MagEdgePDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_); 

      else
	{
	  std::string msg=pdes[i]+" - this type of pdes is unknown";
	  Error(msg.c_str(),__FILE__,__LINE__);
	}     

      // Initialize current PDE
      // -> This step has now moved to method InitPDEs
      //ptpde_[i]->Init();
      Info->FinishProgress();

    }

 

} // end of InitPDE()


  void Domain::CreateCoupledPDE() {
    ENTER_FCN( "Domain::InitCoupledPDE" );
  
    std::string errMsg;

    // check if more than one PDEs are defined
    if (numpde_ <= 1) {
      ptcoupledpde_ = NULL;
      return;
    }

    Info->StartProgress("Creating coupling");

#ifndef XMLPARAMS
    std::string errMsg;
    errMsg  = "Sorry, you are out-dated!\n";
    errMSg += "Coupling is only supported for .xml-files!";
    Error(errMsg.c_str(), __FILE__, __LINE__);
#else
    // ================================
    //   Check for iterative coupling
    // ================================
    
    StdVector<BasePDE*> iterCoupledPDEs;
    StdVector<std::string> iterCoupledPDENames;
    StdVector<std::string> methods;

    

    params->GetIterCoupledPDEList(iterCoupledPDENames);
    
    iterCoupledPDEs.Clear();

    // we have all the names of the PDEs which couple iteratively.
    // Now we have to get the according pointers to the PDEs
    for (Integer i=0; i<iterCoupledPDENames.GetSize(); i++)
	for (Integer j=0; j<ptpde_.GetSize(); j++)
	  if (iterCoupledPDENames[i] == ptpde_[j]->GetName())
	      iterCoupledPDEs.Push_back(ptpde_[j]);
    
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
				       ptgrid_, ptBCs_,InFile_,OutFile_);
    
    delete CouplingDef;	
#endif


    // ================================
    //   Check for direct coupling
    // ================================

    // -- Not implemented yet --

    Info->FinishProgress();
  }

  void Domain::ResetPDEs()
  {
    ENTER_FCN( "Domain::ResetPDEs" );

    for (Integer iPDE=0; iPDE<numpde_; iPDE++)
      if (ptpde_[iPDE])
	delete ptpde_[iPDE];

    CreatePDEs();

    CreateCoupledPDE();
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
