#include <iostream>
#include <fstream>
#include <string>
#include <vector>

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

    //read in the mesh information
    ptgrid_->Read();
 
    if (PrintGridOnly) {
      PrintGrid(0);
      exit(0);
    }

    // allocate an object with an information about boundary condition
    ptBCs_=new BCs(InFile_);

    //read restraints information
    ptBCs_->ReadBCs();

    InitPDEs();

    // Initialize Coupled PDEs
    InitCoupledPDE();
 
    // Set the algebraic systems and read material data
    for (int i=0;i< numpde_;i++)
      ptpde_[i]->SetAlgSys(i);

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
  //   Initialisation of PDEs
  // **************************
  void Domain::InitPDEs() {

    ENTER_FCN( "Domain::InitPDEs" );

    // get numbers of PDEs in domain
    std::vector<std::string> pdes;
#ifndef XMLPARAMS
    conf->getliststr("list_pdes",pdes);
#else
    params->GetPDEList( pdes );
#endif

    numpde_ = pdes.size();
    ptpde_.resize(numpde_);


    // Read dimension from mesh file and perform a consistency check
    Integer dim = InFile_->ReadDim();
#ifdef XMLPARAMS
    std::string probGeo;
    params->Get( "type", probGeo, "geometry" );
    if ( !(probGeo == "3d"          && dim == 3 ||
	   probGeo == "axi"         && dim == 2 ||
	   probGeo == "plainStrain" && dim == 2 ) ) {
      Info->Error( "Dimensions in parameter file and geometry file do not fit",
		   __FILE__, __LINE__ );
    }
#endif
    
    // Allocate all specific PDEs
    // if (!ptalgsys_) Error("You try to allocate object BasePDE with null
    // pointer to AlgSys");


    for (int i=0;i< pdes.size();i++) {
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
	}

      else if (pdes[i] == "piezo")
	ptpde_[i]=new PiezoPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

//      else if (pdes[i] == "acouflownoise")
// 	 ptpde_[i]=new AcouFlowNoise(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

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
    }

 

} // end of InitPDE()


  void Domain::InitCoupledPDE() {
    ENTER_FCN( "Domain::InitCoupledPDE" );
  
    // check if more than one PDEs are defined
    if (numpde_ <= 1) {
      ptcoupledpde_ = 0;
      return;
    }

    // Initialize the definiton of coupled PDEs
    CoupledPDEDef * CouplingDef = new CoupledPDEDef(ptgrid_, ptBCs_);
    CouplingDef->CreateCoupling(orderedpdes_,couplings_, ptpde_);
  
    // create new Coupled PDE
    ptcoupledpde_ = new IterCoupledPDE(orderedpdes_,couplings_,ptgrid_,
				       ptBCs_,InFile_,OutFile_);
    ptcoupledpde_->InitCoupling(numlevel_);

    delete CouplingDef;
  }


  void Domain::PrintGrid(const Integer level) {
    ENTER_FCN( "Domain::PrintGrid" );
    OutFile_->Init(ptgrid_);
    OutFile_->WriteGrid(level);
  }


  void Domain::SetSubdomains() {
    ENTER_FCN( "Domain::SetSubdomains" );
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
	ptpde_[i]->SetAlgSys(i); 
      }
  }

}
