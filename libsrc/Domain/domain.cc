#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "domain.hh"
#include <Domain/grid.hh>
#include <Domain/bcs.hh>
#include <Domain/GridCFS/interface_gridcfs.hh>
#include <AlgebraicSystem/interface_piles.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <PDE/basepde.hh>
#include <CoupledPDE/coupledpdedef.hh>
#include <CoupledPDE/pdecoupling.hh>

#ifdef NETGEN
#include "interface_netgen.hh"
#endif

#ifdef GRIDLIB
#include "interface_gridlib.hh"
#endif

#ifdef ADAPTGRID
#include <Domain/AdaptGrid/interface_adgrid.hh>
#endif

#include <PDE/pdes_header.hh>
#include <CoupledPDE/coupled_pdes_header.hh>

#include <Utils/vector.hh>
#include <Utils/array.hh>

namespace CoupledField
{

Domain:: Domain(FileType * const aptFileType, WriteResults * ptOut, TimeFunc * aptTimeFunc) 
{
#ifdef TRACE
  (*trace) << "entering Domain::Domain" << std::endl;
#endif
  
  numlevel_ = 0;

 InFile_     = aptFileType; 
 OutFile_    = ptOut;
 ptTimeFunc_ = aptTimeFunc;

 Integer i;
 
  // read type of output results from conf-file
 std::string libmesh="cfsgrid";
 conf->ifget("mesh_library",libmesh);

 Integer dim=InFile_->ReadDim();

 // initialize pointer to grid 
 if (dim==2) {

  if (libmesh =="cfsgrid") 
      ptgrid_=new GridInterfaceCFS<2>(InFile_);

#ifdef ADAPTGRID
  else
    if (libmesh== "adaptgrid")
      ptgrid_=new InterfaceAdaptGrid<2>(InFile_);
#endif

    else
      Error("Unknown type of mesh_library in conf-file",__FILE__,__LINE__);
 }
 
 if (dim==3) 
   {
     if (libmesh =="cfsgrid") 
       ptgrid_=new GridInterfaceCFS<3>(InFile_);

#ifdef ADAPTGRID
     else
       if (libmesh== "adaptgrid")
	 ptgrid_=new InterfaceAdaptGrid<3>(InFile_);
#endif
       else
	 Error("Unknown type of mesh_library in conf-file",__FILE__,__LINE__);
   }

 //read in the mesh information
 ptgrid_->Read();

 // allocate an object with an information about boundary condition
 ptBCs_=new BCs(InFile_);

 //read restraints information
 ptBCs_->ReadBCs();

 // Initialize PDEs
 InitPDEs();

 // Initialize Coupled PDEs
 InitCoupledPDE();
 
//set the algebraic systems and read material data
  for (int i=0;i< numpde_;i++)
    {
      ptpde_[i]->SetAlgSys(i); 
      ptpde_[i]->ReadMaterialData();
     } 
}

Domain :: ~Domain()
{
#ifdef TRACE
  (*trace) << "entering Domain::~Domain" << std::endl;
#endif

  if (ptgrid_) delete ptgrid_;
  if (ptBCs_) delete ptBCs_;
  //  if (ptalgsys_) delete ptalgsys_;

}

void Domain :: InitPDEs()
{
#ifdef TRACE
  (*trace) << "entering Domain::InitPDEs" << std::endl;
#endif

  // get numbers of PDEs in domain
  std::vector<std::string> pdes;
  conf->getliststr("list_pdes",pdes);

  numpde_=pdes.size();
  ptpde_.resize(numpde_);


  //allocate all specific PDEs
  //  if (!ptalgsys_) Error("You try to allocate object BasePDE with null pointer to AlgSys");

  for (int i=0;i< pdes.size();i++)
    {
      if (pdes[i] == "acoustic")
  	ptpde_[i]=new AcousticPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

      else if (pdes[i] == "acouflownoise")
 	 ptpde_[i]=new AcouFlowNoise(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

      else if (pdes[i] == "mechanic")
	ptpde_[i]=new MechPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

      else if (pdes[i] == "smooth")
	ptpde_[i]=new SmoothPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);

      else if (pdes[i] == "smoothlaplace") 
	ptpde_[i]=new SmoothLaPlacePDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_); 

      else if (pdes[i] == "electrostatic") 
	ptpde_[i]=new ElecPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_); 

      else if (pdes[i] == "magnetic") 
	ptpde_[i]=new MagEdgePDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_); 

      else
	{
	  std::string msg=pdes[i]+" - this type of pdes is unknown";
	  Error(msg.c_str(),__FILE__,__LINE__);
	}     
    }

 

} // end of InitPDE()


void Domain::InitCoupledPDE()
{
#ifdef TRACE
  (*trace) << "entering Domain::InitCoupledPDE" << std::endl;
#endif
  
  // check if more than one PDEs are defined
  if (numpde_ <= 1)
    {
      ptcoupledpde_ = 0;
      return;
    }


  // Initialize the definiton of coupled PDEs
  CoupledPDEDef * CouplingDef = new CoupledPDEDef(ptgrid_, ptBCs_);
  CouplingDef->CreateCoupling(orderedpdes_,couplings_, ptpde_);
  
  // create new Coupled PDE
  ptcoupledpde_ = new IterCoupledPDE(orderedpdes_,couplings_,ptgrid_, ptBCs_,InFile_,OutFile_);
  ptcoupledpde_->InitCoupling(numlevel_);

  delete CouplingDef;

}



void Domain :: PrintGrid(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Domain::PrintGrid" << std::endl;
#endif

 OutFile_->Init(ptgrid_);
 OutFile_->WriteGrid(level);
}

void Domain :: SetSubdomains()
{
#ifdef TRACE
  (*trace) << "entering Domain::SetSubdomains" << std::endl;
#endif
 ;
}

void Domain::Update(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Domain::Update" << std::endl;
#endif
 
  ptBCs_->Update(ptgrid_);

  // Init AlgSystem
  UpdateAlgSys(level);
   
}

void Domain::UpdateAlgSys(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Domain::UpdateAlgSys" << std::endl;
#endif

  numlevel_ ++;
 
  //set the algebraic systems and read material data
  for (int i=0;i< numpde_;i++)
    {
      ptpde_[i]->DeleteAlgSys(i);
      ptpde_[i]->Reset();
      ptpde_[i]->SetAlgSys(i); 
    } 
  
#ifdef TRACE
  (*trace) << " leaving Domain::UpdateAlgSys " << std::endl;
#endif
}


}
