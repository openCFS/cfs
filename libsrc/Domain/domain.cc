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

void Domain::ArrayTest()
{

  Array<Double> A1(3,6);
  Vector<Double> V1, V2, V3;
  std::vector<Double> v1, v2;

  std::cerr << "Array 1:" << A1 << std::endl;

  V1.Resize(3);
  V1[0] = 3.33;
  V1[1] = 4.5;
  V1[2] = -6.4;

  V2.Resize(3);
  V2[0] = 1;
  V2[1] = 2;
  V2[2] = 3;

  Array<Double> A2(V2);
  
  
  A1 = V1;

  std::cerr << "A1:" << std::endl << A1 << std::endl;
  std::cerr << "A2:" << std::endl << A2 << std::endl;
 
  Array<Double> A3;


  V3 = A1;
  std::cerr << "V3 = " << std::endl << V3 << std::endl;

  A3 = A1 + A2;
  
  std::cerr << "A1 + A2" << std::endl << A3 << std::endl;
  
  A3 = A1 - A2;

  std::cerr << "A1 - A2" << std::endl << A3 << std::endl;

  A3 = A1 * 2;
  
  std::cerr << "A1 * 2 = "<< std::endl << A3 << std::endl;

  A3 = A1 / 2;
  
  std::cerr << "A1 / 2 = " << std::endl << A3 << std::endl;

  std::cerr << "A1 * A2" << A1*A2 << std::endl;

  
  
 
}

Domain:: Domain(FileType * const aptFileType, WriteResults * ptOut, TimeFunc * aptTimeFunc) 
{
#ifdef TRACE
  (*trace) << "entering Domain::Domain" << std::endl;
#endif
  
  //rayTest();
  //ror("End of Arraytest!");

  
  newlevel = 0;

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
#ifdef GRIDLIB
   if (libmesh =="gridlib") ptgrid_=new InterfaceGridlib<2>(InFile_);
  else
#endif
  if (libmesh =="cfsgrid") 
      ptgrid_=new GridInterfaceCFS<2>(InFile_);
#ifdef NETGEN
    else 
  if (libmesh == "netgen")
    ptgrid_=new InterfaceNetGen<2>(InFile_);
#endif
#ifdef ADAPTGRID
     else
   if (libmesh== "adaptgrid")
    ptgrid_=new InterfaceAdaptGrid<2>(InFile_);
#endif
   else
     Error("Unknown type of mesh_library in conf-file",__FILE__,__LINE__);
 }
 
 if (dim==3) {
#ifdef GRIDLIB
   if (libmesh =="gridlib") ptgrid_=new InterfaceGridlib<3>(InFile_);
  else
#endif
  if (libmesh =="cfsgrid") 
      ptgrid_=new GridInterfaceCFS<3>(InFile_);
#ifdef NETGEN
    else 
  if (libmesh == "netgen")
    ptgrid_=new InterfaceNetGen<3>(InFile_);
#endif
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
  if (ptalgsys_) delete ptalgsys_;

  Integer i; 
 
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
//       else if (pdes[i] == "electrostatic3d")
// 	ptpde_[i]=new Elecst3dPDE(ptalgsys_,ptgrid_,ptTimeFunc_,InFile_,OutFile_);
//       else if (pdes[i] == "thermal2d")
// 	ptpde_[i]=new Therm2dPDE(ptalgsys_,ptgrid_,ptTimeFunc_,InFile_,OutFile_);	 
//         else if (pdes[i] == "electrostatic2d") 
//        ptpde_[i]=new Elecst2dPDE(ptalgsys_,ptgrid_,ptTimeFunc_,InFile_,OutFile_); 
//        else if (pdes[i] == "acoustic3d")
// 	 ptpde_[i]=new Acoustic3dPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);
//        else if (pdes[i] == "acouflownoise")
// 	 ptpde_[i]=new AcouFlowNoise(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);
// 	ptpde_[i]=new Elec2dPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_); 
//       else if (pdes[i] == "electric3d") 
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
  ptcoupledpde_->InitCoupling(newlevel);

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
  //  UpdateAlgSys(level);
   
}

// void Domain::UpdateAlgSys(const Integer level)
// {
// #ifdef TRACE
//   (*trace) << "entering Domain::UpdateAlgSys" << std::endl;
// #endif

//   newlevel ++;
//   delete ptalgsys_;

//   ptalgsys_=new AlgSysPILES();
//   if (!ptalgsys_) Error("Can't allocate memory for algebraic system Piles");

//   for (int i=0;i< numpde_;i++) {
//     ptpde_[i]->InitPtAlgSys(ptalgsys_);
//   }
 
//   InitAlgSys(level);
  
// #ifdef TRACE
//   (*trace) << " leaving Domain::UpdateAlgSys " << std::endl;
// #endif
//}


}
