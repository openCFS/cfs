#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "domain.hh"
#include <Domain/GridCFS/interface_gridcfs.hh>
#include <AlgebraicSystem/interface_piles.hh>
#include <DataInOut/GMV/outGMV.hh>

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

namespace CoupledField
{

Domain:: Domain(FileType * const aptFileType, WriteResults * ptOut, TimeFunc * aptTimeFunc) 
{
#ifdef TRACE
  (*trace) << "entering Domain::Domain" << std::endl;
#endif

  newlevel = 0;

 InFile_     = aptFileType; 
 OutFile_    = ptOut;
 ptTimeFunc_ = aptTimeFunc;

 Integer i;
 for (i=0; i< MAXNUMPDE; i++)
   ptpde_[i]=NULL;

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

 InitPDEs();
 
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
  for (i=0; i<MAXNUMPDE; i++)
    if (ptpde_[i]) delete ptpde_[i];
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

  //allocate all specific PDEs
  //  if (!ptalgsys_) Error("You try to allocate object BasePDE with null pointer to AlgSys");

  for (int i=0;i< pdes.size();i++)
    {
      if (pdes[i] == "acoustic2d")
  	ptpde_[i]=new Acoustic2dPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);
//       else if (pdes[i] == "electrostatic3d")
// 	ptpde_[i]=new Elecst3dPDE(ptalgsys_,ptgrid_,ptTimeFunc_,InFile_,OutFile_);
//       else if (pdes[i] == "thermal2d")
// 	ptpde_[i]=new Therm2dPDE(ptalgsys_,ptgrid_,ptTimeFunc_,InFile_,OutFile_);	 
//         else if (pdes[i] == "electrostatic2d") 
//        ptpde_[i]=new Elecst2dPDE(ptalgsys_,ptgrid_,ptTimeFunc_,InFile_,OutFile_); 
//        else if (pdes[i] == "acoustic3d")
// 	 ptpde_[i]=new Acoustic3dPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);
//        else if (pdes[i] == "acou2dflownoise")
// 	 ptpde_[i]=new Acou2dFlowNoise(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);
// 	ptpde_[i]=new Elec2dPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_); 
//       else if (pdes[i] == "electric3d") 
      else if (pdes[i] == "mech")
	ptpde_[i]=new MechPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_);
      else if (pdes[i] == "elec") 
	ptpde_[i]=new ElecPDE(ptgrid_,ptBCs_,ptTimeFunc_,InFile_,OutFile_); 
      else
	{
	  std::string msg=pdes[i]+" - this type of pdes is unknown";
	  Error(msg.c_str(),__FILE__,__LINE__);
	}     
    }

  //set the algebraic systems and read material data
  for (int i=0;i< pdes.size();i++)
    {
      ptpde_[i]->SetAlgSys(i); 
      ptpde_[i]->ReadMaterialData();
     }

} // end of InitPDE()


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
