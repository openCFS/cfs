#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "domain.hh"
#include "interface_gridcfs.hh"
#include "interface_piles.hh"
#include "outGMV.hh"

#ifdef NETGEN
#include "interface_netgen.hh"
#endif

#ifdef GRIDLIB
#include "interface_gridlib.hh"
#endif

#ifdef ADAPTGRID
#include "interface_adgrid.hh"
#endif

#include "therm2dPDE.hh"
#include "acoustic2dPDE.hh"
#include "elecst3dPDE.hh"

namespace CoupledField
{

Domain:: Domain(FileType * const aptFileType, WriteResults * ptOut,  Material * materialdata, TimeFunc * aptTimeFunc) 
{
#ifdef TRACE
  (*trace) << "entering Domain::Domain" << std::endl;
#endif

 InFile_     = aptFileType; 
 OutFile_    = ptOut;
 ptmaterial_ = materialdata;
 ptTimeFunc_ = aptTimeFunc;

  // read type of output results from conf-file
 std::string libmesh;
 conf->get("mesh_library",libmesh);

 Integer dim=InFile_->ReadDim();

 // initialize pointer to grid 
#ifdef GRIDLIB
   if (libmesh =="gridlib") ptgrid_=new InterfaceGridlib<Point2D>(InFile_);
  else
#endif
  if (libmesh =="cfsgrid") 
{
     if (dim==2) ptgrid_=new GridInterfaceCFS<Point2D>(InFile_);
         else ptgrid_=new GridInterfaceCFS<Point3D>(InFile_);
}
#ifdef NETGEN
    else 
  if (libmesh == "netgen")
{
    if (dim==2) ptgrid_=new InterfaceNetGen<Point2D>(InFile_);
        else ptgrid_=new InterfaceNetGen<Point3D>(InFile_);
}
#endif
#ifdef ADAPTGRID
     else
   if (libmesh== "adaptgrid")
{
   if (dim==2) ptgrid_=new InterfaceAdaptGrid<Point2D>(InFile_);
       else ptgrid_=new InterfaceAdaptGrid<Point3D>(InFile_);
}
#endif
   else
     Error("Unknown type of mesh_library in conf-file",__FILE__,__LINE__);

 //read in the mesh information
 ptgrid_->Read();

 // allocate an object with an information about boundary condition
  ptBCs_=new BCs(InFile_);

 //read restraints information
 ptBCs_->ReadBCs();

 //set pointer to Algebraic system
 ptalgsys_ = new AlgSysPILES();
 if (!ptalgsys_) Error("Can't allocate memory for algebraic system Piles");

 // it is important this order of these functions
 InitPDE();
 
 Integer level=0;
 InitAlgSys(level);
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
  if (ptpde_) 
     {
      for (i=0; i<numpde_; i++)
        if (ptpde_[i]) delete ptpde_[i];
      delete [] ptpde_;
     } 
}

void Domain :: InitPDE()
{
#ifdef TRACE
  (*trace) << "entering Domain::InitPDE" << std::endl;
#endif

  // get numbers of PDEs in domain
  std::vector<std::string> pdes;
  conf->getliststr("list_pdes",pdes);

  //allocate all specific PDEs
  if (!ptalgsys_) Error("You try to allocate object BasePDE with null pointer to AlgSys");

  for (int i=0;i< pdes.size();i++)
    {
       if (pdes[i] == "acoustic2d")  { ptpde_[i]=new Acoustic2dPDE(ptalgsys_,ptgrid_,ptmaterial_,ptTimeFunc_,InFile_,OutFile_);}
       else
       if (pdes[i] == "electrostatic3d") {  ptpde_[i]=new Elecst3dPDE(ptalgsys_,ptgrid_,ptmaterial_,ptTimeFunc_,InFile_,OutFile_);
                                         }
       else
    { std::string msg=pdes[i]+" - this type of pdes is unknown";
      Error(msg.c_str(),__FILE__,__LINE__);
    }
      
//       else
//       if (eq == "thermal2d") ptpde_[i]=new Therm2dPDE(ptalgsys_,ptgrid_,ptmaterial_,ptTimeFunc_,InFile_,OutFile_);
   }
}

void Domain :: InitAlgSys(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Domain::InitAlgSys" << std::endl;
#endif

  //check, how much systems are needed and how much matrix graphs
  numsys_   = 1;
  numgraph_ = 1;
  ptalgsys_->InitAlgSys(numsys_, numgraph_);

    //set solver parameters
  Double eps, dampiter;
  Integer maxnumit;
  Integer solvertype;
  Integer precondtype;
  Integer insys;
  Integer numeqcoarse;

  std::cout << " a " << std::endl;

  for (insys=0;insys<numsys_;insys++)
    {
      ptpde_[insys]->SetAlgSys_id(insys);
      ptpde_[insys]->SpecifySolver(solvertype,precondtype,eps,dampiter,maxnumit,numeqcoarse);
      ptalgsys_->SetSolverParameter(insys,eps,dampiter,maxnumit,solvertype,precondtype,  numeqcoarse);
    }

  std::cout << " b " << std::endl;

  //init the algsys-graph
  Integer numnode = ptgrid_->GetMaxnumnodes(level);

  //for each system: first diagonal blocks and then off-diagonalblocks
  for (insys=0;insys<numsys_;insys++)
   {
     ptalgsys_->InitAlgSysGraph(numnode,insys,insys);
   }

 // get the graph - connectivity matrix
  Integer nel, numelem;
  Vector<Integer> connect;

  numelem = ptgrid_->GetMaxnumElem(level);

  std::cout << numelem << std::endl;

  for (insys=0;insys<numsys_;insys++)
    {
      for (nel=0; nel<numelem; nel++)
	{
          ptgrid_->GetConnection(connect, nel, level);
	  ptalgsys_->SetAlgSysGraph(connect.get(),connect.size(),insys,insys);
	}
    }

  //now we can create all the necessary matrices
  Integer matrixtype;
  Integer matrixsystype[5];    
  Integer graphtype; 
  Integer numdofpernode;
  Integer numdirichlets;
  Integer numconstraints;

  for (insys=0;insys<numsys_;insys++)
    {
      ptpde_[insys]->SpecifyMatrices(matrixtype, matrixsystype, graphtype, numdofpernode,  numdirichlets, numconstraints);
      numdirichlets = ptpde_[insys]->GetNumRestraints(ptBCs_, level);

      ptalgsys_->CreateAlgSysMatrices(insys,insys,matrixsystype,matrixtype,graphtype, numdofpernode,  numdirichlets, numconstraints);
    }

  //now reset AlgebraicSystem 
  //matrix_id = 1: system matrix
  Integer matrix_id = 1;
  for (insys=0;insys<numsys_;insys++)
    {
      ptalgsys_->ResetAlgSys(insys,insys,matrix_id);
    }
#ifdef TRACE
  (*trace) << "leaving Domain::IniAlgSys" << std::endl;
#endif
}

void Domain :: PrintGrid(Integer level)
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

void Domain:: Update()
{
#ifdef TRACE
  (*trace) << "entering Domain::Update" << std::endl;
#endif
 
  ptBCs_->Update(ptgrid_);
  
  // Init AlgSystem
  UpdateAlgSys();
    
}

void Domain::UpdateAlgSys()
{
#ifdef TRACE
  (*trace) << "entering Domain::UpdateAlgSys" << std::endl;
#endif

  Integer level=ptgrid_->GetLastLevel();
  std::cout << level << std::endl;

  delete ptalgsys_;
  ptalgsys_ = new AlgSysPILES();
  if (!ptalgsys_) Error("Can't allocate memory for algebraic system Piles");

  InitAlgSys(level);
  
#ifdef TRACE
  (*trace) << " leaving Domain::UpdateAlgSys " << std::endl;
#endif

}

void Domain :: TestGrid()
{
  InterfaceNetGen<Point2D> * ptGrid=new InterfaceNetGen<Point2D>(InFile_); 
  ptGrid->Read();
  Char * name="refine";
  WriteResults * ptInFile=new WriteResultsUnverg(name);

  Vector<Integer> ei;
  ei.Resize(3);
  ei[0]=0;
  ei[1]=1;
  ei[2]=2;
  std::cout << " ok " << std::endl;
 
 Integer e=5; 
  ptGrid->SetRefinementFlag(e);
  std::cout << " we set refinement flags " << std::endl;
  ptGrid->Refine();
  std::cout << " we do refinement " << std::endl;

  ptInFile->Init(ptGrid);
  ptInFile->WriteGrid(0);  

  if (ptInFile) delete ptInFile;
}

}
