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

#include "pdes_header.hh"

namespace CoupledField
{

Domain:: Domain(FileType * const aptFileType, WriteResults * ptOut,  Material * materialdata, TimeFunc * aptTimeFunc) 
{
#ifdef TRACE
  (*trace) << "entering Domain::Domain" << std::endl;
#endif

  newlevel = 0;

 InFile_     = aptFileType; 
 OutFile_    = ptOut;
 ptmaterial_ = materialdata;
 ptTimeFunc_ = aptTimeFunc;

 Integer i;
 for (i=0; i< MAXNUMPDE; i++)
   ptpde_[i]=NULL;

  // read type of output results from conf-file
 std::string libmesh;
 conf->get("mesh_library",libmesh);

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
  for (i=0; i<MAXNUMPDE; i++)
    if (ptpde_[i]) delete ptpde_[i];
}

void Domain :: InitPDE()
{
#ifdef TRACE
  (*trace) << "entering Domain::InitPDE" << std::endl;
#endif

  // get numbers of PDEs in domain
  std::vector<std::string> pdes;
  conf->getliststr("list_pdes",pdes);

  numpde_=pdes.size();

  //allocate all specific PDEs
  if (!ptalgsys_) Error("You try to allocate object BasePDE with null pointer to AlgSys");

  for (int i=0;i< pdes.size();i++)
    {
      if (pdes[i] == "acoustic2d")  { ptpde_[i]=new Acoustic2dPDE(ptalgsys_,ptgrid_,ptmaterial_,ptTimeFunc_,InFile_,OutFile_);}
      else
	if (pdes[i] == "electrostatic3d") {  ptpde_[i]=new Elecst3dPDE(ptalgsys_,ptgrid_,ptmaterial_,ptTimeFunc_,InFile_,OutFile_);
	}     
	else
	  if (pdes[i] == "thermal2d") 
	ptpde_[i]=new Therm2dPDE(ptalgsys_,ptgrid_,ptmaterial_,ptTimeFunc_,InFile_,OutFile_);	 
	  else
	    if (pdes[i]=="electrostatic2d") { ptpde_[i]=new Elecst2dPDE(ptalgsys_,ptgrid_,ptmaterial_,ptTimeFunc_,InFile_,OutFile_); }
	else
	if (pdes[i]=="acoustic3d") { ptpde_[i]=new Acoustic3dPDE(ptalgsys_,ptgrid_,ptmaterial_,ptTimeFunc_,InFile_,OutFile_); }
	else
	 if (pdes[i]=="acou2dflownoise") { ptpde_[i]=new Acou2dFlowNoise(ptalgsys_,ptgrid_,ptmaterial_,ptTimeFunc_,InFile_,OutFile_); }
	    else { std::string msg=pdes[i]+" - this type of pdes is unknown";
	    Error(msg.c_str(),__FILE__,__LINE__);
	    } 	  
    } // end of for

} // end of fnc InitPDE()

void Domain :: InitAlgSys(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering Domain::InitAlgSys" << std::endl;
#endif  

  //check, how much systems are needed and how much matrix graphs
  numsys_   = 1;
  numgraph_ = 1;

  //!
  ptalgsys_->InitAlgSys(numsys_, numgraph_);

    //set solver parameters
  Double eps, dampiter;
  Integer maxnumit;
  Integer solvertype;
  Integer precondtype;
  Integer insys;
  Integer numeqcoarse;
  Double  coarsealpha;

  for (insys=0;insys<numsys_;insys++)
    {
      ptpde_[insys]->SetAlgSys_id(insys);
      ptpde_[insys]->SpecifySolver(solvertype,precondtype,eps,dampiter,maxnumit,numeqcoarse,coarsealpha);
      ptalgsys_->SetSolverParameter(insys,eps,dampiter,maxnumit,solvertype,precondtype,  numeqcoarse, coarsealpha);
    }

  //init the algsys-graph
  Integer numnode = ptgrid_->GetMaxnumnodes(level);
  cout << "numnode:" << numnode << endl;

  Integer matrix_graphtype = NODEGRAPH; //nodal graph

  //for each system: first diagonal blocks and then off-diagonalblocks
  for (insys=0;insys<numsys_;insys++)
   {
     ptalgsys_->InitAlgSysGraph(numnode,insys,insys,matrix_graphtype);
   }

 // get the graph - connectivity matrix
  Integer fe_type;
  Vector<Integer> connect;

  std::vector<Elem*> els;
  std::vector<std::string> * sd=ptgrid_->GetAllSDs();
  Integer isd,iel;

  for (insys=0; insys<numsys_; insys++) {   // loop over systems block
    for (isd=0; isd<sd->size(); isd++) { // loop over subdomains
      ptgrid_->GetElemSD(els,(*sd)[isd],level); 
      for (iel=0; iel < els.size(); iel++) { // loop over elems of subdomains
	connect=els[iel]->connect;
	fe_type=els[iel]->ptElem->feType();
	ptalgsys_->SetAlgSysGraph(connect.get(),connect.size(),fe_type,insys,insys);
      }
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
      numdirichlets = ptpde_[insys]->GetNumRestraints(ptBCs_);
      ptalgsys_->CreateAlgSysMatrices(insys,insys,matrixsystype,matrixtype,graphtype, numdofpernode,  numdirichlets, numconstraints);
    }

  //now reset AlgebraicSystem 
  //matrix_id = 1: system matrix
  Integer matrix_id = 1;
  for (insys=0;insys<numsys_;insys++) {
    ptalgsys_->ResetAlgSys(insys,insys,matrix_id);
  }

#ifdef TRACE
  (*trace) << "leaving Domain::InitAlgSys" << std::endl;
#endif
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

  newlevel ++;
  delete ptalgsys_;

  ptalgsys_=new AlgSysPILES();
  if (!ptalgsys_) Error("Can't allocate memory for algebraic system Piles");

  for (int i=0;i< numpde_;i++) {
    ptpde_[i]->InitPtAlgSys(ptalgsys_);
  }
 
  InitAlgSys(level);
  
#ifdef TRACE
  (*trace) << " leaving Domain::UpdateAlgSys " << std::endl;
#endif
}

}
