#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>

#include "magEdgePDE.hh"
#include <DataInOut/Unverg/outUnverg.hh>
#include <Forms/forms_header.hh>
#include <Forms/elecfieldop.hh>
#include <Forms/elecforceop.hh>
#include <Estimator/spaceerror.hh>
 
namespace CoupledField
{

MagEdgePDE::MagEdgePDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
		 WriteResults *aptOut)
:BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering MagEdgePDE::MagEdgePDE " << std::endl;
#endif

  dofspernode_ =1;
  dofsperedge_ =1;

  SetMatrixFactors();

  pdename_    = "magnetic";
  pdematerialclass_ = "magnetic";
  
  conf->getsubdompde(subdoms_,pdename_);
  ReadBCs(pdename_);

  AssignPDENodeNumbers();
  size_ = NumPDENodes_;

  NumElems_ = ptgrid_->GetMaxnumElem(actlevel_, subdoms_); 
  sol_.Resize(size_);
  sol_.Init(0);

}


void MagEdgePDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering MagEdgePDE::SetMatrixFactors" << std::endl;
#endif
 
  matrix_factor_[0] = 0.0;  // factor for stiffness matrix
  matrix_factor_[1] = 0.0;  // factor for damping matrix
  matrix_factor_[2] = 0.0;  // factor for convection matrix
  matrix_factor_[3] = 0.0;  // factor for mass matrix
}


void MagEdgePDE::DiscreteParamsPDE()
{
#ifdef TRACE
  (*trace) << "entering MagEdgePDE::DiscreteParamsPDE" << std::endl;
#endif

  MatrixType_   = RSCALAREDGE; 
  GraphType_    = EDGEGRAPH;
  SystemMatrix_ = TRUE;

}

void MagEdgePDE::SetAlgSys(const Integer as_sysid)
{
#ifdef TRACE
  (*trace) << "entering  BasePDE::SetAlgSys" << std::endl;
#endif

  as_sysid_ = as_sysid;

  //allocate according algebraic system
  algsys_ = new StandardSystem();

  //set solver parameters  
  SetSolverParameters();
  algsys_->SetCycle(VCYCLE);
  algsys_->SetSmoothType(1);
  algsys_->SetSmoothStepFor(2);
  algsys_->SetSmoothStepBack(2);

  //ask the PDE discrete form
  DiscreteParamsPDE();

  //set the node graph
  SetupMatrixGraph(NumPDENodes_, GraphType_);

  //set the edge graph
  SetupEdgeGraph();

  //allocate the necessary matrices as well as solver and preconditioner
  CreateMatrices_Solver();
}


void MagEdgePDE::SolveStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MagEdgePDE::SolveStepStatic" << std::endl;
#endif

  Integer update = 0;
  Integer job = 1;

  Double * ptsol;

  //compute and assemble element matrices
  // SetupMatrices(level);

  //account for bcs
  //  SetBCs(level,update,0);

  //  algsys_->CalcPrecond(job);
  //  algsys_->Solve();

  //  ptsol = algsys_->GetSolutionVal();

  // save solution
  //  Vector<Double> transsol(NumPDENodes_, ptsol);
  //  sol_=transsol;

#ifdef DEBUG
  //  std::string matFileName = "solMat";
  //  OutFile_->WriteSolMatrix(ptgrid_, level, sol_, matFileName);
#endif
}



void MagEdgePDE::SetupMatrices(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MagEdgePDE::SetupMatrices" << std::endl;
#endif
  
  Matrix<Double> elemmat;  
  Matrix<Double> ptCoord;
  BaseFE         * ptElem;

  if (InfoPrint)
    (*infofile) << " ------------------------- Element matrices --------------- " << std::endl;

  Vector<Integer> connecth, connect_PDE;  
  Integer i, j;

  for (i=0; i<subdoms_.size(); i++)
    {
      //reads eps33 (matrix notation starts with 0)
      Double eps33 = materialData_[i].GetPermittivity(2,2);

      std::vector<Elem*> elemssd;
   
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j < elemssd.size(); j++)
	{  
	  ptElem=elemssd[j]->ptElem;

	  BaseForm * bilinear_stiff = new LaplaceInt(ptElem, eps33);

	  connecth=elemssd[j]->connect;
	  
	  ptgrid_->GetCoordNodesElemMat(connecth, ptCoord, level);

	  // CHANGE connecth
	  Mesh2PDENode(connect_PDE,connecth);

	  // stiffness part
	  bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);
	  
#ifdef DEBUG
	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;
	  (*debug) << elemmat << std::endl;
#endif

	  algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connecth.size(), SYSTEM);
	  
	  delete bilinear_stiff;	  

	}  
      
    }

}


void MagEdgePDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering MagEdgePDE::WriteResultsInFile" << std::endl;
#endif

  Integer step=0;
  Double time=0;

//   Vector<Double> Potentialh, *Eh, *Fh;
//   TransformNodeSolution(Potentialh,sol_);

//   if (OutFile_->IsGMV())
//     OutFile_->WriteSolution(Potentialh,step,time,"magnetic vector potential");
//   else
//     OutFile_->WriteSolution(Potentialh,step,time,"magnetic vector potential");

}


void MagEdgePDE::SetupEdgeGraph()
{
#ifdef TRACE
  (*trace) << "entering MagEdgePDE::SetupEdgeGraph" << std::endl;
#endif

  algsys_->CreateEdgeGraph();

  //////////////////////////////////////////////////////////////////////////////////////////////////

  Integer numedge = algsys_->GetEdgeSize();

  //curently hard coded for tets
  Integer elemsize_edge = 6;

  Integer * epos  = new Integer[elemsize_edge];

  // set the graph - connectivity matrix
  BaseFE * ptElem; 
  Integer nsub, iel, fe_type;
  Vector<Integer> connecth;
  Integer * pos;

  for (nsub=0; nsub<subdoms_.size(); nsub++)
    {
      std::vector<Elem*> elemssd;
      ptgrid_->GetElemSD(elemssd,subdoms_[nsub],actlevel_);

      for (iel=0; iel < elemssd.size(); iel++)
	{  
	  ptElem=elemssd[iel]->ptElem;
	  
	  //Map Mesh Node numbers to PDE node numbers
	  Mesh2PDENode(connecth,elemssd[iel]->connect);
	  fe_type=elemssd[iel]->ptElem->feType();
	  if (fe_type != TET)	
	    Error("Currently just TETs supported for MagEdgePDE",__FILE__,__LINE__);

	  pos = connecth.get();
	  epos[0] = algsys_->GetNode2Edge(pos[3], pos[0]);
	  epos[1] = algsys_->GetNode2Edge(pos[3], pos[1]);
	  epos[2] = algsys_->GetNode2Edge(pos[3], pos[2]);
	  epos[3] = algsys_->GetNode2Edge(pos[1], pos[0]);
	  epos[4] = algsys_->GetNode2Edge(pos[2], pos[0]);
	  epos[5] = algsys_->GetNode2Edge(pos[2], pos[1]);

	  //take the absolute value
	  for (Integer j=0; j<elemsize_edge; j++)
	    epos[j] = abs(epos[j]);

	  algsys_->SetElementPosEdge(epos,elemsize_edge,fe_type);
	}
    }

}


MagEdgePDE::~MagEdgePDE()
{

}

} // end of namespace


