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
#include <General/environment.hh>
 
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
  NumElems_ = ptgrid_->GetMaxnumElem(actlevel_, subdoms_); 

  bField_ = NULL;
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



void MagEdgePDE::CreateMatrices_Solver()
{
#ifdef TRACE
  (*trace) << "entering  MagEdgePDE::CreateMatrices_Solver" << std::endl;
#endif

  Integer matrixclass;
  Integer matrixsystype[5];    
  Integer graphtype; 
  Integer numconstraints;

  //eval number of dirichlet-edges
  EvalNumDirichlet();

  numconstraints = 0;  // currently not handled
  matrixclass    = MatrixType_;
  graphtype      = GraphType_;

  if (SystemMatrix_    != 1)
    Error("One needs at least a system matrix!",__FILE__,__LINE__);

  if (SystemMatrix_     == 1) matrixsystype[0] = SYSTEM;      // memory for the system matrix
  if (StiffnessMatrix_  == 1) matrixsystype[1] = STIFFNESS;   // memory for the stiffness matrix
  if (DampingMatrix_    == 1) matrixsystype[2] = DAMPING;     // memory for the damping matrix
  if (ConvectionMatrix_ == 1) matrixsystype[3] = CONVECTION;  // memory for the convection matrix
  if (MassMatrix_       == 1) matrixsystype[4] = MASS;        // memory for the mass matrix

   //  SpecifyMatrices(matrixclass, matrixsystype, graphtype, numdofpernode,numdirichlets, 
   //  numconstraints);

  //put to algebraic system
  algsys_->CreateMatrix(matrixsystype, matrixclass, graphtype, dofsperedge_,  numEdgedir_,
			numconstraints);

  //create solver and preconditioner
  algsys_->CreateSolver();
  algsys_->CreatePrecond(matrixclass);

  //now reset AlgebraicSystem 
  algsys_->InitRHS();
  algsys_->InitSol();

  for (Integer i=0;i<5;i++)
    {
      if (matrixsystype[i] !=0)
 	algsys_->InitMatrix(i+1);
    }

}


void MagEdgePDE::SolveStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MagEdgePDE::SolveStepStatic" << std::endl;
#endif

  Integer update = 0;
  Integer job = 4;

  Double * ptsol;

  //compute and assemble element matrices
  SetupMatrices(level);

  // calculate source term
  SetupRHS(level);
  
  //account for bcs
  SetBCs(level,update,0);

  algsys_->CalcPrecond(job);
  algsys_->Solve();

  ptsol = algsys_->GetSolutionVal();

  // save solution
  //  Vector<Double> transsol(NumPDENodes_, ptsol);
  for(Integer i=0; i < size_; i++)
    sol_[i] = ptsol[i];

#ifdef DEBUG
  (*debug) << "SolveStepStatic: \n solution \n";
  for (Integer i=0; i < size_; i++)
    (*debug) << ptsol[i] << " ";
  
  (*debug) << std::endl;  
#endif
}



void MagEdgePDE::SetupMatrices(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MagEdgePDE::SetupMatrices" << std::endl;
#endif
  
  Matrix<Double> elemmat;  
  Matrix<Double> elemmat_aux;  

  Matrix<Double> ptCoord;
  BaseFE         * ptElem;

  if (InfoPrint)
    (*infofile) << " ------------------------- Element matrices --------------- " << std::endl;

  Vector<Integer> connecth, connect_PDE;  
  Integer i, j;
  Double reluctivity, conductivity;

  //curently hard coded for tets
  Integer elemsize_edge = 6;
  std::vector<Integer> epos(elemsize_edge);
  std::vector<Integer> esign(elemsize_edge);

  for (i=0; i<subdoms_.size(); i++)
    {
      reluctivity  = 1.0/materialData_[i].GetPermiability();
      conductivity = materialData_[i].GetConductivity();

#ifdef DEBUG
      (*debug) << " Material-Data for Subdomain: " << i+1 << std::endl;
      (*debug) << "   Permeability: " << 1.0/reluctivity << std::endl;
      (*debug) << "   Conductivity: " << conductivity << std::endl;
#endif

      std::vector<Elem*> elemssd;
   
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j < elemssd.size(); j++)
	{  
	  ptElem=elemssd[j]->ptElem;

	  BaseForm * bilinear_stiff = new CurlCurlEdgeInt(ptElem, reluctivity);
	  BaseForm * lapl           = new LaplaceInt(ptElem, reluctivity);

	  connecth=elemssd[j]->connect;
	  
	  ptgrid_->GetCoordNodesElemMat(connecth, ptCoord, level);

	  // CHANGE connecth
	  Mesh2PDENode(connect_PDE, connecth);

	  //get the edge numbers and their signs
	  GetEdgeNumber(connect_PDE.get(), epos, esign);

	  // stiffness part
	  bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);
	  
	  // correct sign of entries in elemmat due to orientation of edge
	  for(Integer ii=0; ii<elemmat.getSize(); ii++)
	    for(Integer jj=0; jj<elemmat.getSize(); jj++)
	      elemmat[ii][jj] *= esign[ii] * esign[jj];
	  
	      
#ifdef DEBUG
	  (*debug) << "Stiffnessmatrix, ElementNumber  " <<   i << std::endl;
	  (*debug) << elemmat << std::endl;
#endif


	  algsys_->SetElementMatrix(elemmat.getinarray(), &epos[0], epos.size(),SYSTEM);

	  lapl->CalcElementMatrix(ptCoord, elemmat_aux);
	  algsys_->SetAuxElementMatrix(elemmat_aux.getinarray(), 
				       connect_PDE.get(), connecth.size());
	  
	  delete bilinear_stiff;
	  delete lapl;

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

  if (!bField_)
    Error("B-Field needs to be calculated before output!",__FILE__,__LINE__);


  if (OutFile_->IsGMV())
    {
      OutFile_->WriteSolution(sol_, step, time, "magnetic vector potential");
      
      // Write Out Vector Data
      for (ShortInt i=0; i<ptgrid_->GetDim(); i++) 
	{
	  std::ostringstream b_fieldname;
	  b_fieldname << "Bfield " << i;
	  OutFile_->WriteDataOnCell(bField_[i], step, time, b_fieldname.str());
	}
    }
  else
    {
      const Integer dim = 3;
      
      std::string b_fieldname = "magnetic field density";
      Matrix<Double> outMat(bField_[0].size(), dim);
      
      for (Integer i=0; i<ptgrid_->GetDim(); i++) 
	for (Integer j=0; j < bField_[0].size(); j++)
	  outMat[i][j] = bField_[i][j];
	  
 	
      OutFile_->WriteDataOnCell(outMat, step, time, b_fieldname);
    }
}


void MagEdgePDE::SetupEdgeGraph()
{
#ifdef TRACE
  (*trace) << "entering MagEdgePDE::SetupEdgeGraph" << std::endl;
#endif

  algsys_->CreateEdgeGraph();

 //number of unknowns (edges)
  size_ = algsys_->GetEdgeSize();

#ifdef DEBUG
   (*debug) << "NumEdges:" << size_ << std::endl;
#endif

  //now we can allocate the memory for the solution vector
  sol_.Resize(size_);
  sol_.Init(0);

  //curently hard coded for tets
  Integer elemsize_edge = 6;

  Integer * epos  = new Integer[elemsize_edge];

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


void MagEdgePDE::EvalNumDirichlet()
{
#ifdef TRACE
  (*trace) << "entering MagEdgePDE::EvalNumDirichlet" << std::endl;
#endif

  Integer i,j,k;

  //curently hard coded for tets
  Integer elemsize_edge = 6;
  Integer * epos  = new Integer[elemsize_edge];
  Integer * dnode = new Integer[size_];

  //set array to zero
  for (i=0; i<size_; i++)
    dnode[i] = 0;
    
  BaseFE  * ptElem;
  std::vector<Elem*>  SurfD;
  std::vector<std::string> surfDirichlet;

  if (conf->ifgetliststr("SurfeDirichlet", surfDirichlet, "magnetic"))
    
      if (surfDirichlet.size())
	SurfD = ptBCs_->getFacesBC(surfDirichlet[0],actlevel_);
      else
	Error("No Surfaces specified as Dirichlet Boundaries",__FILE__,__LINE__);
  else
    Error("No Surfaces specified as Dirichlet Boundaries",__FILE__,__LINE__);

  //hard coded: surface elements of tets are triangles
  Integer surfelemdim = 3;
  Vector<Integer> connecth;
  Integer * pos;

  for (Integer iel=0; iel< SurfD.size(); iel++)
    {
      ptElem=SurfD[iel]->ptElem;

      //Map Mesh Node numbers to PDE node numbers
      Mesh2PDENode(connecth,SurfD[iel]->connect);
      pos = connecth.get();

      epos[0] = algsys_->GetNode2Edge(pos[1], pos[0]);
      epos[1] = algsys_->GetNode2Edge(pos[2], pos[0]);
      epos[2] = algsys_->GetNode2Edge(pos[2], pos[1]);

      for (j=0; j<surfelemdim; j++)
	dnode[abs(epos[j])-1] = abs(epos[j]);

    }

  numEdgedir_=0;
  for (i=0; i<size_; i++)
    {
      if (dnode[i] != 0)
	  numEdgedir_++;
    }

  EdgeDir_ = new Integer[numEdgedir_];

  k=0;
  for (i=0; i<size_; i++)
    {
      if (dnode[i])
	{
	  EdgeDir_[k] = dnode[i];
	  k++;
	}
    }

#ifdef DEBUG
  (*debug) << "MagEdge: dirichlet edges " << numEdgedir_ << std::endl;
  for (i=0; i<numEdgedir_; i++)
    (*debug) <<  "   " << EdgeDir_[i] << std::endl;
#endif

  delete [] dnode;
}

void   MagEdgePDE::SetBCs(const Integer level, const Integer update, const Double time)
{
#ifdef TRACE
  (*trace) << "entering   MagEdgePDE::SetBCs" << std::endl;
#endif

  //currently just homogeneous BC supported

  Double  val;
  Integer edge;

  for (Integer i=0; i<numEdgedir_; i++)
    {  
      edge = EdgeDir_[i];
      val=0; 
      algsys_->SetDirichlet(i+1, edge, val, dofsperedge_, SYSTEM);
    }

#ifdef TRACE
  (*trace) << " leaving  MagEdgePDE::SetBCs " << std::endl;
#endif 
}

void MagEdgePDE::GetEdgeNumber(Integer *pos, std::vector<Integer>& epos, std::vector<Integer>& esign)
{
#ifdef TRACE
  (*trace) << "entering MagEdgePDE::GetEdgeNumber" << std::endl;
#endif

      epos[0] = algsys_->GetNode2Edge(pos[3], pos[0]);
      epos[1] = algsys_->GetNode2Edge(pos[3], pos[1]);
      epos[2] = algsys_->GetNode2Edge(pos[3], pos[2]);
      epos[3] = algsys_->GetNode2Edge(pos[0], pos[1]);
      epos[4] = algsys_->GetNode2Edge(pos[0], pos[2]);
      epos[5] = algsys_->GetNode2Edge(pos[1], pos[2]);

      Integer elemsize_edge = 6;
      for (Integer j=0; j<elemsize_edge; j++)
	{
	  esign[j] = epos[j]/abs(epos[j]);
	  epos[j]  = abs(epos[j]);
	}
#ifdef DEBUG
      (*debug) << "Edge numbers  of Element with sign" << std::endl;
      for (Integer i=0; i<elemsize_edge; i++)
	(*debug) << "Edge: " << epos[i] << "  Sign: " << esign[i] << std::endl;
#endif
}


void MagEdgePDE::SetupRHS(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MagEdgePDE::SetupRHS" << std::endl;
#endif

//   std::vector<std::string> coil;
//   if (conf->ifgetliststr("Coil", coil, "magnetic"))
    
  const Double current = 1;
  const Integer currDirection = 3;
  
  std::vector<Double> elemVec;  
  Matrix<Double> ptCoord;
  BaseFE         * ptElem;

  Vector<Integer> connecth, connect_PDE;  
  Integer i, j;
  Double reluctivity, conductivity;

  //curently hard coded for tets
  Integer elemsize_edge = 6;
  std::vector<Integer> epos(elemsize_edge);
  std::vector<Integer> esign(elemsize_edge);

  for (i=0; i<subdoms_.size(); i++)
    {
      std::vector<Elem*> elemssd;
   
      ptgrid_->GetElemSD(elemssd,subdoms_[i],level);

      for (j=0; j < elemssd.size(); j++)
	{  
	  ptElem=elemssd[j]->ptElem;

	  BaseForm * rhsSource = new LinearEdgeInt(ptElem, current, currDirection);

	  connecth=elemssd[j]->connect;
	  
	  ptgrid_->GetCoordNodesElemMat(connecth, ptCoord, level);

	  // CHANGE connecth
	  Mesh2PDENode(connect_PDE,connecth);

	  //get the edge numbers and their signs
	  GetEdgeNumber(connect_PDE.get(), epos, esign);


#ifdef DEBUG
	  (*debug) << " connecth \n"  << connecth << std::endl;
	  (*debug) << " connect_PDE \n"  << connect_PDE << std::endl
		   << " elemNr " << j << std::endl
		   << " epos " << epos << std::endl;
#endif
	  // stiffness part
	  rhsSource->CalcElemVector(ptCoord, elemVec);
	  
	  // correct sign of entries in elemVec due to orientation of edge
	  for(Integer ii=0; ii<elemVec.size(); ii++)
	    elemVec[ii] *= esign[ii];
	  
	      
#ifdef DEBUG
	  (*debug) << "RHS, ElementNumber  " << j << std::endl;

	  (*debug) << elemVec << std::endl;
#endif

	  algsys_->SetElementRHS(&elemVec[0], &epos[0], epos.size());

	  delete rhsSource;
	} 
    }
}


void MagEdgePDE::PostProcess(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MagEdgePDE::PostProcess" << std::endl;
#endif  

  ShortInt Dim = ptgrid_->GetDim();
  CurlEdgeOp * magField = new CurlEdgeOp(ptgrid_, &Mesh2PDENode_, &sol_, level, algsys_);
 
  // ######### Calculation of the magnetic field #########

  std::vector<Double> lCoord(Dim, 1./4);

  std::vector<Elem*> elemssd;
  Vector<Double> elemMagField;


  bField_ = new Vector<Double>[Dim];

  // Resize solution arrays
  for( Integer i=0; i<Dim; i++)
    bField_[i].Resize(NumElems_);

  
  

  // loop over all subdomains
  for (Integer isd=0; isd<subdoms_.size(); isd++)
    {
      // get vector of Elem of subdomain with color: subdoms[isd]
      ptgrid_->GetElemSD(elemssd, subdoms_[isd], level);
      
      // loop over elements of subdomain
      for (Integer iel=0; iel< elemssd.size(); iel++)
	{

 	  magField->CalcElemCurlEdge(elemMagField, elemssd[iel], lCoord);

#ifdef DEBUG
	  (*debug) << "PostProc elem " << iel << std::endl 
		   << "MagField: \n " << elemMagField << std::endl;
#endif

	  for (Integer k=0; k<Dim; k++)
	    bField_[k][iel] = elemMagField[k];
	}
    }
}




MagEdgePDE::~MagEdgePDE()
{
  ;
}

} // end of namespace


