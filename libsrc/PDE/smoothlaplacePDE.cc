#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>

#include "smoothlaplacePDE.hh" 

namespace CoupledField
{

SmoothLaPlacePDE::SmoothLaPlacePDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
		 WriteResults *aptOut)
:BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering SmoothLaPlacePDE::SmoothLaPlacePDE " << std::endl;
#endif

  firstTurn_ = TRUE;

  SetMatrixFactors();

  pdename_ = "smooth";
  pdematerialclass_ = "piezo";
  
  dofspernode_ = 1;

  conf->getsubdompde(subdoms_,pdename_);
  ReadBCs(pdename_);
  
  AssignPDENodeNumbers(Mesh2PDENode_, PDE2MeshNode_, subdoms_);
  NumPDENodes_ = PDE2MeshNode_.size();
  NumElems_ = ptgrid_->GetMaxnumElem(actlevel_, subdoms_);
  size_ = NumPDENodes_*dofspernode_;

  sol_.reshape(Dim_, NumPDENodes_);
  sol_.init();
  tempSol_.reshape(dofspernode_, NumPDENodes_);
  
  conf->ifget("smooth_factor", smooth_factor_, pdename_ );
  
}


void SmoothLaPlacePDE::DiscreteParamsPDE()
{
#ifdef TRACE
  (*trace) << "entering SmoothLaPlacePDE ::DiscreteParamsPDE" << std::endl;
#endif

  MatrixType_   = RSCALAR;
  GraphType_    = NODEGRAPH; 
  SystemMatrix_ = TRUE;

}


void SmoothLaPlacePDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering SmoothLaPlacePDE::SetMatrixFactors" << std::endl;
#endif

  matrix_factor_[0] = 0.0;  // factor for stiffness matrix
  matrix_factor_[1] = 0.0;  // factor for damping matrix
  matrix_factor_[2] = 0.0;  // factor for convection matrix
  matrix_factor_[3] = 0.0;  // factor for mass matrix
}


void SmoothLaPlacePDE::InitCoupling(PDECoupling * Coupling)
{
#ifdef TRACE
  (*trace) << "entering SmoothLaPlacePDE::InitCoupling" << std::endl;
#endif

  ptCoupling_ = Coupling; 

  // input couplings
  for (Integer i=0; i<ptCoupling_->GetNumInputCouplings(); i++)
    {
      // check for input mechanic displacement
      if (ptCoupling_->GetInputQuantity(i) == "mechdisplacement")
	{
	  numDirichletBCs_ += (dofspernode_ * ptCoupling_->GetInputSize(i));
	  // std::cerr << "Added to SmoothLaPlacePDE::numDirichletBCs " << dofspernode_ * ptCoupling_->GetInputSize(i) << std::endl;
	}
    }
  
  // output couplings
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      // check for output displacement
      if (ptCoupling_->GetOutputQuantity(i) == "smoothdisplacement")
	ptCoupling_->SetOutputDim(i, dofspernode_);
    }
  
  //std::cerr << "SmoothLaPlacePDE: Set NumIDBCs_ to " << numDirichletBCs_ << std::endl;
}
  

void SmoothLaPlacePDE::SolveStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::SolveStepStatic" << std::endl;
#endif

  Integer job = 1;

  Double * ptsol;
  
  
  // Loop over all dimensions
  // as the laplace operator only can
  // solve scalar problems
  
  for (itercount_=0; itercount_<Dim_; itercount_++)
    {
      if (itercount_ != 0)
	CalcInputCoupling();

      //compute and assemble element matrices
      SetupMatrices(level);
      
      //account for bcs
      SetBCs(level,updateBCs_,0);
      
      //updateBCs_ = 1;
      
      algsys_->CalcPrecond(job);
      algsys_->Solve();
      
      firstTurn_ = FALSE;
      
      ptsol = algsys_->GetSolutionVal();
      
      // save solution
      Integer k=0;
      
      for (Integer i=0; i<NumPDENodes_; i++)  
	sol_[itercount_][i] = ptsol[k++];
      
      
      // Initialize matrices in order to get BCs correct
      Integer matrixsystype[5];    
      if (SystemMatrix_     == 1) matrixsystype[0] = SYSTEM;      // memory for the system matrix
      if (StiffnessMatrix_  == 1) matrixsystype[1] = STIFFNESS;   // memory for the stiffness matrix
      if (DampingMatrix_    == 1) matrixsystype[2] = DAMPING;     // memory for the damping matrix
      if (ConvectionMatrix_ == 1) matrixsystype[3] = CONVECTION;  // memory for the convection matrix
      if (MassMatrix_       == 1) matrixsystype[4] = MASS;        // memory for the mass matrix
      
      algsys_->InitRHS();
      algsys_->InitSol();
      
      for (Integer i=0;i<5;i++)
	{
	  if (matrixsystype[i] !=0)
	    algsys_->InitMatrix(i+1);
	}

    } // end for loop

  itercount_ = 0;
      
}


void SmoothLaPlacePDE::CalcOutputCoupling()
{
#ifdef TRACE
  (*trace) << "entering SmoothLaPlacePDE::CalcOutputCoupling" << std::endl;
#endif  

  std::string quantity;
  std::vector<Integer> * couplingnodes;
  Array<Double> * values;

  // loop over all output coupling quantities
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      quantity = ptCoupling_->GetOutputQuantity(i);

      switch(ptCoupling_->GetOutputType(i))
	{

	case NODE:
	  
	  ptCoupling_->GetOutputNodes(i, couplingnodes);
	  ptCoupling_->GetOutputValues(i, values);

	  if (quantity == "smoothdisplacement")
	    {
	      NodeSolutionToCoupling(*values, *couplingnodes);
	    }
	  
	  break;

	case ELEM:
	  Error("No Element coupling output", __FILE__,__LINE__);
	}

    }
}

void SmoothLaPlacePDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering SmoothLaPlacePDE::WriteResultsInFile" << std::endl;
#endif

  Integer laststepcalc=0;
  Double  lasttimecalc=0;
  Array<Double> DispMesh;
 
  TransformNodeSolution(DispMesh, sol_, PDE2MeshNode_);
  OutFile_->WriteNodeSolution(DispMesh, laststepcalc, lasttimecalc,"displacement");
}


void SmoothLaPlacePDE::SetupMatrices(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering SmoothLaPlacePDE::SetupMatrices" << std::endl;
#endif

  
    Matrix<Double> elemmat;
    Matrix<Double> ptCoord;
    Matrix<Double> oldCoord;
    std::vector<Double> LCoord(dofspernode_,0.0);
    Double JDet;
    
    // This is a smaller matrix since it is just for linear 1D elements.
    Matrix<Double> elemmatbnd;
    BaseFE * ptEl;
    Vector<Integer> connecth, connect_PDE;
    std::vector<Elem*> elemssd;
    Integer i, j;
    
    for (i=0; i<subdoms_.size(); i++)
      {	
	ptgrid_->GetElemSD(elemssd,subdoms_[i],level);
	
	for (j=0; j< elemssd.size(); j++)
	  {
	    ptEl = elemssd[j]->ptElem;
	    
	    const Double density = materialData_[i].GetDensity();
	    
	    BaseForm * bilinear_stiff;
	    
	    bilinear_stiff = new LaplaceInt(ptEl, density);
	    
	    connecth=elemssd[j]->connect;

	    GetElemCoords(connecth, ptCoord, level);
 
	    // map connect to PDE node numbers
	    Mesh2PDENode(connect_PDE,connecth,Mesh2PDENode_);

	    
	    // stiffness part
	    bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);

	    algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connecth.size(), SYSTEM);

#ifdef DEBUG
	    (*debug) << "Smooth3d elemmat: " << std::endl << elemmat << std::endl;
#endif


	    delete bilinear_stiff;
	    
	  }
      }
   
#ifdef DEBUG
  algsys_->Print(SYSTEM);
#endif
  
  
#ifdef TRACE
  (*trace) << "Leaving SmoothLaPlacePDE::SetupMatrices" << std::endl;
#endif
}


void SmoothLaPlacePDE::CalcInputCoupling()
{
#ifdef TRACE
  (*trace) << "entering  SmoothLaPlace::CalcInputCoupling" << std::endl;
#endif

  std::vector<Integer> * nodes;
  std::vector<Elem*> * elements;
  Array<Double> * val;
  
  // Reset counter for boundary conditions
  couplingBCsCounter_ = 0;
  
  for (Integer i=0; i<ptCoupling_->GetNumInputCouplings(); i++)
    {

      ptCoupling_ = &ptCoupling_[i];
      ptCoupling_->GetInputValues(i, val);

      //std::cerr << "Calculating Input Coupling" << std::endl;
      //std::cerr << "Number of InputValues: " << val->size() << std::endl;

      switch(ptCoupling_->GetInputType(i))
	{
	  
	case COORD:
	  break;

	case RHS:
	  break;


	case ID_BC:

 	  ptCoupling_->GetInputNodes(i, nodes);
	  for (Integer j=0; j<nodes->size(); j++, couplingBCsCounter_++)
	    {
	      if (updateCouplingBCs_)
		{
		  
		  algsys_->UpdateDirichlet(couplingBCsCounter_+1, (*val)[itercount_][j], SYSTEM);
		}
	      else
		{	
		  algsys_->SetDirichlet(couplingBCsCounter_+1, Mesh2PDENode_[(*nodes)[j] - 1], (*val)[itercount_][j], 1, SYSTEM);
		}
	    }
	  break;
	  
	case MAT:
	  Error("Not implemented yet",__FILE__,__LINE__);
	  break;

	}  // end switch
    } // end for
  
  updateCouplingBCs_ = FALSE;
}


void SmoothLaPlacePDE::SetBCs(const Integer level, const Integer update, const Double atime)
  {
#ifdef TRACE
    (*trace) << "entering SmoothLaPlacePDE::SetBCs" << std::endl;
#endif
    
    Integer node,dof;
    Double val;
    
    // set dirichlet boundary conditions
    Integer i = 0;
    Integer j = couplingBCsCounter_;
    std::list<Integer> nodes;
    Integer sizebc=bcs_hd_.size();
    
    std::vector<Elem*> edgesBC;  // vector of 1D-elements from mesh-file
    Vector<Integer> connecth;
	  

    // homogeneous boundary conditions
    val = 0;
    for (i=0; i< bcs_hd_.size(); i++)
      {
	std::string doftype = bcs_hd_[i]; 
	dof = GetNrBCDof (doftype.substr(0,2));

#ifdef DEBUG
	(*debug) << " Homog. Dirichlet BC" << std::endl;
#endif

	if (dof == (itercount_+1))
	  {
	    nodes=ptBCs_->GetNodesLevel(bcs_hd_[i]);
	    
	    for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	      {
		node=*p;
		
#ifdef DEBUG
		(*debug) << " node: " << Mesh2PDENode_[node-1] << " dof:" << dof << " val: " << val << std::endl;
#endif
		if (update==1)
		  {
		    algsys_->UpdateDirichlet(j+1, val, SYSTEM);
		  }
		else{
		  algsys_->SetDirichlet(j+1, Mesh2PDENode_[node-1], val, 1, SYSTEM);
		}
	      }  
	  }
      }

    
    //inhomogeneous boundary conditions
    for (i=0; i<bcs_id_.size(); i++)
      {
	std::string doftype = bcs_id_[i]; 
	dof = GetNrBCDof (doftype.substr(0,2));
      
#ifdef DEBUG
	(*debug) << " Inhomog. Dirichlet BC" << std::endl;
#endif
	if (dof == (itercount_+1))
	  {
	    nodes = ptBCs_->GetNodesLevel(bcs_id_[i], level);
	    
	    val=val_id_[i];
	    
	    for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	      {
		node=*p;
#ifdef DEBUG
		(*debug) << " node: " << node << " dof:" << dof << " val: " << val << std::endl;
#endif
		//std::cerr << "Smooth: Setting ID_BCs .... " << std::endl;
		
		if (update==1)
		  algsys_->UpdateDirichlet(j+1, val, SYSTEM);
		else	  
	      {
		//std::cerr << "Smooth: Setting ID_BC[" << dof << "][" << node << "] = " << val << std::endl;
		algsys_->SetDirichlet(j+1, Mesh2PDENode_[node-1], val, 1, SYSTEM);
	      }
	      }
	  }
      }

    // load boundary conditions
    for (i=0; i<bcs_loads_.size(); i++)
      {
	std::string doftype = bcs_loads_[i]; 

	dof = GetNrBCDof (doftype.substr(0,2));

#ifdef DEBUG
	(*debug) << " Load BC" << std::endl;
#endif
	if (dof == (itercount_+1))
	  {
	    nodes = ptBCs_->GetNodesLevel(bcs_loads_[i], level);
	    
	    val = val_loads_[i];
	    
	    for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	      {
		node=*p;
		
#ifdef DEBUG
		(*debug) << " node: " << Mesh2PDENode_[node-1] << " dof:" << dof << " val: " << val << std::endl;
#endif
		
		val = val_loads_[i];
		algsys_->SetNodeRHS(val, Mesh2PDENode_[node-1], 1);
		
	      }
	  }    
      }
#ifdef TRACE
    (*trace) << "leaving SmoothLaPlacePDE::SetBCs" << std::endl;
#endif
  }



Integer SmoothLaPlacePDE::GetNrBCDof(const std::string & dofStartString)
  {
#ifdef TRACE
    (*trace) << "entering SmoothLaPlacePDE::GetNrBCDof" << std::endl;
#endif

   Integer nrActDof;
   
   if (dofStartString == "ux")
      nrActDof = 1;
    else 
      if (dofStartString == "uy")
	nrActDof = 2;
      else
	if (dofspernode_ == 3)
	  if (dofStartString == "uz")
	    nrActDof = 3;
	  else
	    {
	    Error("Unknown dof-type in homog. BC; substring must start with ux, uy or uz!!",__FILE__,__LINE__);
	    std::cerr << dofStartString << std::endl;
	    }
	else
	  {
	    std::cerr << dofStartString << std::endl;
	    Error("Unknown dof-type in homog. BC; substring must start with ux or uy!!",__FILE__,__LINE__);
	  }
    
    return nrActDof;
  }

Boolean SmoothLaPlacePDE::HasOutput(std::string output)
{
#ifdef TRACE
  (*trace) << "entering SmoothLaPlacePDE::HasOutput" << std::endl;
#endif
  
  if (output == "smoothdisplacement")
    return TRUE;

  return FALSE;
}


} // end namespace CoupledField
