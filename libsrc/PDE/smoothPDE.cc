#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include <Estimator/actimeerror.hh>
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>

#include "smoothPDE.hh" 

namespace CoupledField
{

SmoothPDE::SmoothPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
		 WriteResults *aptOut)
:BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::SmoothPDE " << std::endl;
#endif

  SetMatrixFactors();

  pdename_ = "smooth";
  pdematerialclass_ = "piezo";
  
  
  conf->getstr("subtype", subType_, pdename_ );
  
  if (subType_ == "3d")
    dofspernode_ = 3;
  else
    dofspernode_ = 2;

  conf->getsubdompde(subdoms_,pdename_);
  ReadBCs(pdename_);
  
  AssignPDENodeNumbers(Mesh2PDENode_, PDE2MeshNode_, subdoms_);
  NumPDENodes_ = PDE2MeshNode_.size();
  size_ = NumPDENodes_*dofspernode_;

  sol_.reshape(dofspernode_, NumPDENodes_);
  sol_.init();
  
  smooth_factor_ = 1.0;
  conf->get("smooth_factor", smooth_factor_, pdename_ );
}


void SmoothPDE::DiscreteParamsPDE()
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE ::DiscreteParamsPDE" << std::endl;
#endif

  MatrixType_   = RBLOCK;
  GraphType_    = NODEGRAPH; 
  SystemMatrix_ = TRUE;

}


void SmoothPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::SetMatrixFactors" << std::endl;
#endif

  matrix_factor_[0] = 0.0;  // factor for stiffness matrix
  matrix_factor_[1] = 0.0;  // factor for damping matrix
  matrix_factor_[2] = 0.0;  // factor for convection matrix
  matrix_factor_[3] = 0.0;  // factor for mass matrix
}


void SmoothPDE::InitCoupling(PDECoupling * Coupling)
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::InitCoupling" << std::endl;
#endif

  ptCoupling_ = Coupling; 


  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    if (ptCoupling_->GetOutputQuantity(i) == "smoothdisplacement")
      ptCoupling_->SetOutputDim(i, Dim_);
}
  

void SmoothPDE::SolveStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::SolveStepStatic" << std::endl;
#endif

  Integer update = 0;
  Integer job = 1;

  Double * ptsol;

  //compute and assemble element matrices
  SetupMatrices(level);

  //account for bcs
  SetBCs(level,update,0);

  algsys_->CalcPrecond(job);
  algsys_->Solve();

  ptsol = algsys_->GetSolutionVal();

  // save solution
  Integer k=0;
  
  for (Integer i=0; i<NumPDENodes_; i++)  
    for (Integer dim=0; dim<dofspernode_; dim++)
       sol_[dim][i] = ptsol[k++];
  
}


void SmoothPDE::CalcOutputCoupling()
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::CalcOutputCoupling" << std::endl;
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

void SmoothPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::WriteResultsInFile" << std::endl;
#endif

  Integer laststepcalc=0;
  Double  lasttimecalc=0;
  Array<Double> DispMesh;
 
  TransformNodeSolution(DispMesh, sol_, PDE2MeshNode_);
  OutFile_->WriteNodeSolution(DispMesh, laststepcalc, lasttimecalc,"displacement");
}


  void SmoothPDE::SetupMatrices(const Integer level)
  {
#ifdef TRACE
    (*trace) << "entering SmoothPDE::SetupMatrices" << std::endl;
#endif

    Matrix<Double> elemmat;
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

	    BaseForm * bilinear_mass  = new MassInt(ptEl, density);
	    BaseForm * bilinear_stiff;
	    
	    if (subType_ == "plainStrain")
	      bilinear_stiff = new mechPlainStrainInt(ptEl, materialData_[i]);
	    else if (subType_ == "3d")
	      bilinear_stiff = new mech3DInt(ptEl, materialData_[i]);
	    else 
	      {
		std::string errMessg;
		errMessg = "Unknown subtype \"" + subType_ + "\" in SmoothPDE!";		
		Error(errMessg.c_str(),__FILE__,__LINE__);
	      }


	    connecth=elemssd[j]->connect;

	    Matrix<Double> ptCoord;
	    GetElemCoords(connecth, ptCoord, level);
 
	    // map connect to PDE node numbers
	    Mesh2PDENode(connect_PDE,connecth,Mesh2PDENode_);

	    std::vector<Double> LCoord(dofspernode_,0.0);
	    Double JDet;
	    // stiffness part
	    bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);

	    // Divide elemmat by JDet
	    JDet = ptEl->CalcJacobianDet(LCoord,ptCoord);
	    elemmat = elemmat * (smooth_factor_ * 1.0 / JDet);

	    algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connecth.size(), SYSTEM);

#ifdef DEBUG
	    (*debug) << "Smooth3d elemmat: " << std::endl << elemmat << std::endl;
#endif

 	    // mass part
 	    bilinear_mass->CalcElementMatrix(ptCoord, elemmat);

 	    Matrix <Double> elemMatMultDof;
 	    MassMultiDof(elemMatMultDof, elemmat, dofspernode_);

 	    algsys_->SetElementMatrix(elemMatMultDof.getinarray(), connect_PDE.get(), connecth.size(), MASS);

  
	    delete bilinear_stiff;
	    delete bilinear_mass;
	  }
      }

#ifdef DEBUG
    algsys_->Print(SYSTEM);
#endif
    
    
#ifdef TRACE
    (*trace) << "Leaving SmoothPDE::SetupMatrices" << std::endl;
#endif
  }


  void SmoothPDE::SetBCs(const Integer level, const Integer update, const Double atime)
  {
#ifdef TRACE
    (*trace) << "entering SmoothPDE::SetBCs" << std::endl;
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
	nodes=ptBCs_->GetNodesLevel(bcs_hd_[i]);
   
	for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	  {
	    node=*p;
	    
#ifdef DEBUG
	    (*debug) << " node: " << Mesh2PDENode_[node-1] << " dof:" << dof << " val: " << val << std::endl;
#endif
	    if (update==1)
	      algsys_->UpdateDirichlet(j+1, val, SYSTEM);
	    else
	      algsys_->SetDirichlet(j+1, Mesh2PDENode_[node-1], val, dof, SYSTEM);
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

	nodes = ptBCs_->GetNodesLevel(bcs_id_[i], level);

	val=val_id_[i];

	for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	  {
	    node=*p;
#ifdef DEBUG
	    (*debug) << " node: " << node << " dof:" << dof << " val: " << val << std::endl;
#endif
	    std::cerr << "Smooth: Setting ID_BCs .... " << std::endl;

	    if (update==1)
	      algsys_->UpdateDirichlet(j+1, val, SYSTEM);
	    else	  
	      {
		std::cerr << "Smooth: Setting ID_BC[" << dof << "][" << node << "] = " << val << std::endl;
	      algsys_->SetDirichlet(j+1, Mesh2PDENode_[node-1], val, dof, SYSTEM);
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

	nodes = ptBCs_->GetNodesLevel(bcs_loads_[i], level);

	val = val_loads_[i];

	for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	  {
	    node=*p;

#ifdef DEBUG
	    (*debug) << " node: " << Mesh2PDENode_[node-1] << " dof:" << dof << " val: " << val << std::endl;
#endif
	  
	    val = val_loads_[i];
	    algsys_->SetNodeRHS(val, Mesh2PDENode_[node-1], dof);
	  
	  }
      }    
#ifdef TRACE
    (*trace) << "leaving SmoothPDE::SetBCs" << std::endl;
#endif
  }



Integer SmoothPDE::GetNrBCDof(const std::string & dofStartString)
  {
#ifdef TRACE
    (*trace) << "entering SmoothPDE::GetNrBCDof" << std::endl;
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

bool SmoothPDE::HasOutput(std::string output)
{
#ifdef TRACE
  (*trace) << "entering SmoothPDE::HasOutput" << std::endl;
#endif
  
  if (output == "smoothdisplacement")
    return true;
}


} // end namespace CoupledField
