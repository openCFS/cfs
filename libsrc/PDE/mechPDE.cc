#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include <Estimator/actimeerror.hh>
#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>

#include "mechPDE.hh" 
#include <Forms/nLinElastInt.hh>

namespace CoupledField
{

MechPDE::MechPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
		 WriteResults *aptOut)
  :BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc), nonLin_(FALSE)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::MechPDE " << std::endl;
#endif

  SetMatrixFactors();

  pdename_ = "mechanic";
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

  sol_.reshape(dofspernode_,NumPDENodes_);

  if (conf->get_option("nonlin",  pdename_ ))
    {
      std::cout << "NONLIN option set !!" << std::endl << std::endl;
      nonLin_=TRUE;
    }
  
}

void MechPDE::DiscreteParamsPDE()
{
#ifdef TRACE
  (*trace) << "entering MechPDE ::DiscreteParamsPDE" << std::endl;
#endif

  MatrixType_   = RBLOCK;
  GraphType_    = NODEGRAPH; 
  SystemMatrix_ = TRUE;

}


void MechPDE::SetMatrixFactors()
{
#ifdef TRACE
  (*trace) << "entering MechPDE::SetMatrixFactors" << std::endl;
#endif

  matrix_factor_[0] = 0.0;  // factor for stiffness matrix
  matrix_factor_[1] = 0.0;  // factor for damping matrix
  matrix_factor_[2] = 0.0;  // factor for convection matrix
  matrix_factor_[3] = 0.0;  // factor for mass matrix
}


void MechPDE::InitCoupling(PDECoupling * Coupling)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::InitCoupling" << std::endl;
#endif
  
  ptCoupling_ = Coupling;
  
  Array<Double> * val;
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      if (ptCoupling_->GetOutputQuantity(i) == "mechdisplacement")
	{
	  //std::cerr << "MechPDE::InitCoupling" << std::endl;
	  ptCoupling_->SetOutputDim(i, Dim_);
	  ptCoupling_->GetOutputValues(i, val);
	  //std::cerr << "mechdisplacement size = " << (*val).size() << " dim = " << val->dim() << std::endl;
	}
    }
}

  
void MechPDE::SolveStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::SolveStepStatic" << std::endl;
#endif

  Integer job = 1;

  Double * ptsol;

  //compute and assemble element matrices
  if (updateBCs_ != 1)
    SetupMatrices(level);
  

  //account for bcs
  SetBCs(level,updateBCs_,0);
  algsys_->CalcPrecond(job);

  
  updateBCs_ = 1;

   algsys_->Solve();

  ptsol = algsys_->GetSolutionVal();

  // save solution
  Integer k=0;

  for (Integer i=0; i<NumPDENodes_; i++)   
    for (Integer dim=0; dim<dofspernode_; dim++)
      sol_[dim][i] = ptsol[k++];
  
 //  //Initialize matrices in order to get BCs correct
  Integer matrixsystype[5];    
  if (SystemMatrix_     == 1) matrixsystype[0] = SYSTEM;      // memory for the system matrix
  if (StiffnessMatrix_  == 1) matrixsystype[1] = STIFFNESS;   // memory for the stiffness matrix
  if (DampingMatrix_    == 1) matrixsystype[2] = DAMPING;     // memory for the damping matrix
  if (ConvectionMatrix_ == 1) matrixsystype[3] = CONVECTION;  // memory for the convection matrix
  if (MassMatrix_       == 1) matrixsystype[4] = MASS;        // memory for the mass matrix
  
  algsys_->InitRHS();
  algsys_->InitSol();
  
  //for (Integer i=0;i<5;i++)
  //{
       //     if (matrixsystype[i] !=0)
  //algsys_->InitMatrix(i+1);
//   }
}

void MechPDE::CalcOutputCoupling()
{
#ifdef TRACE
  (*trace) << "entering MechPDE::CalcOutputCoupling" << std::endl;
#endif  

  std::string quantity;
  std::vector<Integer> * couplingnodes;
  Array<Double> * values;
  
  // std::cerr << "MechPDE has " << ptCoupling_->GetNumOutputCouplings() << "couplings " << std::endl;

  // loop over all output coupling quantities
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      quantity = ptCoupling_->GetOutputQuantity(i);

      switch(ptCoupling_->GetOutputType(i))
	{

	case NODE:
	  
	  ptCoupling_->GetOutputNodes(i, couplingnodes);
	  ptCoupling_->GetOutputValues(i, values);

	  if (quantity == "mechdisplacement")
	    {
	      NodeSolutionToCoupling(*values, *couplingnodes);
	      //std::cerr << "CouplingTerms: " << *values << std::endl;
	    }
	  
	  break;

	case ELEM:
	  Error("No Element coupling output", __FILE__,__LINE__);
	}

    }

}


void MechPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering MechPDE::WriteResultsInFile" << std::endl;
#endif

  Integer laststepcalc=0;
  Double  lasttimecalc=0;
  Array<Double> DispMesh;
 
  TransformNodeSolution(DispMesh, sol_, PDE2MeshNode_);
   OutFile_->WriteNodeSolution(DispMesh, laststepcalc, lasttimecalc,"displacement");
}


void MechPDE::SetupMatrices(const Integer level)
  {
#ifdef TRACE
    (*trace) << "entering MechPDE::SetupMatrices" << std::endl;
#endif

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
	    connecth = elemssd[j]->connect;


	    Matrix<Double> ptCoord;
	    GetElemCoords(connecth, ptCoord, level);
	    // map connect to PDE node numbers
	    Mesh2PDENode(connect_PDE, connecth, Mesh2PDENode_);


	    // =================================================================
	    // stiffness matrix 
	    // =================================================================
	    AssembleStiffness(ptEl, connect_PDE, ptCoord, materialData_[i]);


	    // =================================================================
	    // mass matrix 
	    // =================================================================
	    AssembleMass(ptEl, connect_PDE, ptCoord, materialData_[i].GetDensity());
	  }
      }

#ifdef DEBUG
    algsys_->Print(SYSTEM);
#endif
    
    
#ifdef TRACE
    (*trace) << "Leaving MechPDE::SetupMatrices" << std::endl;
#endif
  }






// assemble stiffness part of FE-equation
void MechPDE::AssembleStiffness(BaseFE * ptEl, Vector<Integer>& connect_PDE,  Matrix<Double>& ptCoord, MaterialData& actMatData)
{
  Matrix<Double> elemmat;
  BaseForm * bilinear_stiff;
  
  
  if (subType_ == "plainStrain")
    bilinear_stiff = new mechPlainStrainInt(ptEl, actMatData);
  else if (subType_ == "3d")
    bilinear_stiff = new mech3DInt(ptEl, actMatData);
  else 
    {
      std::string errMessg;
      errMessg = "Unknown subtype \"" + subType_ + "\" in mech PDE!";		
      Error(errMessg.c_str(),__FILE__,__LINE__);
    }
  
  // stiffness part
  bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);
  algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), SYSTEM);
  
  delete bilinear_stiff;
  

  if (nonLin_)
    {
      if (subType_ != "3d")
	Error("For nonlin mechanics, up to now only 3d sims supported! ",__FILE__,__LINE__);

      Matrix<Double> elDisp;
      std::cout << "For nonlin Sim: no elem disp defined yet!" << std::endl;
      exit(1);

      nLinMech3dInt * stiff_nonLin1 = NULL;      
      stiff_nonLin1 = new nLinMech3dInt(ptEl, actMatData);

      stiff_nonLin1->setActElemDispl(elDisp);
      stiff_nonLin1->CalcElementMatrix(ptCoord, elemmat);
      algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), SYSTEM);

      /*
	nLinMech3dInt * stiff_nonLin2 = NULL;
	stiff_nonLin2 = new nLinMech3dInt_part2(ptEl, actMatData);      
	stiff_nonLin1->setActElemDispl(elDisp);
	bilinear_stiff_nonLin2->CalcElementMatrix(ptCoord, elemmat);
	algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connecth.size(), SYSTEM);
      */
    }

}






// assemble mass part of FE-equation
void MechPDE::AssembleMass(BaseFE * ptEl, Vector<Integer>& connect_PDE,  Matrix<Double>& ptCoord, Double density)
{

  Matrix<Double> elemmat;
  BaseForm * bilinear_mass  = new MassInt(ptEl, density);
  
  
  // mass part
  bilinear_mass->CalcElementMatrix(ptCoord, elemmat);
  
  Matrix <Double> elemMatMultDof;
  MassMultiDof(elemMatMultDof, elemmat, dofspernode_);
  
  // 	    algsys_->SetElementMatrix(elemMatMultDof.getinarray(), connect_PDE.get(), connecth.size(), MASS);
  algsys_->SetElementMatrix(elemMatMultDof.getinarray(), connect_PDE.get(), connect_PDE.size(), MASS);  

#ifdef DEBUG
	    (*debug) << "Mech3d elemmat: " << std::endl << elemmat << std::endl;
#endif
  
  delete bilinear_mass;
}





void MechPDE::SetBCs(const Integer level, const Integer update, const Double atime)
  {
#ifdef TRACE
    (*trace) << "entering MechPDE::SetBCs" << std::endl;
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
	  
	    if (update==1)
	      algsys_->UpdateDirichlet(j+1, val, SYSTEM);
	    else	    
	      algsys_->SetDirichlet(j+1, Mesh2PDENode_[node-1], val, dof, SYSTEM);
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
    (*trace) << "leaving MechPDE::SetBCs" << std::endl;
#endif
  }



Integer MechPDE::GetNrBCDof(const std::string & dofStartString)
  {
#ifdef TRACE
    (*trace) << "leaving MechPDE::GetNrBCDof" << std::endl;
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
	    Error("Unknown dof-type in homog. BC; substring must start with ux, uy or uz!!",__FILE__,__LINE__);
	else
	  Error("Unknown dof-type in homog. BC; substring must start with ux or uy!!",__FILE__,__LINE__);
    
    return nrActDof;
  }

Boolean MechPDE::HasOutput(std::string output)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::HasOutput" << std::endl;
#endif

  if (output == "mechdisplacement")
    return TRUE;

  return FALSE;
}


} // end namespace CoupledField
