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
  :BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc), nonLin_(FALSE),
   incStopCrit_(1e-2), residualStopCrit_(1e-3)
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


  if (conf->get_option("reducedInt",  pdename_ ))
    {
      std::cout << "REDUCED INTEGRATION set !!" << std::endl << std::endl;
      reducedInt_=TRUE;
    }

  if (conf->get_option("nonlin",  pdename_ ))
    {
      std::cout << "NONLIN option set !!" << std::endl << std::endl;
      nonLin_=TRUE;

      conf->ifget("incStopCrit", incStopCrit_, pdename_); // incremental stopping criterion
      conf->ifget("residualStopCrit", residualStopCrit_, pdename_); // residual stopping criterion
    }

  conf->ifget("preStressVal", preStressVal_, pdename_);

  if (preStressVal_)
    GetDirection(preStressDir_, "preStressDir");
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
	  ptCoupling_->SetOutputDim(i, Dim_);
	  ptCoupling_->GetOutputValues(i, val);
	}
    }
}


void MechPDE::SolveStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::SolveStepStatic" << std::endl;
#endif

  if (nonLin_)
    SolveStepStaticNonLin(level);
  else
    SolveStepStaticLin(level);
}



  
void MechPDE::SolveStepStaticLin(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::SolveStepStaticLin" << std::endl;
#endif

  const Integer job = 1;

  //compute and assemble element matrices
  if (updateBCs_ != 1)
    SetupMatrices(level);  

  //account for bcs
  SetBCs(level,updateBCs_,0);

  algsys_->CalcPrecond(job);
  
  updateBCs_ = 1;

  algsys_->Solve();

  StoreToSolArray(algsys_->GetSolutionVal());

  // initialize for (eventual) new setup
  // if condition is needed here !!!!!
  InitMatrices();
  algsys_->InitRHS();
  algsys_->InitSol();  
}




void MechPDE::SolveStepStaticNonLin(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::SolveStepStaticNonLin" << std::endl;
#endif

  const Integer job = 1;
  Double extForcesL2Norm;  
  Double residualL2Norm;  
  Double solIncrL2Norm;
  Double actSolL2Norm;
  Double incrementalErr;
  Double residualErr;
  Boolean performOneMoreStep;
  Integer iterationCounter=0;
  
  std::vector<Double> actSol(NumPDENodes_ * dofspernode_,0);
  std::vector<Double> solIncrement(NumPDENodes_ * dofspernode_,0);


  extForcesL2Norm = SetExternalForces(level);


  do
    {
      iterationCounter++;
      std::cout << std::endl << "Nonlinear Mechanics: Perform internal loop nr. " 
		<< iterationCounter << std::endl;      
#ifdef DEBUG
      *debug << std::endl << "====================================================== " << std::endl
	     <<	"Nonlinear Mechanics: Perform internal loop nr. " << iterationCounter << std::endl;      
#endif      
      // setup and solve new system (rhs is already set) =====================
      SetupMatrices(level);
      algsys_->CalcPrecond(job);
      algsys_->Solve();


      // new solution is only an increment of the full solution =============
      StoreAlgsysToVec(solIncrement, algsys_->GetSolutionVal() );
      actSol += solIncrement;
      StoreVecToSolArray(actSol);


      // calculate incremental error ========================================
      solIncrL2Norm = L2Norm(solIncrement);
      actSolL2Norm = L2Norm(actSol);
      if (actSolL2Norm)
	incrementalErr = solIncrL2Norm / actSolL2Norm;
      else
	Error("Zero solution vector!! ", __FILE__,__LINE__);      


      // recalculate RHS with new values to get new residual (f^(k+1))========
      algsys_->InitRHS();
      SetBCs(level, updateBCs_, 0);
      SetupRHS(level);      

      std::vector<Double> actRHS;
      StoreAlgsysToVec(actRHS, algsys_->GetRHSVal() );

      // calculation of residual error =======================================
      residualL2Norm = RhsL2Norm(actRHS); // L2Norm of  ( f_i^(k+1) - f_a )
      residualErr = residualL2Norm / extForcesL2Norm;



      std::cout << " === Residual error    " << residualErr << std::endl;
      std::cout << " === Incremental error " << incrementalErr << std::endl;      

      *infofile << std::endl << " ======================================================= " << std::endl
		<< " NONLINEAR ITERATION " << iterationCounter << std::endl
		<< " ======================================================= " << std::endl;
      //      *infofile << " === actSol: " << std::endl << actSol << std::endl << std::endl;
//       *infofile << " === incrementalSol: " << std::endl << solIncrement << std::endl << std::endl;
//       *infofile << " === actRHS (after eliminating penalty entries): " << std::endl 
// 		<< actRHS << std::endl << std::endl;
      *infofile << " === Residual norm:          " << residualL2Norm << std::endl;
      *infofile << "     Norm of ext. forces:    " << extForcesL2Norm << std::endl;
      *infofile << "     Residual error          " << residualErr << std::endl;

      *infofile << " === Incremental sol L2Norm: " << solIncrL2Norm << std::endl;
      *infofile << "     Actual solution L2Norm: " << actSolL2Norm << std::endl;
      *infofile << "     Incremental error       " << incrementalErr << std::endl;
      


      // boolean variable, holds condition if another iteration step is necessary
      performOneMoreStep = 
	(incrementalErr > incStopCrit_)||(residualErr > residualStopCrit_);      
      
      
      if(performOneMoreStep)
	{
	  InitMatrices();
	  algsys_->InitSol();
	}
    }while(performOneMoreStep && iterationCounter < maxnumiter_);  
}





// sets external forces and returns L2Norm of them
Double MechPDE::SetExternalForces(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::SetExternalForces" << std::endl;
#endif

  Double extForcesL2Norm;  
  std::vector<Double> extForces;


  // account for bcs before first solving step =======================
  SetBCs(level, updateBCs_, 0);


  // stores rhs vector into extForces and returns that L2-norm
  StoreAlgsysToVec(extForces, algsys_->GetRHSVal() );


//   *infofile << " === extForces: " << std::endl << extForces 
// 	    << std::endl << std::endl;


  extForcesL2Norm = L2Norm(extForces);
 
  //  if extForcesL2Norm is 0, no residual norm can be calculated
  if (!extForcesL2Norm)
    Error("Zero external force vector!! ", __FILE__,__LINE__);

  
  return extForcesL2Norm;
}






// calculates L2-norm of RHS regarding dirichlet entries due to penalty formulation by setting them 0
Double MechPDE::RhsL2Norm(std::vector<Double>& actRHS)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::RhsL2Norm" << std::endl;
#endif  
  Integer node, dof;
  
  std::list<Integer> nodes;
  
  
  // Eliminate dirichlet node from RHS (due to penalty formulation)
  for (Integer i=0; i< bcs_hd_.size(); i++)
    {
      std::string doftype = bcs_hd_[i]; 
      dof = GetNrBCDof (doftype.substr(0,2));      
      nodes=ptBCs_->GetNodesLevel(bcs_hd_[i]);
      
      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++)
	{
	  node=*p;

	  actRHS[(Mesh2PDENode_[node-1]-1)*dofspernode_ + dof-1] = 0;
	}
    }

  return L2Norm(actRHS);
}





// stores an algsys_ vector into a std::vector and returns that L2-norm
void MechPDE::StoreAlgsysToVec(std::vector<Double>& stdVec, Double * pt)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::StoreAlgsysToVec" << std::endl;
#endif  

  const Integer numElems = NumPDENodes_ * dofspernode_;
  
  stdVec.resize(numElems);

  for (Integer i=0; i<numElems; i++)   
    stdVec[i] = pt[i];
}


  

// returns that L2-norm of an algsys vector
Double MechPDE::AlgsysL2Norm(Double * pt)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::AlgsysL2Norm" << std::endl;
#endif  

  const Integer numElems = NumPDENodes_ * dofspernode_;
  Double quadSum = 0;
  
  for (Integer i=0; i<numElems; i++)   
    quadSum += pt[i]*pt[i];

  return sqrt(quadSum);
}
  




void MechPDE::CalcOutputCoupling()
{
#ifdef TRACE
  (*trace) << "entering MechPDE::CalcOutputCoupling" << std::endl;
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

	  if (quantity == "mechdisplacement")
	    NodeSolutionToCoupling(*values, *couplingnodes);
	  
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

	MaterialData lambdaMat(materialData_[i]);
	MaterialData mueMat ( materialData_[i]);

	if (reducedInt_)
	  CalcReducedMat(lambdaMat, mueMat, materialData_[i]);
	    

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
	    if (reducedInt_)
	      AssembleStiffness(ptEl, connect_PDE, ptCoord, mueMat);
	    else
	      AssembleStiffness(ptEl, connect_PDE, ptCoord, materialData_[i]);


	    // =================================================================
	    // mass matrix 
	    // =================================================================
	    // AssembleMass(ptEl, connect_PDE, ptCoord, materialData_[i].GetDensity());
	  }

      if (reducedInt_)
	{
	  IntegrationType origIntType = elemssd[0]->ptElem->GetIntType();
	  ptgrid_->SetIntTypeAllElems(GaussOrder1);
	  

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
	      AssembleStiffness(ptEl, connect_PDE, ptCoord, lambdaMat);
	      
	      
	      // =================================================================
	      // mass matrix 
	      // =================================================================
	      // AssembleMass(ptEl, connect_PDE, ptCoord, materialData_[i].GetDensity());
	    }

	  ptgrid_->SetIntTypeAllElems(origIntType);
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
void MechPDE::AssembleStiffness(BaseFE * ptEl, Vector<Integer>& connect_PDE,  
				Matrix<Double>& ptCoord, MaterialData& actMatData)
{
#ifdef TRACE
    (*trace) << "entering MechPDE::AssembleStiffness" << std::endl;
#endif

  Matrix<Double> elemmat;
  BaseForm * bilinear_stiff;
  
  
  if (subType_ == "plainStrain")
    bilinear_stiff = new mechPlainStrainInt(ptEl, actMatData);
  else if (subType_ == "3d")
    bilinear_stiff = new mech3DInt(ptEl, actMatData);
  else 
    Error("Unknown subtype in mech PDE! ",__FILE__,__LINE__);   
  
  
  bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);
  algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), SYSTEM);
  
  delete bilinear_stiff;



  if (nonLin_)
    {
      if (subType_ != "3d")
	Error("For nonlin mechanics, up to now only 3d sims supported! ",__FILE__,__LINE__);

      // returns the solution vector belonging to all nodes of the actual element      
      Matrix<Double> elDisp;
      GetSolOfElement(elDisp, connect_PDE);
      
      // assemble prestress, if in config-file given
      if (preStressVal_)
	AssemblePreStressMat(ptEl, connect_PDE, ptCoord, actMatData, elDisp);

      
      nLinMech3dInt_BNonLin stiff_nonLin1(ptEl, actMatData);
      stiff_nonLin1.setActElemDispl(elDisp);
      stiff_nonLin1.CalcElementMatrix(ptCoord, elemmat);
      algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), SYSTEM);
      

      nLinMech3dInt_PiolaStress stiff_nonLin2(ptEl, actMatData);      
      stiff_nonLin2.setActElemDispl(elDisp);
      stiff_nonLin2.CalcElementMatrix(ptCoord, elemmat);
      algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), SYSTEM);
    }

}

void MechPDE::
GetDirection(Directions& dir, const std::string keyword)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::GetDirection " << std::endl;
#endif

  std::string dirString;
  conf->getstr(keyword, dirString, pdename_); 
  
  if (dirString == "x")
    dir = X;
  else
    if (dirString == "y")
      dir = Y;
    else
      if (dirString == "y")
	dir = Z;
      else
	if (dirString == "radXY")
	  dir = radXY;
	else
	  if (dirString == "radXZ")
	    dir = radXZ;
	  else
	    if (dirString == "radXZ")
	      dir = radXZ;
	    else
	      Error("The direction mentioned in the config-file is not implemented! ",__FILE__,__LINE__);
}






void MechPDE::
AssemblePreStressMat(BaseFE * ptEl, Vector<Integer>& connect_PDE, Matrix<Double>& ptCoord, MaterialData& actMatData, Matrix<Double>& elDisp)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::AssemblePreStressMat " << std::endl;
#endif

  // calc matrix part of prestress ==========================================
  PreStressInt preStressInt(ptEl, actMatData, preStressVal_, preStressDir_);
  
  Matrix<Double> elemmat;
  preStressInt.setActElemDispl(elDisp);
  preStressInt.CalcElementMatrix(ptCoord, elemmat);
  
  algsys_->SetElementMatrix(elemmat.getinarray(), connect_PDE.get(), connect_PDE.size(), SYSTEM);
}


void MechPDE::
AssemblePreStressRHS(BaseFE * ptEl, Vector<Integer>& connect_PDE, Matrix<Double>& ptCoord, MaterialData& actMatData, Matrix<Double>& elDisp)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::AssemblePreStressRHS " << std::endl;
#endif

  // calc RHS prestress =====================================================
  PreStressLinFormInt rhsPreStress(ptEl, actMatData, preStressVal_, preStressDir_);
  
  std::vector<Double> elemVec;  
  rhsPreStress.setActElemDispl(elDisp);
  rhsPreStress.CalcElemVector(ptCoord, elemVec);
  
  // subtract internal forces on rhs from external forces 
  elemVec *= -1;
  
  algsys_->SetElementRHS(&elemVec[0], connect_PDE.get(), connect_PDE.size());
}






void MechPDE::
GetSolOfElement( Matrix<Double>& elDisp, Vector<Integer>& connect_PDE)
{
#ifdef TRACE
    (*trace) << "entering MechPDE::GetSolOfElement" << std::endl;
#endif

  // displacement of element nodes
  elDisp.Resize(dofspernode_, connect_PDE.size());

  for (Integer dim=0; dim<dofspernode_; dim++)
    for(Integer actNode=0; actNode<connect_PDE.size(); actNode++)
      elDisp[dim][actNode] = sol_[dim][connect_PDE[actNode]-1];
}





// assemble mass part of FE-equation
void MechPDE::AssembleMass(BaseFE * ptEl, Vector<Integer>& connect_PDE, 
			   Matrix<Double>& ptCoord, Double density)
{
#ifdef TRACE
    (*trace) << "entering MechPDE::AssembleMass" << std::endl;
#endif

  Matrix<Double> elemmat;
  BaseForm * bilinear_mass  = new MassInt(ptEl, density);
  
  
  // mass part
  bilinear_mass->CalcElementMatrix(ptCoord, elemmat);
  
  Matrix <Double> elemMatMultDof;
  MassMultiDof(elemMatMultDof, elemmat, dofspernode_);
  
  algsys_->SetElementMatrix(elemMatMultDof.getinarray(), connect_PDE.get(), connect_PDE.size(), MASS);

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
      
// #ifdef DEBUG
// 	(*debug) << std::endl << " Homog. Dirichlet BC" << std::endl;
// #endif
	nodes=ptBCs_->GetNodesLevel(bcs_hd_[i]);
   
	for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	  {
	    node=*p;

// #ifdef DEBUG
// 	    (*debug) << " node: " << Mesh2PDENode_[node-1] << " dof:" << dof << " val: " << val << "    global node nr: " << node << std::endl;
// #endif
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
      
// #ifdef DEBUG
// 	(*debug) << " Inhomog. Dirichlet BC" << std::endl;
// #endif

	nodes = ptBCs_->GetNodesLevel(bcs_id_[i], level);

	val=val_id_[i];

	for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	  {
	    node=*p;
// #ifdef DEBUG
// 	      (*debug) << " node: " << node << " dof:" << dof << " val: " << val << std::endl;
// #endif
	  
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

// #ifdef DEBUG
// 	(*debug) << " Load BC" << std::endl;
// #endif

	nodes = ptBCs_->GetNodesLevel(bcs_loads_[i], level);

	val = val_loads_[i];

	for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	  {
	    node=*p;

// #ifdef DEBUG
// 	    (*debug) << " node: " << Mesh2PDENode_[node-1] << " dof:" << dof << " val: "
// 		     << val << "    global node nr: " << node << std::endl;
// #endif
	  
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



void MechPDE::SetupRHS(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::SetupRHS" << std::endl;
#endif
  
  std::vector<Double> elemVec;  
  Matrix<Double> ptCoord;
  BaseFE         * ptElem;
  
  Vector<Integer> connecth, connect_PDE;
  
  
  if (nonLin_)
    for (Integer i=0; i<subdoms_.size(); i++)
      {	
	std::vector<Elem*> elemssd;
	
	ptgrid_->GetElemSD(elemssd,subdoms_[i],level);
	
	for (Integer j=0; j < elemssd.size(); j++)
	  {  
	    ptElem   = elemssd[j]->ptElem;
	    connecth = elemssd[j]->connect;
	    GetElemCoords(connecth, ptCoord, level);
	    
	    // get node numbers belonging to PDE
	    Mesh2PDENode(connect_PDE,connecth,Mesh2PDENode_); 
	    
	    
	    // fetch solution at element nodes
	    Matrix<Double> elDisp;
	    GetSolOfElement(elDisp, connect_PDE);
	    
	    // RHS Integrator
	    nLinMech_linFormInt rhsSource(ptElem, materialData_[i]);
	    rhsSource.setActElemDispl(elDisp);
	    rhsSource.CalcElemVector(ptCoord, elemVec);
	    
	    
	    // subtract internal forces on rhs from external forces 
	    elemVec *= -1;

	    
	    algsys_->SetElementRHS(&elemVec[0], connect_PDE.get(), connect_PDE.size());

	    
	    // assemble prestress, if in config-file given
	    if (preStressVal_)
	      AssemblePreStressRHS(ptElem, connect_PDE, ptCoord, materialData_[i], elDisp);
	  }
      }
}


void MechPDE::CalcReducedMat(MaterialData& lambdaMat, MaterialData& mueMat, MaterialData& mat)
{
  Double lambda, mue;
  mat.GetMatrixData(1,2, lambda);
  mat.GetMatrixData(4,4, mue);
  
  Matrix<Double> * lMechMat = lambdaMat.GetMatrix();
  lMechMat -> Init();

  Matrix<Double> * mueMechMat = mueMat.GetMatrix();
  mueMechMat -> Init();


  for(Integer actRow=0; actRow<3; actRow++)
    {
      for(Integer actCol=0; actCol<3; actCol++)
	(*lMechMat)[actRow][actCol] = lambda;

      (*mueMechMat)[actRow][actRow] = 2*mue;
    }

  for(Integer actRow=3; actRow<6; actRow++)
      (*mueMechMat)[actRow][actRow] = mue;  
}



} // end namespace CoupledField
