#ifdef NEWBASEPDE

#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include <DataInOut/Unverg/outUnverg.hh>
#include <DataInOut/GMV/outGMV.hh>
#include <Forms/forms_header.hh>
#include <Estimator/spaceerror.hh>

#include <Forms/nLinElastInt.hh>
#include <DataInOut/WriteInfo.hh>
#include <Driver/analysis.hh>
#include "newmark.hh"

#include "newMechPDE.hh" 

namespace CoupledField
{

MechPDE::MechPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
		 WriteResults *aptOut)
  :BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc), 
   nonLin_(FALSE),
   incStopCrit_(1e-2), 
   residualStopCrit_(1e-3)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::MechPDE " << std::endl;
#endif

  pdename_          = "mechanic";
  pdematerialclass_ = "piezo";
  
  conf->getstr("subtype", subType_, pdename_ );
  
  if (subType_ == "3d")
    dofspernode_ = 3;
  else if (subType_ == "axi")
    {
      isaxi_ = TRUE;
      dofspernode_ = 2;
    }
  else
    // default is planeStrain 
    dofspernode_ = 2;

  conf->getsubdompde(subdoms_,pdename_);
  ReadBCs(pdename_);
  
  AssignPDENodeNumbers(Mesh2PDENode_, PDE2MeshNode_, subdoms_);
  numPDENodes_ = PDE2MeshNode_.size();
  size_        = numPDENodes_ * dofspernode_;

  sol_.reshape(dofspernode_, numPDENodes_);


  if (conf->get_option("reducedInt",  pdename_ ))
    {
      std::cout << "REDUCED INTEGRATION set !!" << std::endl << std::endl;
      reducedInt_=TRUE;
    }

  if (conf->get_option("nonlin",  pdename_ ))
    {
      nonLin_=TRUE;

      conf->ifget("incStopCrit", incStopCrit_, pdename_); // incremental stopping criterion
      conf->ifget("residualStopCrit", residualStopCrit_, pdename_); // residual stopping criterion
    }

  preStressVal_ = 0;
  conf->ifget("preStressVal", preStressVal_, pdename_);

  if (preStressVal_)
    GetDirection(preStressDir_, "preStressDir");


  std::string dampstr;
  conf->ifget("damping",dampstr,pdename_);
  if (dampstr == "rayleigh")
    damping_type_ = RAYLEIGH;
  else
   damping_type_ = NONE;

  conf->ifgetliststr("homoBCdof", homDirichDof_, pdename_);
  conf->ifgetliststr("inhomoBCdof", inhomDirichDof_, pdename_);

  //check for b.c. input data
  if (bcs_hd_.size() != homDirichDof_.size()) 
     Error("Inconsistent definition of homogeneous Dirichlet Boundary Conditions");
  if (bcs_id_.size() != inhomDirichDof_.size()) 
     Error("Inconsistent definition of inhomogeneous Dirichlet Boundary Conditions");
 
  conf->ifgetliststr("loadDof", loadDof_, pdename_);

  //check for load data
  if (bcs_loads_.size() != loadDof_.size()) 
     Error("Inconsistent definition of loads");


  // set analysis parameters
  assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, subdoms_);
  assemble_->SetGraphType(NODEGRAPH);
  assemble_->SetMesh2PDENode(&Mesh2PDENode_);
  assemble_->SetMatrixType(RBLOCK);
  assemble_->SetNumDirichlet(GetNumRestraints(actlevel_));
  assemble_->SetPtrBCs(ptBCs_);

  if (damping_type_)
    assemble_->NeedDampingMatrix();

  ReadMaterialData();
   
  DefineIntegrators(actlevel_);  
}




  void MechPDE::DefineIntegrators(const Integer level)
  {
#ifdef TRACE
  (*trace) << "entering MechPDE::DefineIntegerators" << std::endl;
#endif

  Boolean nonLin = FALSE;

  if (nonLin_ && subType_ != "3d")
    Error("For nonlin mechanics, up to now only 3d sims supported! ",__FILE__,__LINE__);

  for (int actSD = 0; actSD < subdoms_.size(); actSD++)
    {
      BaseForm * bilinearStiff;
      
      if (subType_ == "plainStrain")
	bilinearStiff = new mechPlainStrainInt(materialData_[actSD]);
      else if (subType_ == "3d")
	bilinearStiff = new mech3DInt(materialData_[actSD]);
      else 
	Error("Unknown subtype in mech PDE! ",__FILE__,__LINE__);

	  
      //assemble_->AddIntegrator(bilinearStiff, subdoms_[actSD], STIFFNESS, nonLin);
      assemble_->AddIntegrator(bilinearStiff, subdoms_[actSD], STIFFNESS, nonLin);


      // add mass integrator
      BaseForm * bilinear_mass  = new MassInt(materialData_[actSD].GetDensity(), dofspernode_, isaxi_);
      assemble_->AddIntegrator(bilinear_mass, subdoms_[actSD], MASS, nonLin);
    }
  
  }



Integer MechPDE::
GetBCDof(const std::string dofString)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::GetBCDof " << std::endl;
#endif

  if (dofString == "ux")
    return 1;
  else
    if (dofString == "uy")
      return 2;
    else
      if (dofString == "uz")
	return 3;
      else
	Error("The direction mentioned in the config-file is not implemented! ",__FILE__,__LINE__);
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
      if (dirString == "z")
	dir = Z;
      else
	if (dirString == "radXY")
	  dir = radXY;
	else
	  if (dirString == "radXZ")
	    dir = radXZ;
	  else
	    if (dirString == "radYZ")
	      dir = radYZ;
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


void MechPDE::SetupRHS(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::SetupRHS" << std::endl;
#endif
  
  std::vector<Double> elemVec;  
  Matrix<Double> ptCoord;
  BaseFE         * ptElem;
  std::vector<Elem*> elemssd;
	
  Vector<Integer> connecth, connect_PDE;
  Double pressureVal=0;
  
  conf->ifget("pressureLoad", pressureVal, pdename_); // incremental stopping criterion

  if (pressureVal)
    {
      Directions pressureDir;
      GetDirection(pressureDir, "pressureDir");

      std::vector <std::string> pressureDomain;  //!< name of all subdomains containing coils
      conf->ifgetliststr("pressureDomain", pressureDomain, pdename_);

      std::cout << "pressureDomain " << pressureDomain[0] << std::endl;
      
      ptgrid_->GetElemSD(elemssd, pressureDomain[0], level);
	
      for (Integer j=0; j < elemssd.size(); j++)
	{  
	  ptElem   = elemssd[j]->ptElem;
	  connecth = elemssd[j]->connect;
	  GetElemCoords(connecth, ptCoord, level);
	  
	  // get node numbers belonging to PDE
	  Mesh2PDENode(connect_PDE,connecth,Mesh2PDENode_);
 
	  
	  SurfaceIntLinForm pressure(ptElem, pressureDir);
	  pressure.SetMultiplier(pressureVal);
	  pressure.CalcElemVector(ptCoord, elemVec);
	  std::cout << "!!!!!!!!!!!!! belemNr= " << j+1 <<std::endl;

	  algsys_->SetElementRHS(&elemVec[0], connect_PDE.get(), connect_PDE.size());
	}
    }
  
  
  if (nonLin_)
    for (Integer i=0; i<subdoms_.size(); i++)
      {	
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



void MechPDE::SetBCs(const Integer level, const Integer update, const Double atime)
  {
#ifdef TRACE
    (*trace) << "entering MechPDE::SetBCs" << std::endl;
#endif
    
    Integer node,dof;
    Double val;
    
    // set dirichlet boundary conditions
    Integer i = 0;
    Integer j;
    if (PDEisCoupled_)
      j = couplingBCsCounter_;
    else
      j = 0;

    std::list<Integer> nodes;
    Integer sizebc=bcs_hd_.size();
    
    std::vector<Elem*> edgesBC;  // vector of 1D-elements from mesh-file
    Vector<Integer> connecth;
	  

    // homogeneous boundary conditions
    val = 0;
    for (i=0; i< bcs_hd_.size(); i++) 
     {
	std::string doftype = bcs_hd_[i]; 
	//	dof = GetNrBCDof (doftype.substr(0,2));
	dof = GetBCDof(homDirichDof_[i]);
	Info->WriteHomBC(pdename_, bcs_hd_[i], dof);
	
      
// #ifdef DEBUG
// 	(*debug) << std::endl << " Homog. Dirichlet BC" << std::endl;
// #endif
	nodes=ptBCs_->GetNodesLevel(bcs_hd_[i]);
   
	for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	  {
	    node=*p;

#ifdef DEBUG
 	    (*debug) << " node: " << Mesh2PDENode_[node-1] << " dof:" << dof << " val: " << val << "    global node nr: " << node << " running nr: " << j << std::endl;
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
	//	dof = GetNrBCDof (doftype.substr(0,2));
	dof = GetBCDof(inhomDirichDof_[i]);
      
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

	//	dof = GetNrBCDof (doftype.substr(0,2));
	dof = GetBCDof(loadDof_[i]);
	

// #ifdef DEBUG
// 	(*debug) << " Load BC" << std::endl;
// #endif

	nodes = ptBCs_->GetNodesLevel(bcs_loads_[i], level);

	val = val_loads_[i];
	Info->WriteLoad(pdename_, bcs_loads_[i], val, dof);

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




// ======================================================
// COUPLING SECTION
// ======================================================


void MechPDE::InitCoupling(PDECoupling * Coupling)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::InitCoupling" << std::endl;
#endif
  
  PDEisCoupled_ = TRUE;
  ptCoupling_   = Coupling;
  
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



Boolean MechPDE::HasOutput(std::string output)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::HasOutput" << std::endl;
#endif

  if (output == "mechdisplacement")
    return TRUE;

  return FALSE;
}





// ======================================================
// SOLVING SECTION
// ======================================================

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
  //  if (updateBCs_ != 1)
    //    SetupMatrices(level);  

  assemble_->AssembleMatrices(level);


  //account for bcs
  SetBCs(level,updateBCs_,0);
  SetupRHS(level);
  
  algsys_->CalcPrecond(job);
  
  updateBCs_ = 1;

  algsys_->Solve();

  StoreToSolArray(algsys_->GetSolutionVal());

  // initialize for (eventual) new setup
  // if condition is needed here !!!!!
  //InitMatrices();
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
  
  std::vector<Double> actSol(numPDENodes_ * dofspernode_,0);
  std::vector<Double> solIncrement(numPDENodes_ * dofspernode_,0);


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
      //      SetupMatrices(level);
      assemble_->AssembleMatrices(level);
      
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




      *cla << std::endl << " ======================================================= " << std::endl
		<< " NONLINEAR ITERATION " << iterationCounter << std::endl
		<< " ======================================================= " << std::endl;
      *cla << " === Residual norm:          " << residualL2Norm << std::endl;
      *cla << "     Norm of ext. forces:    " << extForcesL2Norm << std::endl;
      *cla << "     Residual error          " << residualErr << std::endl;

      *cla << " === Incremental sol L2Norm: " << solIncrL2Norm << std::endl;
      *cla << "     Actual solution L2Norm: " << actSolL2Norm << std::endl;
      *cla << "     Incremental error       " << incrementalErr << std::endl;

      
      Info->WriteNonLinIter(pdename_, iterationCounter, residualErr, incrementalErr);
      

      // boolean variable, holds condition if another iteration step is necessary
      performOneMoreStep = 
	(incrementalErr > incStopCrit_)||(residualErr > residualStopCrit_);      
      
      
      if(performOneMoreStep)
	{
	  assemble_->InitMatrices();
	  algsys_->InitSol();
	}
    }while(performOneMoreStep && iterationCounter < maxnumiter_);  
}



void MechPDE::SolveStepTrans(const Integer kstep, const Double asteptime, 
				   const Integer level, const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::SolveStepTrans" << std::endl;
#endif

  lasttimecalc_= asteptime;
  Boolean Recalc=FALSE;

  if (laststepcalc_==kstep && kstep!=0) 
    Recalc=TRUE;
  else 
    laststepcalc_= kstep;

  Double * ptsol;
  Integer update,job;

  //current problem with Array class
  Array<Double> solhelp;
  solhelp.reshape(1, numPDENodes_*dofspernode_);
  Integer k=0;
  for (Integer i=0; i< numPDENodes_; i++)
    for (Integer j=0; j<dofspernode_; j++)
      {
	solhelp[0][k] = sol_[j][i];
	k++;
      }

  //perform predictor step
  TS_alg_->Predictor(solhelp);

  if (kstep==0)
    {
      update = 0;
      job = 1;
      //      SetupMatrices(level);
      assemble_->AssembleMatrices(level);
      
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
      algsys_->InitRHS();
      ComputeRHS(lasttimecalc_);
      TS_alg_->UpdateRHS();
    }
  else if (reset)
    {
      update = 1;
      job    = 1;

      algsys_->InitRHS();
      algsys_->InitMatrix(SYSTEM);
      algsys_->InitMatrix(STIFFNESS);
      algsys_->InitMatrix(MASS);
      if (damping_type_)
        algsys_->InitMatrix(DAMPING);

      algsys_->ConstructEffectiveMatrix(matrix_factor_);
      ComputeRHS(lasttimecalc_);
      TS_alg_->UpdateRHS();
    }
  else
    {
      update = 1;
      job    = 3;
      algsys_->InitRHS();
      ComputeRHS(lasttimecalc_);
      TS_alg_->UpdateRHS();
    };

  SetBCs(level,update,lasttimecalc_);
  algsys_->CalcPrecond(job);
  algsys_->Solve();
  ptsol = algsys_->GetSolutionVal();

  k=0;
  for (Integer i=0; i< numPDENodes_; i++)
    for (Integer dim=0; dim<dofspernode_; dim++)
      {
	sol_[dim][i]  = ptsol[k];
	solhelp[0][k] = ptsol[k];
	k++;
      }

  //perform corrector step  
  TS_alg_->Corrector(solhelp);
}



void MechPDE :: InitTimeStepping(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering AcousticPDE::InitTimeStepping" << std::endl;
#endif

  TS_alg_ = new Newmark(pdename_, algsys_, 1, numPDENodes_*dofspernode_, damping_type_);
  TS_alg_->Init(matrix_factor_, dt);

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
  // to incorporate pressure loads
  SetupRHS(level);      


  // stores rhs vector into extForces and returns that L2-norm
  StoreAlgsysToVec(extForces, algsys_->GetRHSVal() );


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

  const Integer numElems = numPDENodes_ * dofspernode_;
  
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

  const Integer numElems = numPDENodes_ * dofspernode_;
  Double quadSum = 0;
  
  for (Integer i=0; i<numElems; i++)   
    quadSum += pt[i]*pt[i];

  return sqrt(quadSum);
}
  



// ======================================================
// POSTPROCESSING SECTION
// ======================================================


void MechPDE::WriteResultsInFile()
{
#ifdef TRACE
  (*trace) << "entering MechPDE::WriteResultsInFile" << std::endl;
#endif

  Array<Double> DispMesh;
 
  TransformNodeSolution(DispMesh, sol_, PDE2MeshNode_);
  OutFile_->WriteNodeSolution(DispMesh, laststepcalc_, lasttimecalc_,"displacement");
}


} // end namespace CoupledField

#endif //#ifdef NEWBASEPDE
