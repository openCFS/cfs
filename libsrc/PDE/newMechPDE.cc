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
#include <Driver/assemble.hh>
#include "newmark.hh"

#include "newMechPDE.hh" 

namespace CoupledField
{

MechPDE::MechPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
		 WriteResults *aptOut)
  :BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc), 
   lambdaMat(NULL),
   mueMat(NULL)
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

  reducedInt_ = FALSE;
  if (conf->get_option("reducedInt",  pdename_ ))
    {
      std::cout << "REDUCED INTEGRATION set !!" << std::endl << std::endl;
      reducedInt_=TRUE;
    }

  nonLin_ = FALSE;
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


  //check for damping model
  std::string dampstr;
  conf->ifget("damping",dampstr,pdename_);
  if (dampstr == "rayleigh")
    damping_type_ = RAYLEIGH;
  else
   damping_type_ = NONE;

  if (damping_type_)
    assemble_->NeedDampingMatrix();


  conf->ifgetliststr("homogenBCDof", homDirichDof_, pdename_);  
  conf->ifgetliststr("inhomogenBCDof", inhomDirichDof_, pdename_);


  // just for consistency with old script
  conf->ifgetliststr("homoBCDof", homDirichDof_, pdename_);
  conf->ifgetliststr("homoBCdof", homDirichDof_, pdename_);
  conf->ifgetliststr("inhomoBCDof", inhomDirichDof_, pdename_);
  conf->ifgetliststr("inhomoBCdof", inhomDirichDof_, pdename_);

  //check for b.c. input data
  if (bcs_hd_.size() != homDirichDof_.size()) 
     Error("Inconsistent definition of homogeneous Dirichlet Boundary Conditions");
  if (bcs_id_.size() != inhomDirichDof_.size()) 
     Error("Inconsistent definition of inhomogeneous Dirichlet Boundary Conditions");
 
  // set assemble parameters
  assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, subdoms_, surfdoms_);
  assemble_->SetGraphType(NODEGRAPH);
  assemble_->SetMesh2PDENode(&Mesh2PDENode_);
  assemble_->SetMatrixType(RBLOCK);
  assemble_->SetNumDirichlet(GetNumRestraints(actlevel_));
  assemble_->SetPtrBCs(ptBCs_);
  assemble_->SetPtr2Sol(&sol_);
  assemble_->SetPtr2TimeFnc(ptTimeFunc_);
  
  ReadMaterialData();
   
  DefineIntegrators(actlevel_);  

  ReadSavings();

}


  MechPDE::~MechPDE()
  {
#ifdef TRACE
    (*trace) << "entering MechPDE::~MechPDE " << std::endl;
#endif

    if  (lambdaMat)
      delete lambdaMat;
    
    if  (mueMat)
      delete mueMat;
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

      // ==============  add stiffness ===========================================

      MaterialData actSDMat(materialData_[actSD]);

      // ==============  add "standard" stiffness ===============================
      if (!reducedInt_)  
	{
	  BaseForm * bilinearStiff = GetStiffIntegrator(actSDMat);
	  IntegratorDescriptor * actIntDescr = new IntegratorDescriptor(bilinearStiff, STIFFNESS);

	  //check for damping
	  if (damping_type_ == RAYLEIGH)    
	      actIntDescr->SetSecondaryMat(DAMPING, actSDMat.GetDampingBeta());
	  
	  assemble_->AddIntegrator(actIntDescr, subdoms_[actSD]);
	}
	  
	
      // ==================  add reduced stiffness ===============================
      else 
	{
	  
	  lambdaMat = new MaterialData(actSDMat);
	  mueMat = new MaterialData(actSDMat);
	  
	  // use a "different" material data set for reduced integration
	  CalcReducedMat(*lambdaMat, *mueMat, actSDMat);	

	  BaseForm * bilinearStiff1 = GetStiffIntegrator(*mueMat);
	  assemble_->AddIntegrator(bilinearStiff1, subdoms_[actSD], STIFFNESS, nonLin);


	  BaseForm * bilinearStiff2 = GetStiffIntegrator(*lambdaMat);
	  assemble_->AddIntegrator(bilinearStiff2, subdoms_[actSD], STIFFNESS, nonLin);
	}


      // ==============  add nonlinear stiffness ===============================
      if (nonLin_)
	{
	  if (subType_ != "3d")
	    Error("For nonlin mechanics, up to now only 3d sims supported! ",__FILE__,__LINE__);
	  
	  // assemble prestress, if in config-file given
// 	  if (preStressVal_)
// 	    AssemblePreStressMat(ptEl, connect_PDE, ptCoord, actMatData, elDisp);
	  
      
	  BaseForm * nLinPart1 =  new nLinMech3dInt_BNonLin(actSDMat);
	  assemble_->AddIntegrator(nLinPart1, subdoms_[actSD], STIFFNESS, nonLin_);
	  

	  BaseForm * nLinPart2 = new nLinMech3dInt_PiolaStress(actSDMat);      
	  assemble_->AddIntegrator(nLinPart2, subdoms_[actSD], STIFFNESS, nonLin_);
	}



      // ==============  add mass ================================================
      double density = actSDMat.GetDensity();
      BaseForm * bilinearMass  = new MassInt(density, dofspernode_, isaxi_);

      IntegratorDescriptor * actIntDescr = new IntegratorDescriptor(bilinearMass, MASS);

      //check for damping (mass part)
      if (damping_type_ == RAYLEIGH)    
	  actIntDescr->SetSecondaryMat(DAMPING, actSDMat.GetDampingAlfa());

      assemble_->AddIntegrator(actIntDescr, subdoms_[actSD]);


      // ==================== RHS ================================================
      if (nonLin_)
	{
	  BaseForm * rhsSource = new nLinMech_linFormInt(actSDMat);
	  assemble_->AddRhsIntegrator(rhsSource, subdoms_[actSD], nonLin_);
	}
      
    }
  }


BaseForm *
MechPDE::GetStiffIntegrator(MaterialData& actSDMat, Boolean reducedInt)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::GetStiffIntegrator " << std::endl;
#endif
  
  BaseForm * bilinearStiff;

  if (subType_ == "plainStrain")
    bilinearStiff = new mechPlainStrainInt(actSDMat);
  else if (subType_ == "axi")
    bilinearStiff = new mechAxiInt(actSDMat);
  else if (subType_ == "3d")
    bilinearStiff = new mech3DInt(actSDMat);
  else 
    Error("Unknown subtype in mech PDE! ",__FILE__,__LINE__);

  return bilinearStiff;
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

  iterCoupledCounter_ = 0;

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

void MechPDE:: PreStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MechPDE:: PreStepStatic" << std::endl;
#endif

  if (PDEisCoupled_)
    {
      algsys_->InitRHS();
      algsys_->InitSol();
    }

}


void MechPDE::StepStaticNonLin(const Integer level)
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

      //      SetupRHS(level);      
      assemble_->AssembleRHS(level);
      

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



void MechPDE :: PostStepStatic(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::PostStepStatic" << std::endl;
#endif

  if (PDEisCoupled_)
    iterCoupledCounter_++;

}


void MechPDE :: InitTimeStepping(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::InitTimeStepping" << std::endl;
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
  //  SetupRHS(level);      
  assemble_->AssembleRHS(level);
  

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

  Array<Double> sol_mesh, solder1_mesh, solder2_mesh;
  Array<Double> sol_der1Array, sol_der2Array;

  if (savesol_ == TRUE && (analysistype_== STATIC || analysistype_== TRANSIENT))
    {
      TransformNodeSolution(sol_mesh, sol_, PDE2MeshNode_);
      OutFile_->WriteNodeSolution(sol_mesh, laststepcalc_, lasttimecalc_,"displacement");
    }

  if (analysistype_== TRANSIENT)
    {
      if (savederiv1_ = TRUE)
	{
	  sol_der1Array = getS1();
	  TransformNodeSolution(solder1_mesh,sol_der1Array,PDE2MeshNode_);
	  OutFile_->WriteNodeSolution(solder1_mesh,laststepcalc_,lasttimecalc_,"velocity");
	}

      if (savederiv2_ == TRUE)
	{
	  sol_der2Array = getS2();
	  TransformNodeSolution(solder2_mesh,sol_der2Array,PDE2MeshNode_);
	  OutFile_->WriteNodeSolution(solder2_mesh,laststepcalc_,lasttimecalc_,"acceleration");
	}
    }
}


} // end namespace CoupledField

#endif //#ifdef NEWBASEPDE
