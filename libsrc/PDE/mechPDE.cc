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
#include <Elements/basefe.hh>

#include "DataInOut/ParamHandling/BaseParamHandler.hh"

#include "mechPDE.hh"

namespace CoupledField
{

MechPDE::MechPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, 
		 WriteResults *aptOut)
  :BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc), 
   lambdaMat(NULL),
   mueMat(NULL),
   reducedInt_(FALSE),
   preStressVal_(0.0)

{
    ENTER_FCN( "MechPDE::MechPDE" );

    pdename_          = "mechanic";
    pdematerialclass_ = "piezo";

#ifndef XMLPARAMS
    conf->getstr("subtype", subType_, pdename_ );
#else
    params->Get( "subtype", subType_, pdename_ );
#endif

    if (subType_ == "3d")
      {
	dofspernode_ = 3;
	Info->PrintF("", "=== 3D PROBLEM\n");
      }
    else if (subType_ == "axi")
      {
	isaxi_ = TRUE;
	dofspernode_ = 2;
	Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
      }
    else if (subType_ == "plainStrain" )
      {
	dofspernode_ = 2;
	Info->PrintF("", "=== PLAIN STRAIN PROBLEM\n");
      }
    else
      {
	std::string errmsg = "Subtype " + subType_ + " is not defined for";
	errmsg += " PDEs of type " + pdename_ + '\n';
	Info->Error( errmsg, __FILE__, __LINE__ );
      }


#ifndef XMLPARAMS
    conf->getsubdompde(subdoms_,pdename_);
#else
    params->GetList( "name", subdoms_, pdename_, "region" );
#endif

    ReadBCs(pdename_);
  

    AssignPDENodeNumbers(Mesh2PDENode_, PDE2MeshNode_, subdoms_);
    numPDENodes_ = PDE2MeshNode_.size();
    size_        = numPDENodes_ * dofspernode_;

    // Initialize solution class
    sol_->SetNumSolutions(1);
    sol_->SetSolutionType(MECH_DISPLACEMENT);
    sol_->SetNumNodes(numPDENodes_);
    sol_->SetDof(dofspernode_);
    sol_->Init(0.0);

    effectiveMass_ = FALSE;
#ifndef XMLPARAMS
    lineSearch_ = FALSE;
    if (conf->get_option("lineSearch",  pdename_ ))
      lineSearch_ = TRUE;

    effectiveMass_ = FALSE;
    if (conf->get_option("effMass",  pdename_ ))
      effectiveMass_ = TRUE;
#endif


    // Check whether reduced integration is to be performed
#ifndef XMLPARAMS
    reducedInt_ = conf->get_option( "reducedInt",  pdename_ );
#else
    reducedInt_ = params->IsSet( "reducedInt", pdename_ );
#endif

    if ( reducedInt_ == TRUE )
      {
	std::cout << "REDUCED INTEGRATION set !!" << std::endl << std::endl;
	reducedInt_=TRUE;
      }

    // Check whether we have to perform a non-linear computation
#ifndef XMLPARAMS
    nonLin_ = conf->get_option( "nonlin",  pdename_ );
#else
    // ==============================================================
    // NOTE: Currently we can only treat geometric non-linearity and
    //       we assume that for a mechanic PDE all regions either
    //       are linear or non-linear!
    // ==============================================================
    std::vector<std::string> nonLinRegion;
    params->GetList( "nonLinear", nonLinRegion, pdename_, "region" );
    // Should not happen with validating parser, but beware!
    if ( nonLinRegion.size() == 0 ) {
      nonLin_ = FALSE;
    }
    else {
      for ( Integer k = 1; k < nonLinRegion.size(); k++ ) {
	if ( nonLinRegion[k] != nonLinRegion[0] ) {
	  Info->Error( "Non-linearity should be the same for all regions!",
		       __FILE__, __LINE__ );
	}
      }
      nonLin_ = nonLinRegion[0] == "geo" ? TRUE : FALSE;
    }
#endif

    if( nonLin_ == TRUE )
      {
#ifndef XMLPARAMS
	// incremental stopping criterion
	conf->ifget("incStopCrit", incStopCrit_, pdename_);

	// residual stopping criterion
	conf->ifget("residualStopCrit", residualStopCrit_, pdename_);
#else
	// incremental stopping criterion
	params->Get( "incStopCrit", incStopCrit_, pdename_, "nonLinear" );

	// residual stopping criterion
	params->Get( "resStopCrit", incStopCrit_, pdename_, "nonLinear" );
#endif
      }

#ifndef XMLPARAMS
    preStressVal_ = 0;
    conf->ifget("preStressVal", preStressVal_, pdename_);
#else
    params->Get( "preStressVal", preStressVal_, pdename_ );
#endif

    if (preStressVal_)
      GetDirection(preStressDir_, "preStressDir");


    // check for damping model
#ifndef XMLPARAMS
    std::string dampstr;
    conf->ifget("damping",dampstr,pdename_);
    if (dampstr == "rayleigh")
      damping_type_ = RAYLEIGH;
    else
      damping_type_ = NONE;
#else
    if( params->HasValue( "type", "rayleigh", pdename_, "damping" ) )
      {
	damping_type_ = RAYLEIGH;
      }
    else
      {
	damping_type_ = NONE;
      }
#endif

    if (damping_type_)
      assemble_->NeedDampingMatrix();

#ifndef XMLPARAMS
    conf->ifgetliststr("homogenBCDof", homDirichDof_, pdename_);  
    conf->ifgetliststr("inhomogenBCDof", inhomDirichDof_, pdename_);

    // just for consistency with old script
    conf->ifgetliststr("homoBCDof", homDirichDof_, pdename_);
    conf->ifgetliststr("homoBCdof", homDirichDof_, pdename_);
    conf->ifgetliststr("inhomoBCDof", inhomDirichDof_, pdename_);
    conf->ifgetliststr("inhomoBCdof", inhomDirichDof_, pdename_);
#else
    params->GetList( "dof", homDirichDof_  , pdename_, "dirichletHom" );  
    params->GetList( "dof", inhomDirichDof_, pdename_, "dirichletInhom" );  
#endif

    //check for b.c. input data
    if (bcs_hd_.size() != homDirichDof_.size())
      {
	std::string errmsg = "Inconsistent definition of homogeneous ";
	errmsg += "Dirichlet Boundary Conditions\n";
	errmsg += " bcs_hd_.size() = " + Info->GenStr( bcs_hd_.size() );
	errmsg += "\n homDirichDof_.size() = "
	  + Info->GenStr( homDirichDof_.size() ) + '\n';
	Info->Error( errmsg, __FILE__, __LINE__ );
      }
    if (bcs_id_.size() != inhomDirichDof_.size()) 
      {
	std::string errmsg = "Inconsistent definition of inhomogeneous ";
	errmsg += "Dirichlet Boundary Conditions";
	Info->Error( errmsg, __FILE__, __LINE__ );
      }
 

  // set assemble parameters
  assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, subdoms_, surfdoms_);
  assemble_->SetGraphType(NODEGRAPH);
  assemble_->SetMesh2PDENode(&Mesh2PDENode_);

#ifdef USE_OLAS
  assemble_->SetMatrixEntryType(DOUBLE);
  assemble_->SetMatrixStorageType(SPARSE_NONSYM);
#else
  assemble_->SetMatrixType(RBLOCK);
#endif 

  assemble_->SetNumDirichlet(GetNumRestraints(actlevel_));

  assemble_->SetPtrBCs(ptBCs_);
  assemble_->SetPtr2Sol(sol_);
  assemble_->SetPtr2TimeFnc(ptTimeFunc_);

  ReadMaterialData();
    
  DefineIntegrators(actlevel_);  
  
  ReadSavings();

}


  MechPDE::~MechPDE()
  {
    ENTER_FCN( "MechPDE::~MechPDE" );

    if  (lambdaMat)
      delete lambdaMat;
    
    if  (mueMat)
      delete mueMat;
  }
  


  void MechPDE::DefineIntegrators(const Integer level)
  {
    ENTER_FCN( "MechPDE::DefineIntegerators" );

    Boolean nonLin = FALSE;


    for (int actSD = 0; actSD < subdoms_.size(); actSD++)
      {

	// ==============  add stiffness ======================================

	MaterialData actSDMat(materialData_[actSD]);

	// ==============  add "standard" stiffness ===========================
	if (!reducedInt_)  
	  {
	    BaseForm * bilinearStiff = GetStiffIntegrator(actSDMat);
	    IntegratorDescriptor * actIntDescr =
	      new IntegratorDescriptor(bilinearStiff, STIFFNESS);

	    //check for damping
	    if (damping_type_ == RAYLEIGH)    
	      actIntDescr->SetSecondaryMat(DAMPING, actSDMat.GetDampingBeta());
	
	    assemble_->AddIntegrator(actIntDescr, subdoms_[actSD]);
	  }
	  
	
	// ==================  add reduced stiffness ==========================
	else 
	  {
	  
	    lambdaMat = new MaterialData(actSDMat);
	    mueMat = new MaterialData(actSDMat);
	  
	    // use a "different" material data set for reduced integration
	    CalcReducedMat(*lambdaMat, *mueMat, actSDMat);	

	    BaseForm * bilinearStiff1 = GetStiffIntegrator(*mueMat);
	    assemble_->AddIntegrator(bilinearStiff1, subdoms_[actSD],
				     STIFFNESS, nonLin);

	    BaseForm * bilinearStiff2 = GetStiffIntegrator(*lambdaMat);
	    assemble_->AddIntegrator(bilinearStiff2, subdoms_[actSD],
				     STIFFNESS, nonLin);
	  }


	// ==============  add nonlinear stiffness ============================
	if (nonLin_)
	  {
	    BaseForm *nLinPart1 = NULL;
	    BaseForm *nLinPart2 = NULL;
	  
	    if (subType_ == "3d")
	      {	  
		nLinPart1 = new nLinMech3dInt_BNonLin(actSDMat);    
		nLinPart2 = new nLinMech3dInt_PiolaStress(actSDMat);
	      }
	    else if (subType_ == "plainStrain")
	      {
		nLinPart1 = new nLinMechPlaneStrainInt_BNonLin(actSDMat);    
		nLinPart2 = new nLinMechPlaneStrainInt_PiolaStress(actSDMat);

	      }
	    else if (subType_ == "axi")
	      {
		nLinPart1 = new nLinMechAxiInt_BNonLin(actSDMat);    
		nLinPart2 = new nLinMechAxiInt_PiolaStress(actSDMat);

	      }
	  
	    assemble_->AddIntegrator(nLinPart1, subdoms_[actSD], STIFFNESS,
				     nonLin_);
	    assemble_->AddIntegrator(nLinPart2, subdoms_[actSD], STIFFNESS,
				     nonLin_);

	    // assemble prestress, if in config-file given
	    // 	  if (preStressVal_)
	    // 	    AssemblePreStressMat(ptEl, connect_PDE, ptCoord,
	    //      actMatData, elDisp);
	  }


	// ==============  add mass ===========================================
	double density = actSDMat.GetDensity();
	BaseForm * bilinearMass  = new MassInt(density, dofspernode_, isaxi_);

	IntegratorDescriptor * actIntDescr =
	  new IntegratorDescriptor(bilinearMass, MASS);

	//check for damping (mass part)
	if (damping_type_ == RAYLEIGH)    
	  actIntDescr->SetSecondaryMat(DAMPING, actSDMat.GetDampingAlfa());

	assemble_->AddIntegrator(actIntDescr, subdoms_[actSD]);


	// ==================== RHS ===========================================
	if (nonLin_)
	  {
	    BaseForm * rhsSource = new nLinMech_linFormInt(actSDMat, isaxi_);
	    assemble_->AddRhsIntegrator(rhsSource, subdoms_[actSD], nonLin_);
	  }      
      }
  }


  BaseForm *
  MechPDE::GetStiffIntegrator(MaterialData& actSDMat, Boolean reducedInt)
  {

    ENTER_FCN( "MechPDE::GetStiffIntegrator" );
  
    BaseForm * bilinearStiff = NULL;

    if (subType_ == "plainStrain")
      bilinearStiff = new mechPlainStrainInt(actSDMat);
    else if (subType_ == "axi")
      bilinearStiff = new mechAxiInt(actSDMat);
    else if (subType_ == "3d")
      bilinearStiff = new mech3DInt(actSDMat);
    else 
      Info->Error("Unknown subtype in mech PDE! ",__FILE__,__LINE__);

    return bilinearStiff;
  }



  Integer MechPDE::GetBCDof(const std::string dofString)
  {
    ENTER_FCN( "MechPDE::GetBCDof" );

    Integer retVal = 0;

    if (dofString == "ux") retVal = 1;
    else if (dofString == "uy") retVal = 2;
    else if (dofString == "uz") retVal = 3;
    else
      {
#ifndef XMLPARAMS
	std::string errmsg = "The direction '" + dofString;
	errmsg += "' mentioned in the config-file is not implemented";
	Info->Error( errmsg, __FILE__, __LINE__ );
#else
	// According to the Schema definition of the parameter file this cannot
	// happen. Did the parser not perform validation?
	std::string errmsg = "Direction should be one of ux, uy, uz\n";
	  errmsg += "and not" + dofString;
	Info->Error( errmsg, __FILE__, __LINE__ );
#endif
      }

    return retVal;
  }
  


  void MechPDE::
  GetDirection(Directions& dir, const std::string keyword)
  {
    ENTER_FCN( "MechPDE::GetDirection" );

    std::string dirString;

#ifndef XMLPARAMS
    conf->getstr(keyword, dirString, pdename_);
#else
    params->Get( keyword, dirString, pdename_);
#endif

    if (dirString == "x")
      dir = X;
    else if (dirString == "y")
      dir = Y;
    else if (dirString == "z")
      dir = Z;
    else if (dirString == "radXY")
      dir = radXY;
    else if (dirString == "radXZ")
      dir = radXZ;
    else if (dirString == "radYZ")
      dir = radYZ;
    else
      {
#ifndef XMLPARAMS
	std::string errmsg = "The direction '" + dirString;
	errmsg += "' mentioned in the config-file is not implemented";
	Info->Error( errmsg, __FILE__, __LINE__ );
#else
	// According to the Schema definition of the parameter file this cannot
	// happen. Did the parser not perform validation?
	std::string errmsg = "Direction should be one of x, y, z, radXY, ";
	  errmsg += "radXZ, radYZ\nand not" + dirString;
	Info->Error( errmsg, __FILE__, __LINE__ );
#endif
      }
  }




void MechPDE::CalcReducedMat(MaterialData& lambdaMat, MaterialData& mueMat, MaterialData& mat)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::CalcReducedMat" << std::endl;
#endif

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
  
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      if (ptCoupling_->GetOutputQuantity(i) == "mechdisplacement")
	{
	  // Intialize the memory of the coupling values
	  ptCoupling_->CreateStoreSol(i,MECH_DISPLACEMENT,isComplex_);
	}

      if (ptCoupling_->GetOutputQuantity(i) == "mechforce")
	{
	  // Intialize the memory of the coupling values
	  ptCoupling_->CreateStoreSol(i,MECH_FORCE,isComplex_); 
	}
    }

  iterCoupledCounter_ = 0;
}





void MechPDE::CalcOutputCoupling()
{
#ifdef TRACE
  (*trace) << "entering MechPDE::CalcOutputCoupling" << std::endl;
#endif  

  Integer dof;
  std::string quantity;
  std::vector<Integer> * couplingnodes = NULL;
  std::vector<Elem*> * couplingElems = NULL;
  std::vector<Elem*> * neighbours = NULL;
  std::vector<MaterialData*> * couplingMaterials = NULL;
  BaseStoreSol * values = NULL;
  
  TRY_CAST
  

  // loop over all output coupling quantities
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
    {
      quantity = ptCoupling_->GetOutputQuantity(i);
      ptCoupling_->GetOutputValues(i, values);

      PTRCAST(values,StoreSol<Double>,temp);
	
      switch(ptCoupling_->GetOutputType(i))
	{
	case NODE:
	  
	  if (quantity == "mechdisplacement")
	    {
	      ptCoupling_->GetOutputNodes(i, couplingnodes);
	    	      
	      sol_->NodeSolutionToCoupling(*values, *couplingnodes,Mesh2PDENode_);
	    }
	  

	  if (quantity == "mechforce")
	    {
	      ptCoupling_->GetOutputNodes(i, couplingnodes);
	      ptCoupling_->GetOutputElements(i, couplingElems);

	      ptCoupling_->GetOwnMaterials(i, couplingMaterials);
	      ptCoupling_->GetOutputNeighbourElems(i, neighbours);
	      dof = ptCoupling_->GetOutputDof(i);

	      if (!neighbours->size())
		{
		  std::string errMsg = "In mechanic PDE: No neighbour elements ";
		  errMsg += "for acoustic-coupling at output interface ";
		  errMsg += ptCoupling_->GetOutputRegion(i);
		  Error(errMsg.c_str(),  __FILE__,__LINE__);  
		}
	  
	      
	      CalcAcousticCouplingRHS(couplingElems, *couplingnodes, 
				      couplingMaterials, *temp, dof, neighbours);	      
	    } 
	  break;

	case ELEM:
	  Error("No Element coupling output", __FILE__,__LINE__);
	}
    }

  CATCH_CAST
}


void MechPDE::CalcAcousticCouplingRHS(std::vector<Elem*> * couplingElems, 
				      std::vector<Integer>& couplingNodes,
				      std::vector<MaterialData*>* couplingMaterials,
				      StoreSol<Double> & elemCouplingSols,
				      Integer couplingdofM,
				      std::vector<Elem*> * neighbours)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::CalcAcousticCouplingRHS" << std::endl;
#endif

  Double density ;  
  Integer nrNodesperEl;

  elemCouplingSols.Init(0.0);
  
  for (Integer actElem=0; actElem<couplingElems->size(); actElem++)
    {
      BaseFE * ptElem = (*couplingElems)[actElem]->ptElem;
      Vector<Integer> connecth = (*couplingElems)[actElem]->connect;
      
      Matrix<Double> ptCoord; 
      GetElemCoords(connecth, ptCoord, actlevel_);
      
      // get correct density belonging to the the neighbouring element
      // in the fluid subdomain
      density = (*couplingMaterials)[actElem]->GetDensity();
      
      BaseForm * bilinear_mass = new MassInt(ptElem, density, isaxi_);
      Matrix<Double> elemmat;
      bilinear_mass->CalcElementMatrix(ptCoord, elemmat);
      delete bilinear_mass;	  


      Vector<Integer> connect_PDE;
      Mesh2PDENode(connect_PDE, connecth, Mesh2PDENode_);
      
      Vector<Double> sol;
      GetDerivSolVecOfElement(sol, connect_PDE);	 

      Vector<Double> nSol(connecth.GetSize());   // solution in normal direction
      nSol.Init();
      

      // the normal vector points outwards of the mechanical domain
      // (see. Kaltenbacher, "Num. Sim. of Mech. Act. & Sens." chapter 8.2)
      Vector<Double> n;
      CalcLineNormalVec(n, *(*couplingElems)[actElem], *(*neighbours)[actElem]);


      for (Integer actNode=0; actNode < connecth.GetSize(); actNode++)
	for (Integer actDof=0; actDof<dofspernode_; actDof++)
	  nSol[actNode] += sol[actDof + actNode*dofspernode_] * n[actDof];


      Vector<Double> forceOnElem = elemmat * nSol;  
      
      for (Integer actNode=0; actNode<ptCoord.GetSizeRow(); actNode++)
	{
	  Integer nodePos = 0;
	  
	  while(connecth[actNode] != couplingNodes[nodePos] && nodePos < couplingNodes.size()) 
	    nodePos++;
	  elemCouplingSols(nodePos,0) += forceOnElem[actNode];
	}      
    }
} 





void MechPDE::GetSolOfElement( Matrix<Double>& sol, Vector<Integer>& connect_PDE)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::GetSolOfElement" << std::endl;
#endif

  // displacement of element nodes
  sol.Resize(dofspernode_, connect_PDE.GetSize());
  
  for(Integer actNode=0; actNode<connect_PDE.GetSize(); actNode++)
    for(Integer actDof=0; actDof < dofspernode_; actDof++)
      // sol[actDof][actNode] = sol_[actDof][connect_PDE[actNode]-1];
      sol_->Get(connect_PDE[actNode]-1,actDof,sol[actDof][actNode]);
}







Boolean MechPDE::HasOutput(std::string output)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::HasOutput" << std::endl;
#endif

  if (output == "mechdisplacement" || output == "mechforce")
    return TRUE;

  return FALSE;
}





// ======================================================
// STATIC SOLVING SECTION
// ======================================================

void MechPDE:: PreStepStatic(const Integer kstep, const Double asteptime,
			     const Integer level, const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering MechPDE:: PreStepStatic" << std::endl;
#endif

  if (PDEisCoupled_)
    // init RHS at this place, because forces of other PDEs are added to RHS afterwards
    algsys_->InitRHS();     
  else
    // if PDE is coupled, the solution of the prior outer loops must be kept
    algsys_->InitSol();

}


void MechPDE::StepStaticNonLin(const Integer kstep, const Double aTime,
			       const Integer level, const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::SolveStepStaticNonLin" << std::endl;
#endif

  const Integer job = 1;
  Boolean performOneMoreStep;
  Integer iterationCounter=0;
  
  Vector<Double> actSol(numPDENodes_ * dofspernode_);
  
  Integer k=0;
  for (Integer i=0; i< numPDENodes_; i++)
    for (Integer j=0; j<dofspernode_; j++)
      {
	//actSol[k] = sol_[j][i];
	sol_->Get(i,j,actSol[k]);
	k++;
      }


  Vector<Double> solIncrement(numPDENodes_ * dofspernode_);

  SetBCs(level, updateBCs_, 0);

  // setup right hand side ==============================================      

  Double extForcesL2Norm = SetExternalForces(level);  
  assemble_->AssembleNLRHS(level);


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
      assemble_->InitNonLinMatrices();
      assemble_->AssembleMatrices(level);
      
#ifdef USE_OLAS
      algsys_->BuildInDirichlet();
      algsys_->SetupPrecond(job);
#else
      algsys_->CalcPrecond(job);
#endif

      algsys_->Solve();

      

      // new solution is only an increment of the full solution =============
      StoreAlgsysToVec(solIncrement, algsys_->GetSolutionVal() );

      Double residualL2Norm;
      Double etaLineSearch=0;

      if (!lineSearch_)
	actSol += solIncrement;
      else
	// TRUE is for transient simulation
	residualL2Norm = LineSearch(solIncrement, actSol, etaLineSearch, level);
      
      //StoreVecToSolArray(actSol);
      sol_->SetCompleteVector(actSol);


      // recalculate RHS with new values to get new residual (f^(k+1))========
#ifndef USE_OLAS    
      std::vector<Double>  help;
      extForces_.ToStdVector(help);
      algsys_->InitRHS(help);
#endif

      assemble_->AssembleNLRHS(level);  // inner forces due to nonlin formulation


      // =====================================================================
      // calculation of error norms
      // =====================================================================
      if (!lineSearch_)
	{
	  Vector<Double> actRHS;
	  StoreAlgsysToVec(actRHS, algsys_->GetRHSVal() );       
	  
	  // calculation of residual error =======================================
	  residualL2Norm = RhsL2Norm(actRHS); // L2Norm of  ( f_i^(k+1) - f_a )
	}


      // calculation of residual error =======================================
      Double residualErr = residualL2Norm / extForcesL2Norm;


      // calculate incremental error ========================================
      Double solIncrL2Norm = solIncrement.NormL2();
      Double actSolL2Norm  = actSol.NormL2();
      
      Double incrementalErr;
      
      if (actSolL2Norm)
	incrementalErr = solIncrL2Norm / actSolL2Norm;
      else
	{
	  incrementalErr = solIncrL2Norm;
	  Warning("Zero solution vector!! ", __FILE__,__LINE__);      
	}
      


      // =====================================================================
      // output of norms and data
      // =====================================================================


      WriteClaNlNorms(iterationCounter, residualL2Norm, extForcesL2Norm, residualErr,  
		      solIncrL2Norm, actSolL2Norm, incrementalErr);

      
      Info->WriteNonLinIter(pdename_, iterationCounter, residualErr, incrementalErr, etaLineSearch);


      // boolean variable, holds condition if another iteration step is necessary
      performOneMoreStep = 
	(incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);      
      
      
      if (!(performOneMoreStep && iterationCounter < maxnumiter_))
	mycout << "incrementalErr " << incrementalErr << myendl 
	       << "incStopCrit_ " << incStopCrit_ << myendl
	       << "residualErr " << residualErr  << myendl
	       << "residualStopCrit_ " << residualStopCrit_ << myendl;

    }while(performOneMoreStep && iterationCounter < maxnumiter_);  

}




Double MechPDE::LineSearch(Vector<Double>& solIncrement, Vector<Double>& actSol, 
			   Double& etaLineSearch, Integer level, Boolean trans)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::LineSearch" << std::endl;
#endif

  Vector<Double> solOld(actSol);
  const Integer nrEtas = 4;
  const Double eta[nrEtas] = {1, 0.5, 0.25, 0.125};
  Double etaOpt;
  Double residualL2NormOpt = 1e15;
  
  for(Integer i=0; i<nrEtas; i++)
    {
      actSol = solIncrement * eta[i];
      actSol += solOld;

      sol_->SetCompleteVector(actSol);
      //StoreVecToSolArray(actSol);

      // recalculate RHS with new values to get new residual (f^(k+1))========
#ifndef USE_OLAS      
      std::vector<Double>  help;
      extForces_.ToStdVector(help);
      algsys_->InitRHS(help);
#endif

      if(trans)
	{
	  assemble_->AssembleNLRHS(level, lasttimecalc_);
	  TS_alg_->UpdateRHS(actSol);
	}
      else
	assemble_->AssembleNLRHS(level);



      // =====================================================================
      // calculation of error norms
      // =====================================================================
      Vector<Double> actRHS;
      StoreAlgsysToVec(actRHS, algsys_->GetRHSVal() );

      // calculation of residual error =======================================
      Double residualL2Norm = RhsL2Norm(actRHS); // L2Norm of  ( f_i^(k+1) - f_a )

      if (residualL2Norm < residualL2NormOpt)
	 {
	   residualL2NormOpt = residualL2Norm;
	   etaOpt = eta[i];
	 }
    }

  etaLineSearch = etaOpt;
  
  myCout << "EtaOpt = " << etaOpt << myEndl;
  
  actSol =  solIncrement * etaOpt;
  actSol += solOld;

  return residualL2NormOpt;
}



void MechPDE::WriteClaNlNorms(const Integer iterationCounter, const Double residualL2Norm, const Double extForcesL2Norm,
			      const Double residualErr, const Double solIncrL2Norm, const Double actSolL2Norm, 
			      const Double incrementalErr)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::WriteClaNlNorms" << std::endl;
#endif
  
      *cla << std::endl << " ======================================================= " << std::endl
		<< " NONLINEAR ITERATION " << iterationCounter << std::endl
		<< " ======================================================= " << std::endl;
      *cla << " === Residual norm:          " << residualL2Norm << std::endl;
      *cla << "     Norm of ext. forces:    " << extForcesL2Norm << std::endl;
      *cla << "     Residual error          " << residualErr << std::endl;

      *cla << " === Incremental sol L2Norm: " << solIncrL2Norm << std::endl;
      *cla << "     Actual solution L2Norm: " << actSolL2Norm << std::endl;
      *cla << "     Incremental error       " << incrementalErr << std::endl;      
}





void MechPDE :: PostStepStatic(const Integer kstep, const Double asteptime,
			    const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::PostStepStatic" << std::endl;
#endif

  if (PDEisCoupled_)
    iterCoupledCounter_++;

}


// ======================================================
// TRANSIENT SOLVING SECTION
// ======================================================


void MechPDE :: InitTimeStepping(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::InitTimeStepping" << std::endl;
#endif

  if (effectiveMass_)  
    TS_alg_ = new NewmarkEffMass(pdename_, algsys_, 1, numPDENodes_*dofspernode_, damping_type_);
  else
    TS_alg_ = new Newmark(pdename_, algsys_, 1, numPDENodes_*dofspernode_, damping_type_);

  TS_alg_->Init(matrix_factor_, dt);

}






void MechPDE::StepTransNonLin(const Integer kstep, const Double asteptime,
			    const Integer level, const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::StepTransNonLin" << std::endl;
#endif

  const Integer job = 1;
  const Integer update = 0;  
  
  static Integer timeStepCounter=1;
  Double * ptsol;
  Boolean performOneMoreStep;
  Integer iterationCounter=0;

  Vector<Double> actSol(numPDENodes_ * dofspernode_);
  Vector<Double> solIncrement(numPDENodes_ * dofspernode_);
  
  algsys_->InitRHS();

  // Cast BaseStoreSol into StoreSol<Double>,
  // since this function is only called
  // in the transient case
 TRY_CAST
  PTRCAST(sol_,StoreSol<Double>,solhelp)
  

  actSol = solhelp->GetCompleteVector();
  TS_alg_->Predictor(solhelp->GetCompleteVector());

  Double extForcesL2Norm = SetExternalForces(level);

  timeStepCounter++;

  //perform predictor step
  TS_alg_->UpdateRHS(actSol);
  SetBCs(level, update, lasttimecalc_);

  
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
      assemble_->InitNonLinMatrices();
      assemble_->AssembleMatrices(level);
      algsys_->ConstructEffectiveMatrix(matrix_factor_);

      TS_alg_->UpdateRHS(actSol);
      SetBCs(level, update, lasttimecalc_);

#ifdef USE_OLAS
      algsys_->BuildInDirichlet();
      algsys_->SetupPrecond(job);
#else
      algsys_->CalcPrecond(job);
#endif

      algsys_->Solve();

      // new solution is only an increment of the full solution =============
      StoreAlgsysToVec(solIncrement, algsys_->GetSolutionVal() );

      Double residualL2Norm;
      Double etaLineSearch = 0;
      
      if (!lineSearch_)
	actSol += solIncrement;
      else
	// TRUE is for transient simulation
	residualL2Norm = LineSearch(solIncrement, actSol, etaLineSearch, level, TRUE);


      sol_->SetCompleteVector(actSol);
      //StoreVecToSolArray(actSol);


      // recalculate RHS with new values to get new residual (f^(k+1))========
      algsys_->InitRHS();
      assemble_->AssembleSrcRHS(level,lasttimecalc_); 
      assemble_->AssembleNLRHS(level, lasttimecalc_);  // inner forces due to nonlin formulation
      TS_alg_->UpdateRHS(actSol);




      // =====================================================================
      // calculation of error norms
      // =====================================================================

      if (!lineSearch_)
	{
	  Vector<Double> actRHS;
	  StoreAlgsysToVec(actRHS, algsys_->GetRHSVal() );       
	  
	  // calculation of residual error =======================================
	  residualL2Norm = RhsL2Norm(actRHS); // L2Norm of  ( f_i^(k+1) - f_a )
	}
      
      Double residualErr    = residualL2Norm / extForcesL2Norm;
      
      if (extForcesL2Norm < NORM_EPS)
	residualErr = residualL2Norm;  // take the absolute error instead of the relative


      // calculate incremental error ========================================
      Double solIncrL2Norm = solIncrement.NormL2();
      Double actSolL2Norm = actSol.NormL2();
      Double incrementalErr;

      if (actSolL2Norm)
	incrementalErr = solIncrL2Norm / actSolL2Norm;
      else
	{
	  Warning("Zero solution vector!! ", __FILE__,__LINE__);
	  incrementalErr = -EPS; // don't block the iteration loop
	}
      
      if (actSolL2Norm < NORM_EPS)
	incrementalErr = actSolL2Norm; // take the absolute error instead of the relative
       



      // =====================================================================
      // output of norms and data
      // =====================================================================


      WriteClaNlNorms(iterationCounter, residualL2Norm, extForcesL2Norm, residualErr,  
		      solIncrL2Norm, actSolL2Norm, incrementalErr);
      
      Info->WriteNonLinIter(pdename_, iterationCounter, residualErr, incrementalErr, etaLineSearch);
      

      // boolean variable, holds condition if another iteration step is necessary
      performOneMoreStep = 
	(incrementalErr > incStopCrit_)||(residualErr > residualStopCrit_);

    } while(performOneMoreStep && iterationCounter < maxnumiter_);  

  
    //perform corrector step  
  TS_alg_->Corrector(solhelp->GetCompleteVector());
  CATCH_CAST
}







// sets external forces and returns L2Norm of them
Double MechPDE::SetExternalForces(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::SetExternalForces" << std::endl;
#endif

  Double extForcesL2Norm;  

  // account for bcs before first solving step =======================
  SetBCs(level, updateBCs_, 0);

  // to incorporate loads
  assemble_->AssembleSrcRHS(level, lasttimecalc_);
  

  // stores rhs vector into extForces and returns that L2-norm
  StoreAlgsysToVec(extForces_, algsys_->GetRHSVal() );


  extForcesL2Norm = extForces_.NormL2();
 
  //  if extForcesL2Norm is 0, no residual norm can be calculated
  if (!extForcesL2Norm)
    Warning("Zero external force vector!! ", __FILE__,__LINE__);

  
  
  return extForcesL2Norm;
}





// calculates L2-norm of RHS regarding dirichlet entries due to penalty formulation by setting them 0
Double MechPDE::RhsL2Norm(Vector<Double>& actRHS)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::RhsL2Norm" << std::endl;
#endif  
  Integer node, dof;
  
  std::list<Integer> nodes;

  //  myCout << " RHS " << actRHS << myEndl;
  
  
  // Eliminate dirichlet node from RHS (due to penalty formulation)
  for (Integer i=0; i< bcs_hd_.size(); i++)
    {
      dof = GetNrBCDof ( homDirichDof_[i] );
      nodes=ptBCs_->GetNodesLevel(bcs_hd_[i]);
      
      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++)
	{
	  node=*p;

	  actRHS[(Mesh2PDENode_[node-1]-1)*dofspernode_ + dof-1] = 0;
	}
    }

  return actRHS.NormL2();
}




// stores an algsys_ vector into a std::vector and returns that L2-norm
void MechPDE::StoreAlgsysToVec(Vector<Double>& vec, Double * pt)
{
#ifdef TRACE
  (*trace) << "entering MechPDE::StoreAlgsysToVec" << std::endl;
#endif  

  const Integer numElems = numPDENodes_ * dofspernode_;
  
  vec.Resize(numElems);

  for (Integer i=0; i<numElems; i++)   
    vec[i] = pt[i];
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

  StoreSol<Double> sol_mesh, solder1_mesh, solder2_mesh;
  StoreSol<Double> sol_der1Array, sol_der2Array;
  
  if (savesol_ == TRUE && (analysistype_== STATIC || analysistype_== TRANSIENT))
    {
      sol_->TransformNodeSolution(sol_mesh,PDE2MeshNode_,ptgrid_,actlevel_);
      //TransformNodeSolution(sol_mesh, sol_, PDE2MeshNode_);
      OutFile_->WriteNodeSolution(sol_mesh, laststepcalc_, lasttimecalc_,"displacement");
    }

  if (analysistype_== TRANSIENT)
    {
      if (savederiv1_ == TRUE)

	{
	  sol_der1Array.SetNumSolutions(1);
	  sol_der1Array.SetNumNodes(numPDENodes_);
	  sol_der1Array.SetSolutionType(MECH_VELOCITY);
	  sol_der1Array.SetDof(Dim_);
	  sol_der1Array.Init(0);
	  sol_der1Array.SetSolVector(MECH_VELOCITY,getS1());
	  
	  sol_der1Array.TransformNodeSolution(solder1_mesh,PDE2MeshNode_,ptgrid_,actlevel_);
	  OutFile_->WriteNodeSolution(solder1_mesh,laststepcalc_,lasttimecalc_,"velocity");
	}

      if (savederiv2_ == TRUE)
	{
	  sol_der2Array.SetNumSolutions(1);
	  sol_der2Array.SetNumNodes(numPDENodes_);
	  sol_der2Array.SetSolutionType(MECH_ACCELERATION);
	  sol_der2Array.SetDof(Dim_);
	  sol_der2Array.Init(0);
	  sol_der2Array.SetSolVector(MECH_ACCELERATION,getS2());
	  sol_der2Array.TransformNodeSolution(solder2_mesh,PDE2MeshNode_,ptgrid_,actlevel_);
	  OutFile_->WriteNodeSolution(solder2_mesh,laststepcalc_,lasttimecalc_,"acceleration");
	}
    }
}


} // end namespace CoupledField

