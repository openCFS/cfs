#include <stdio.h>

#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "Forms/forms_header.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "Driver/assemble.hh"
#include "newmark.hh"
#include "Elements/basefe.hh"
#include "blocknodeEQN.hh"
#include "superblockEQN.hh"
#include "scalarblockEQN.hh"

#include "piezoPDE.hh" 

namespace CoupledField {

  PiezoPDE::PiezoPDE( Grid *aptgrid, BCs *aptbcs,  TimeFunc *aptTimeFunc,
		      FileType *aptFileType, WriteResults *aptOut ) :
    BasePDE( aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc ) {

    ENTER_FCN( "PiezoPDE::PiezoPDE" );

    // Set name of the PDE and its material class
    pdename_ = "piezo";
    pdematerialclass_ = "piezo";

#ifndef XMLPARAMS
    // Determine dimension of problem
    conf->getstr( "subtype", subType_, pdename_ ); 
    if (subType_ == "3d") {
      dofspernode_ = 4;
      Info->PrintF("", "=== 3D PROBLEM\n");
    }
    else if (subType_ == "axi") {
      isaxi_ = TRUE;
      dofspernode_ = 3;
      Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
    }
    else {
      // default is planeStrain 
      dofspernode_ = 3;
      Info->PrintF("", "=== PLAIN STRAIN PROBLEM\n");
    }

    conf->getsubdompde(subdoms_,pdename_);
#else

    // Get problem geometry and PDE subtype
    params->Get( "subtype", subType_, pdename_ );
    std::string probGeo;
    params->Get( "type", probGeo, "geometry" );

    // Set number of degrees of freedom and
    // ensure that subtype fits to problem geometry
    if ( subType_ == "3d" && probGeo == "3d" ) {
      dofspernode_ = 4;
      Info->PrintF("", "=== 3D PROBLEM\n");
    }
    else if ( subType_ == "axi" && probGeo == "axi" ) {
      isaxi_ = TRUE;
      dofspernode_ = 3;
      Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
    }
    else if ( subType_ == "plainStrain" && probGeo == "plane" ) {
	dofspernode_ = 3;
	Info->PrintF("", "=== PLAIN STRAIN PROBLEM\n");
    }
    else {
      std::string errmsg = "Subtype '" + subType_;
      errmsg += "' of PDE '" + pdename_ + "' does not fit to problem ";
      errmsg += "geometry '" + probGeo + "'\n";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // Get list of subdomains
    params->GetList( "name", subdoms_, pdename_, "region" );
    Info->PrintF( pdename_, " PiezoPDE lives on regions:" );
    for ( Integer k = 0; k < subdoms_.GetSize(); k++ ) {
      Info->PrintF( pdename_, " %s", subdoms_[k].c_str() );
    }

#endif

    ReadBCs(pdename_);
    
   

#ifndef XMLPARAMS
   effectiveMass_ = FALSE;
   if (conf->get_option("effMass",  pdename_ ))
     effectiveMass_ = TRUE;

   //check for damping model
   std::string dampstr;
   conf->ifget("damping",dampstr,pdename_);
   if (dampstr == "rayleigh")
     {
       //       Error("Currenrly Rayleigh damping not woirking for PiezoPDE",__FILE__,__LINE__);       
       dampingType_ = RAYLEIGH;
     }
   else
     dampingType_ = NONE;
   
   if (dampingType_)
     assemble_->NeedDampingMatrix();
#else

   // Use effective mass approach?
   effectiveMass_ = params->IsSet( "effMass" );

   // Do we use damping?
   if( params->HasValue( "type", "rayleigh", pdename_, "damping" ) ) {
     dampingType_ = RAYLEIGH;
     Info->PrintF( pdename_, " Using RAYLEIGH damping\n" );
   }
   else {
     dampingType_ = NONE;
     Info->PrintF( pdename_, " Using no damping\n" );
   }
   if( dampingType_ != NONE ) {
     assemble_->NeedDampingMatrix();
   }

#endif

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
   if (bcs_hd_.GetSize() != homDirichDof_.GetSize()) 
     Error("Inconsistent definition of homogeneous Dirichlet Boundary Conditions");
   if (bcs_id_.GetSize() != inhomDirichDof_.GetSize()) 
     Error("Inconsistent definition of inhomogeneous Dirichlet Boundary Conditions");
   
   // Two possible numberings of equations possible
   // Simply unquote the desired line
   // 1.) Super-blocknumbering
   SuperBlockEQN dummy = SuperBlockEQN(ptgrid_, ptBCs_, subdoms_, actlevel_, dofspernode_);
   //eqnData_ = new SuperBlockEQN(ptgrid_, ptBCs_, subdoms_, actlevel_, dofspernode_);
   
   // 2.) Normal blocknumbering
   //eqnData_ = new BlockNodeEQN(ptgrid_, ptBCs_, subdoms_, actlevel_, dofspernode_);
   
   // 3.) Scalar blocknumbering
   eqnData_ = new ScalarBlockEQN(ptgrid_, ptBCs_, subdoms_, actlevel_, dofspernode_);

   eqnData_->SetHomoDirichletBCs(bcs_hd_, homDirichDof_);
   eqnData_->CalcMapping();
   //eqnData_->Print(std::cerr);
   
   numPDENodes_ = eqnData_->GetNumLocalNodes();
   numElems_ = eqnData_->GetNumLocalElems();
   
   size_        = numPDENodes_ * dofspernode_;
   sol_->SetNumSolutions(2);
   sol_->SetNumNodes(numPDENodes_);
   sol_->SetSolutionType(MECH_DISPLACEMENT,0);
   sol_->SetSolutionType(ELEC_POTENTIAL,1);
   sol_->SetNumDofs(dim_,MECH_DISPLACEMENT); // displacements have dof of mesh-dimension
   sol_->SetNumDofs(1,ELEC_POTENTIAL);  // electric potential
   sol_->SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
   sol_->Init(0.0);

   // set assemble parameters
   assemble_->SetPtr2EQNData(eqnData_); 
   assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, subdoms_, surfdoms_);
   assemble_->SetGraphType(NODEGRAPH);

#ifdef USE_OLAS
   assemble_->SetMatrixEntryType(OLAS::DOUBLE);
   assemble_->SetMatrixStorageType(OLAS::SPARSE_NONSYM);
#else
   if (eqnData_->GetNumDofsPerEQN() == 1)
     assemble_->SetMatrixType(RSCALAR);
   else
     assemble_->SetMatrixType(RBLOCK);
#endif 
   
   assemble_->SetNumDirichlet(GetNumRestraints(actlevel_));
   
   assemble_->SetPtrBCs(ptBCs_);
   assemble_->SetPtr2Sol(sol_);
   assemble_->SetPtr2TimeFnc(ptTimeFunc_);
   
   ReadMaterialData();
   
   DefineIntegrators(actlevel_);  
   
#ifndef XMLPARAMS
   ReadSavings();
#else
   ReadStoreResults();
#endif

  }


  // *********************
  //   DefineIntegrators
  // *********************
  void PiezoPDE::DefineIntegrators( const Integer level ) {
    ENTER_FCN( "PiezoPDE::DefineIntegerators" );

    Boolean nonLin = FALSE;

    for ( int actSD = 0; actSD < subdoms_.GetSize(); actSD++ ) {

      // ==============  add stiffness ========================================

      MaterialData actSDMat(materialData_[actSD]);

      // ==============  add "standard" stiffness =============================

      BaseForm * bilinearStiff = GetStiffIntegrator(actSDMat);
      IntegratorDescriptor *actIntDescrStiff =
	new IntegratorDescriptor(bilinearStiff, STIFFNESS);
      assemble_->AddIntegrator(actIntDescrStiff, subdoms_[actSD]);
      
      //check for damping
      if ( dampingType_ == RAYLEIGH ) {
	Boolean isdamping = TRUE;
	Boolean reducedIntegration = FALSE; //is currently not supported
	BaseForm * dampStiff = GetStiffIntegrator( actSDMat,reducedIntegration,
						   isdamping );
	IntegratorDescriptor *actIntDescrDamp =
	  new IntegratorDescriptor(dampStiff, DAMPING);
	assemble_->AddIntegrator(actIntDescrDamp, subdoms_[actSD]);
      }

      // ==============  add mass =============================================
      Double density = actSDMat.GetDensity();
      Integer electricPot = 4;
      
      BaseForm * bilinearMass  = new MassInt(density, dofspernode_,
					     electricPot, isaxi_);

      IntegratorDescriptor * actIntDescrMass =
	new IntegratorDescriptor(bilinearMass, MASS);

      //check for damping (mass part)
      if (dampingType_ == RAYLEIGH)    
	actIntDescrMass->SetSecondaryMat(DAMPING, actSDMat.GetDampingAlfa(),
					 analysistype_);
      assemble_->AddIntegrator(actIntDescrMass, subdoms_[actSD]);
    }
  }


  // **********************
  //   GetStiffIntegrator
  // **********************
  BaseForm *PiezoPDE::GetStiffIntegrator( MaterialData& actSDMat,
					  Boolean reducedInt,
					  Boolean isdamping ) {

    ENTER_FCN( "PiezoPDE::GetStiffIntegrator" );
  
    BaseForm * bilinearStiff;
    
    if ( subType_ == "plainStrain" ) {
      Error( "plainStrain for PiezoPDE currently not supported",
	     __FILE__, __LINE__ );
      // bilinearStiff = new mechPlainStrainInt(actSDMat);
    }

    else if ( subType_ == "axi" ) {
      Error( "axi for PiezoPDE currently not supported", __FILE__, __LINE__ );
      // bilinearStiff = new mechAxiInt(actSDMat);
    }

    else if ( subType_ == "3d" ) {
      bilinearStiff = new linPiezo3DInt( actSDMat, isdamping );
    }

    else {
      Error( "Unknown subtype in piezo PDE!", __FILE__, __LINE__ );
    }

    return bilinearStiff;
  }


// ======================================================
// TRANSIENT SOLVING SECTION
// ======================================================


  void PiezoPDE::InitTimeStepping( const Double dt ) {
    ENTER_FCN( "PiezoPDE::InitTimeStepping" );

    if (effectiveMass_)
      TS_alg_ = new NewmarkEffMass(pdename_, algsys_, eqnData_, dampingType_);
    else
      TS_alg_ = new Newmark(pdename_, algsys_, eqnData_, dampingType_);

    TS_alg_->Init(matrix_factor_, dt);
  }


  void PiezoPDE::WriteResultsInFile()
  {
    ENTER_FCN( "PiezoPDE::WriteResultsInFile" );

    NodeStoreSol<Double> * solTransient;
    NodeStoreSol<Complex> * solHarmonic;

    if (analysistype_ == STATIC ||
	analysistype_ == HARMONIC) {
      
      if (savesol_) {
	solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);
	outFile_->WriteNodeSolutionTransient(*solTransient, laststepcalc_, lasttimecalc_);
      }
      
      //element results
      if (calcStress_.GetSize() !=0 ) {
	outFile_->WriteElemSolutionTransient(Stress_, laststepcalc_, lasttimecalc_);
      }
    } else if (analysistype_ == HARMONIC) {
      
      if (savesol_)
	{
	  solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);
	  outFile_->WriteNodeSolutionHarmonic(*solHarmonic,  actFreqStep_, 
					      actFrequency_, complexFormat_);
	}
    }
  }
    


  Integer PiezoPDE::GetBCDof(const std::string dofStartString)
  {
    ENTER_FCN( "PiezoPDE::GetBCDof" );

    Integer nrActDof = 0;
    
    if ( dofStartString == "ux" )
      nrActDof = 1;
    if ( dofStartString == "uy" )
      nrActDof = 2;
    if ( dofStartString == "uz" )
      nrActDof = 3;
    if ( dofStartString == "ep" )
      nrActDof = 4;
    if ( nrActDof == 0 )
      Error("Unknown dof-type in homog. BC; substring must start with ux, uy, uz or ep!!",__FILE__,__LINE__);

    return nrActDof;
  }


// ***********************************************************************
//   Obtain information on desired output quantities from parameter file
// ***********************************************************************
#ifdef XMLPARAMS
void PiezoPDE::ReadStoreResults() {

  ENTER_FCN( "PiezoPDE::ReadStoreResults" );

  // By default we only save the solution at nodal values
  savesol_    = TRUE;
  savederiv1_ = FALSE;
  savederiv2_ = FALSE;

  // Determine what solution values the user wants to be stored
  StdVector<std::string> nodeValues;
  params->GetList( "type", nodeValues, pdename_, "nodeHistory" );

  for ( Integer i = 0;  i < nodeValues.GetSize(); i++ ) {
    if ( nodeValues[i] == "displacement" ) savesol_    = TRUE;
    if ( nodeValues[i] == "velocity"     ) savederiv1_ = TRUE;
    if ( nodeValues[i] == "acceleration" ) savederiv2_ = TRUE;
  }

  // -----------
  //   Stress
  // -----------

  // Determine regions for which stress must be computed
  params->CGetList( "region", calcStress_, "type", "stress", 0, pdename_,
		      "elemResults" );

  // If the symbolic name is "all" compute electric field for all regions
  if ( calcStress_.GetSize() == 1 && calcStress_[0] == "all" ) {
    calcStress_ = subdoms_;
  }

  // Log to info file
  if ( calcStress_.GetSize() > 0 ) {
    Info->PrintF( pdename_,
		  " Computing mechanical stress for regions:");
    for ( Integer k = 0; k < calcStress_.GetSize(); k++ ) {
      Info->PrintF( pdename_, " %s", calcStress_[k].c_str() );
    }
  }
  

    //check for charge computation
    params->GetList( "region", chargeNeighborRegion_, pdename_, "charge" );
    params->GetList( "element", calcCharge_, pdename_, "charge" );

  }

#endif

// ************************************************************
//   PostProcess
// ************************************************************

void PiezoPDE::PostProcess(const Integer level) {

  ENTER_FCN( "PiezoPDE::PostProcess" );

  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
  
  if (calcStress_.GetSize() !=0 ) {

    //get the correct bilinearform
    ShortInt stressElecDim, stressDim, elecDim;
    Vector<Double> intPoint;

#ifndef XMLPARAMS 
      if (subType_ == "plainStrain") 
#else
      if (subType_ == "planeStrain") 
#endif
	{
	  
	  stressDim = 3;
	  elecDim   = 2;
	  intPoint.Resize(2); 
	  intPoint.Init(0);
	}

      else if (subType_ == "axi") {
	stressDim = 4;
	elecDim   = 2;
	intPoint.Resize(2); 
	intPoint.Init(0);
      }

      else if (subType_ == "3d") {
	stressDim = 6;
	elecDim   = 3;
	intPoint.Resize(3); 
	intPoint.Init(0);
      }

      else 
	Info->Error("StressOp: Unknown subtype in mech PDE! ",__FILE__,__LINE__);  
      
      
      // Resize solution arrays
      Stress_.SetNumSolutions(1);
      Stress_.SetSolutionType(MECH_STRESS);
      Stress_.SetNumElems(numElems_);
      Stress_.SetNumDofs(stressDim);
      Stress_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
      Stress_.Init(0);
      
      Vector<Double> elemElecStress, elemStress;
      elemElecStress.Resize(stressDim+elecDim);
      elemStress.Resize(stressDim);
      elemElecStress.Init(0);
      elemStress.Init(0);
      
      // loop over all subdomains
      for (Integer isd=0; isd<subdoms_.GetSize(); isd++) {
	
	MaterialData actSDMat(materialData_[isd]);
	PiezoStressStrain *stress;
	
// 	if (subType_ == "planeStrain") 
// 	  stress = new PiezoStressStrainPlaneStrain(actSDMat);

// 	else if (subType_ == "axi") 
// 	  stress = new PiezoStressStrainAxi(actSDMat);

	if (subType_ == "3d") 
	  stress = new PiezoStressStrain3D(actSDMat);
	else
	  Error("piezoPDE: Stress computation for 2D plane and axi not supported");

	
	// get vector of Elements of subdomains
	StdVector<Elem*> elemssd;     
	ptgrid_->GetElemSD(elemssd,subdoms_[isd],level);
	
	// loop over elements of subdomain
	for (Integer iel=0; iel< elemssd.GetSize(); iel++) {
	  Integer pdeElem = eqnData_->Mesh2PDEElem(elemssd[iel]->elemNum);
	  
	  //set element pointer
	  BaseFE * ptEl = elemssd[iel]->ptElem;
	  stress->SetElemPtr(ptEl);
	  
	  //set element solution	
	  Matrix<Double> elSol;
	  StdVector<Integer> connecth = elemssd[iel]->connect;
	  sol_->GetElemSolutionAsMatrix(elSol, connecth);
	  stress->SetActElemSol(elSol);
	  
	  //get coordinates of element
	  Matrix<Double> ptCoord;
	  GetElemCoords(connecth, ptCoord, level);
	  
	  Vector<Double> actStress;	
	  
	  //set the integration point
	  stress->SetIntPoint(intPoint);

	  //calculates the stress
	  stress->CalcStressVec(elemElecStress,1,ptCoord);
	  elemStress = elemElecStress.Part(0,stressDim-1);
	  Stress_.SetElemResult(pdeElem-1, elemStress);
	}
      }
  }

  //calc charges
  if (calcCharge_.GetSize() !=0 ) {
    for (Integer i=0; i<calcCharge_.GetSize(); i++)
      std::cout << "Calc charge for:" << calcCharge_[i] << std::endl;
  }
  

          
  }


} // end of namespace
