#include <stdio.h>

#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "Forms/forms_header.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "Driver/assemble.hh"
#include "newmark.hh"
#include "Elements/basefe.hh"

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
      dofspernode_ = 3;
      Info->PrintF("", "=== 3D PROBLEM\n");
    }
    else if ( subType_ == "axi" && probGeo == "axi" ) {
      isaxi_ = TRUE;
      dofspernode_ = 2;
      Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
    }
    else if ( subType_ == "plainStrain" && probGeo == "plane" ) {
	dofspernode_ = 2;
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
    for ( Integer k = 0; k < subdoms_.size(); k++ ) {
      Info->PrintF( pdename_, " %s", subdoms_[k].c_str() );
    }

#endif

    ReadBCs(pdename_);
   
   // Map global numeration of element and nodes to local one
   AssignPDENodeNumbers(mesh2PDENode_, pde2MeshNode_, subdoms_);  
   AssignPDEElemNumbers(mesh2PDEElem_, pde2MeshElem_, subdoms_);
   numPDENodes_ = pde2MeshNode_.size();
   numElems_ = pde2MeshElem_.size(); 
   
   size_        = numPDENodes_ * dofspernode_;

   sol_->SetNumSolutions(2);
   sol_->SetNumNodes(numPDENodes_);
   sol_->SetSolutionType(MECH_DISPLACEMENT,0);
   sol_->SetSolutionType(ELEC_POTENTIAL,1);
   sol_->SetNumDofs(dim_,MECH_DISPLACEMENT); // displacements have dof of mesh-dimension
   sol_->SetNumDofs(1,ELEC_POTENTIAL);  // electric potential
   sol_->Init(0.0);

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
     Info->PrintF(pdename_, " Using RAYLEIGH damping\n" );
   }
   else {
     dampingType_ = NONE;
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
   if (bcs_hd_.size() != homDirichDof_.size()) 
     Error("Inconsistent definition of homogeneous Dirichlet Boundary Conditions");
   if (bcs_id_.size() != inhomDirichDof_.size()) 
     Error("Inconsistent definition of inhomogeneous Dirichlet Boundary Conditions");
 
  // set assemble parameters
  assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, subdoms_, surfdoms_);
  assemble_->SetGraphType(NODEGRAPH);
  assemble_->SetMesh2PDENode(&mesh2PDENode_);
#ifdef USE_OLAS
  assemble_->SetMatrixEntryType(OLAS::DOUBLE);
  assemble_->SetMatrixStorageType(OLAS::SPARSE_NONSYM);
#else
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


  void PiezoPDE::DefineIntegrators(const Integer level)
  {
    ENTER_FCN( "PiezoPDE::DefineIntegerators" );

  Boolean nonLin = FALSE;


  for (int actSD = 0; actSD < subdoms_.size(); actSD++)
    {

      // ==============  add stiffness ===========================================

      MaterialData actSDMat(materialData_[actSD]);

      // ==============  add "standard" stiffness ===============================
      BaseForm * bilinearStiff = GetStiffIntegrator(actSDMat);
      IntegratorDescriptor * actIntDescrStiff = new IntegratorDescriptor(bilinearStiff, STIFFNESS);
      assemble_->AddIntegrator(actIntDescrStiff, subdoms_[actSD]);
      
      //check for damping
      if (dampingType_ == RAYLEIGH)    
	{
	  Boolean isdamping = TRUE;
	  Boolean reducedIntegration = FALSE; //is currently not supported
	  BaseForm * dampStiff = GetStiffIntegrator(actSDMat,reducedIntegration,
						    isdamping);
	  IntegratorDescriptor * actIntDescrDamp = new IntegratorDescriptor(dampStiff, DAMPING);
	  assemble_->AddIntegrator(actIntDescrDamp, subdoms_[actSD]);
	}

	  
      // ==============  add mass ================================================
      Double density = actSDMat.GetDensity();
      Integer electricPot = 4;
      
      BaseForm * bilinearMass  = new MassInt(density, dofspernode_, electricPot, isaxi_);

      IntegratorDescriptor * actIntDescrMass = new IntegratorDescriptor(bilinearMass, MASS);

      //check for damping (mass part)
      if (dampingType_ == RAYLEIGH)    
	actIntDescrMass->SetSecondaryMat(DAMPING, actSDMat.GetDampingAlfa(),analysistype_);

      assemble_->AddIntegrator(actIntDescrMass, subdoms_[actSD]);

    }
  }


  BaseForm *
  PiezoPDE::GetStiffIntegrator(MaterialData& actSDMat, Boolean reducedInt, Boolean isdamping)
  {
    ENTER_FCN( "PiezoPDE::GetStiffIntegrator" );
  
    BaseForm * bilinearStiff;
    
    if (subType_ == "plainStrain")
      Error("plainStrain for PiezoPDE currently not supported",__FILE__,__LINE__);
    //bilinearStiff = new mechPlainStrainInt(actSDMat);
    else if (subType_ == "axi")
      Error("axi for PiezoPDE currently not supported",__FILE__,__LINE__);
    //    bilinearStiff = new mechAxiInt(actSDMat);
    else if (subType_ == "3d")
      bilinearStiff = new linPiezo3DInt(actSDMat,isdamping);
    else 
      Error("Unknown subtype in mech PDE! ",__FILE__,__LINE__);
    
    return bilinearStiff;
  }


// ======================================================
// TRANSIENT SOLVING SECTION
// ======================================================


void PiezoPDE :: InitTimeStepping(const Double dt)
{
  ENTER_FCN( "PiezoPDE::InitTimeStepping" );

  if (effectiveMass_)  
    TS_alg_ = new NewmarkEffMass(pdename_, algsys_, 1, numPDENodes_*dofspernode_, dampingType_);
  else
    TS_alg_ = new Newmark(pdename_, algsys_, 1, numPDENodes_*dofspernode_, dampingType_);

  TS_alg_->Init(matrix_factor_, dt);

}


  void PiezoPDE::WriteResultsInFile()
  {
    ENTER_FCN( "PiezoPDE::WriteResultsInFile" );

    StoreSol<Double> DispMesh;
    StoreSol<Double> PotentialMesh;
    StoreSol<Double> DispPDE, PotentialPDE;

    sol_->GetSolution(MECH_DISPLACEMENT,DispPDE);
    sol_->GetSolution(ELEC_POTENTIAL,PotentialPDE);

    DispPDE.TransformNodeSolution( DispMesh,pde2MeshNode_,ptgrid_,actlevel_);
    PotentialPDE.TransformNodeSolution( PotentialMesh,pde2MeshNode_,ptgrid_,actlevel_);

    outFile_->WriteNodeSolution(DispMesh, laststepcalc_, lasttimecalc_,
				"displacement");
    outFile_->WriteNodeSolution(PotentialMesh, laststepcalc_, lasttimecalc_,
				"E-Potential");
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

} // end of namespace
