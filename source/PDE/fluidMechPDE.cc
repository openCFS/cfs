// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "fluidMechPDE.hh"

#include <sstream>
#include <iomanip>

#include "Forms/forms_header.hh"
#include "Forms/linElastInt.hh"
#include "Forms/massInt.hh"
#include "Forms/linPressureInt.hh"
#include "Forms/singleEntryInt.hh"
#include "Forms/linSurfStressInt.hh"
#include "Driver/assemble.hh"
#include "bdf2.hh"
#include "trapezoidal.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/resultHandler.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"
#include "Domain/ansatzFct.hh"
#include "Driver/solveStepFluidMech.hh"
#include "Utils/SmoothSpline.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/SIMP.hh"

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

namespace CoupledField {

DECLARE_LOG(fluidmechpde)
DEFINE_LOG(fluidmechpde, "fluidmechpde")


    FluidMechPDE::FluidMechPDE(Grid * aptgrid, ParamNode* paramNode )
      :SinglePDE( aptgrid, paramNode ) {

      pdename_          = "fluidMech";
      pdematerialclass_ = FLOW;
      maxTimeDerivOrder_ = 2;

      nonLin_      = false;

      needSolPrev_ = true;

      //set coordinate Update true to obtain the updated Lagrange formulation
      coordUp_ = true;
   
      gridSol_ = new NodeStoreSol<Double>;
      // ****************************
      // DETERMINE GEOMETRY
      // ****************************

      // Get problem geometry and PDE subtype
      myParam_->Get("subType", subType_ );

      std::string probGeo;
      param->Get("domain")->Get("geometryType", probGeo );

      // Set number of degrees of freedom and
      // ensure that subtype fits to problem geometry
      if ( subType_ == "3d" && probGeo == "3d" ) {
        stressDim_ = 6;
        Info->PrintF("", "=== 3D PROBLEM\n");
      }
      else if ( subType_ == "axi" && probGeo == "axi" ) {
        isaxi_ = true;
        stressDim_ = 4;
        Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
      }
      else if ( subType_ == "plane" && probGeo == "plane" ) {
        stressDim_ = 3;
        Info->PrintF("", "=== PLANE PROBLEM\n");
      }
      

      else {
        EXCEPTION( "Subtype '" <<  subType_ << "' of PDE '"
                   <<  pdename_ <<  "' does not fit to problem  geometry '"
                   << probGeo << "'"; );
      }

      approxType_ = myParam_->Get("type")->AsString();
      stabilizationType_ = myParam_->Get("stabilizationType")->AsString();
      movingMesh_ = myParam_->Get("movingMesh")->AsBool();

      forceFac_ = 1.0;
      
      // Register functions for scripting
      //RegisterFunctions();

    }

  FluidMechPDE::~FluidMechPDE()
  {

  }

  void FluidMechPDE::ReadSpecialBCs() {
    
    // read pressure information
    ReadPressureLoads();

  }

  void FluidMechPDE::InitNonLin()
  {
   // ==============================================================
   // NOTE: We assume that for a fluid mechanical PDE all regions
   //       either are non-linear due to the convective term!
   // ==============================================================
    
    nonLin_ = true;

  }
  

  void FluidMechPDE::DefineIntegrators() 
  {
    
    Double density, dynamicViscosity, kinematicViscosity;

    // help variables for parameter checking
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;
    
    //volume integrators

    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {

      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;

      actSDMat->GetScalar(density,DENSITY,REAL);
      actSDMat->GetScalar(dynamicViscosity,DYNAMIC_VISCOSITY,REAL);
      actSDMat->GetScalar(kinematicViscosity,KINEMATIC_VISCOSITY,REAL);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );

      // get current region name and get grip on paramNode
      std::string actRegionName;
      actRegionName = ptgrid_->RegionIdToName( actRegion );

      // ==============  add stiffness matrix ===============

      //****************************************************************
      //******** stabilized ********************************************
      //****************************************************************
      BaseForm * bilinearStiff_UV = GetStiffIntegrator_UV(density,kinematicViscosity);
      BiLinFormContext * actStiffContext_UV = new BiLinFormContext(bilinearStiff_UV, STIFFNESS );
      actStiffContext_UV->SetPtPdes(this, this);
      actStiffContext_UV->SetResults( results_[0],results_[0],actSDList, actSDList );
      assemble_->AddBiLinearForm( actStiffContext_UV );
      //******************************************************************************************************************
      if( stabilizationType_ != "none"){	
    	  BaseForm * bilinearStiff_PQ = GetStiffIntegrator_PQ(density,kinematicViscosity);
    	  BiLinFormContext * actStiffContext_PQ = new BiLinFormContext(bilinearStiff_PQ, STIFFNESS );
    	  actStiffContext_PQ->SetPtPdes(this, this);
    	  actStiffContext_PQ->SetResults( results_[1],results_[1],actSDList, actSDList );
    	  assemble_->AddBiLinearForm( actStiffContext_PQ );
      }

      //******************************************************************************************************************
      if (approxType_!= "taylorHood"){
    	  BaseForm * bilinearStiff_UQ = GetStiffIntegrator_UQ(density,kinematicViscosity);
    	  BiLinFormContext * actStiffContext_UQ = new BiLinFormContext(bilinearStiff_UQ, STIFFNESS );
    	  actStiffContext_UQ->SetPtPdes(this, this);
    	  actStiffContext_UQ->SetResults( results_[1],results_[0],actSDList, actSDList );
    	  assemble_->AddBiLinearForm( actStiffContext_UQ );
      }
      //******************************************************************************************************************
      BaseForm * bilinearStiff_PV = GetStiffIntegrator_PV(density,kinematicViscosity);
      BiLinFormContext * actStiffContext_PV = new BiLinFormContext(bilinearStiff_PV, STIFFNESS );
      // We also need to set the transposed of the coupling
      // matrix to the lower diagonal side
      if (approxType_== "taylorHood")
    	  actStiffContext_PV->SetCounterPart( true );
      actStiffContext_PV->SetPtPdes(this, this);
      actStiffContext_PV->SetResults( results_[0],results_[1],actSDList, actSDList );
      assemble_->AddBiLinearForm( actStiffContext_PV );
   	  
          // ==============  add mass matrix ====================
      if (analysistype_ == TRANSIENT) {
    	  isInstationary_=true;
    	  //****************************************************************
    	  //******** stabilized ********************************************
    	  //****************************************************************
    	  //****************************************************************
    	  BaseForm * bilinearMass_UV = GetMassIntegrator_UV(density,kinematicViscosity);
    	  BiLinFormContext * actMassContext_UV = new BiLinFormContext(bilinearMass_UV, MASS );
    	  actMassContext_UV->SetPtPdes(this, this);
    	  actMassContext_UV->SetResults(results_[0],results_[0],actSDList, actSDList );
    	  assemble_->AddBiLinearForm( actMassContext_UV );

    	  if(stabilizationType_ != "none"){	
    		  //****************************************************************
    		  BaseForm * bilinearMass_UQ = GetMassIntegrator_UQ(density,kinematicViscosity);
    		  BiLinFormContext * actMassContext_UQ = new BiLinFormContext(bilinearMass_UQ, MASS );
    		  actMassContext_UQ->SetPtPdes(this, this);
    		  actMassContext_UQ->SetResults(results_[1],results_[0],actSDList,actSDList);
    		  assemble_->AddBiLinearForm( actMassContext_UQ );
    		  //****************************************************************
    	  }
      }
      eqnMap_->AddResult( *results_[0], actSDList );
      eqnMap_->AddResult( *results_[1], actSDList );
    }
    
    //surface integrators
    //RHS-part
    for (UInt actSF = 0; actSF < pressSurf_.GetSize(); actSF++) {

    	LinearSurfForm * rhsSrcSurf = 
    		new PressureLinForm(pressVals_[actSF], pressPhase_[actSF], 
    				subType_, isaxi_ );
    	rhsSrcSurf->SetVoluInfo( materials_ );

    	LinearFormContext * pressRhs = 
    		new LinearFormContext( rhsSrcSurf );
    	pressRhs->SetPtPde( this );
    	pressRhs->SetResult( results_[0], pressSurf_[actSF] );
    	assemble_->AddLinearForm( pressRhs );

    	// Give entities and result to equation numbering class
    	// and solution class
    	eqnMap_->AddResult( *results_[0], pressSurf_[actSF] );
    }
  }


    BaseForm * FluidMechPDE::GetStiffIntegrator(BaseMaterial* actSDMat, 
                                       RegionIdType regionId,
                                       bool reducedInt) {
        // Get region name
        std::string regionName = ptgrid_->RegionIdToName( regionId );
    
        BaseForm * bilinearStiff = NULL;
     
        //transform the type
        SubTensorType tensorType;
        String2Enum(subType_,tensorType);
       
        bilinearStiff = new linElastInt(actSDMat, tensorType);

        return bilinearStiff;
      }


    BaseForm * FluidMechPDE::GetStiffIntegrator_UV(Double density,
                                                   Double kinematicViscosity) {
      BaseForm * bilinearStiff = NULL;
      if (subType_ == "plane" || subType_ == "3d"){
    	  if( stabilizationType_ != "none" && approxType_ == "lagrange" )
    		  bilinearStiff = new FluidMechPlaneStiffInt_UV(density,kinematicViscosity,stabilParams_,movingMesh_,stabilizationType_);
    	  else if( stabilizationType_ == "none" && approxType_ == "taylorHood" )
    		  bilinearStiff = new FluidMechPlaneMixedStiffInt_UV(density,kinematicViscosity,stabilParams_,movingMesh_);
    	  else
    		  Info->Error("Unknown discretizationType in fluidMech PDE! ",__FILE__,__LINE__);
      }
      else
        Info->Error("Unknown subtype in fluidMech PDE! ",__FILE__,__LINE__);
      bilinearStiff->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
      bilinearStiff->SetGridSolution( dynamic_cast<NodeStoreSol<Double>&>(*gridSol_ ));
      return bilinearStiff;
    }

    BaseForm * FluidMechPDE::GetStiffIntegrator_PQ(Double density,Double kinematicViscosity) {
      BaseForm * bilinearStiff = NULL;
      if ((subType_ == "plane" || subType_ == "3d") )
        bilinearStiff = new FluidMechPlaneStiffInt_PQ(density,kinematicViscosity,stabilParams_,movingMesh_,stabilizationType_);
      else
        Info->Error("Unknown subtype in fluidMech PDE! ",__FILE__,__LINE__);
      bilinearStiff->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
      bilinearStiff->SetGridSolution( dynamic_cast<NodeStoreSol<Double>&>(*gridSol_ ));
      return bilinearStiff;
    }

    BaseForm * FluidMechPDE::GetStiffIntegrator_UQ(Double density,Double kinematicViscosity) {
      BaseForm * bilinearStiff = NULL;
      if (subType_ == "plane" || subType_ == "3d")
    	  if( stabilizationType_ != "none" && approxType_ == "lagrange" )
    		  bilinearStiff = new FluidMechPlaneStiffInt_UQ(density,kinematicViscosity,stabilParams_,movingMesh_,stabilizationType_);
    	  else if( stabilizationType_ == "none" && approxType_ == "taylorHood" )
    		  bilinearStiff = new FluidMechPlaneMixedStiffInt_UQ(density,kinematicViscosity,stabilParams_,movingMesh_);
    	  else
    		  Info->Error("Unknown discretizationType in fluidMech PDE! ",__FILE__,__LINE__);
      else
        Info->Error("Unknown subtype in fluidMech PDE! ",__FILE__,__LINE__);
      bilinearStiff->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
      bilinearStiff->SetGridSolution( dynamic_cast<NodeStoreSol<Double>&>(*gridSol_ ));
      return bilinearStiff;
    }

    BaseForm * FluidMechPDE::GetStiffIntegrator_PV(Double density,Double kinematicViscosity) {
      BaseForm * bilinearStiff = NULL;
      if (subType_ == "plane" || subType_ == "3d")
    	  if( stabilizationType_ != "none" && approxType_ == "lagrange" )
    		  bilinearStiff = new FluidMechPlaneStiffInt_PV(density, kinematicViscosity,stabilParams_,movingMesh_,stabilizationType_);
    	  else if( stabilizationType_ == "none" && approxType_ == "taylorHood" )
    		  bilinearStiff = new FluidMechPlaneMixedStiffInt_PV(density,kinematicViscosity,stabilParams_,movingMesh_);
    	  else
    		  Info->Error("Unknown discretizationType in fluidMech PDE! ",__FILE__,__LINE__);
      else
        Info->Error("Unknown subtype in fluidMech PDE! ",__FILE__,__LINE__);
      bilinearStiff->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
      bilinearStiff->SetGridSolution( dynamic_cast<NodeStoreSol<Double>&>(*gridSol_ ));
      return bilinearStiff;
    }

    BaseForm * FluidMechPDE::GetMassIntegrator_UV(Double density,Double kinematicViscosity) {

      BaseForm * bilinearMass = NULL;

      if (subType_ == "plane" || subType_ == "3d") {
    	  if( stabilizationType_ != "none" && approxType_ == "lagrange" )
    		  bilinearMass = new FluidMechPlaneMassInt_UV(density,kinematicViscosity,stabilParams_,movingMesh_,stabilizationType_);
    	  else if( stabilizationType_ == "none" && approxType_ == "taylorHood" )
    		  bilinearMass = new FluidMechPlaneMixedMassInt_UV(density,kinematicViscosity,stabilParams_,movingMesh_);
    	  else
    		  Info->Error("Unknown discretizationType in fluidMech PDE! ",__FILE__,__LINE__);
      }
      else 
        Info->Error("Unknown subtype in fluidMech PDE! ",__FILE__,__LINE__);

      bilinearMass->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
      bilinearMass->SetGridSolution( dynamic_cast<NodeStoreSol<Double>&>(*gridSol_ ));
      return bilinearMass;
    }

    BaseForm * FluidMechPDE::GetMassIntegrator_UQ(Double density,Double kinematicViscosity) {

      BaseForm * bilinearMass = NULL;

      if (subType_ == "plane" || subType_ == "3d") {
        bilinearMass = new FluidMechPlaneMassInt_UQ(density,kinematicViscosity,stabilParams_,movingMesh_,stabilizationType_);
      }
      else 
        Info->Error("Unknown subtype in fluidMech PDE! ",__FILE__,__LINE__);

      bilinearMass->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
      bilinearMass->SetGridSolution( dynamic_cast<NodeStoreSol<Double>&>(*gridSol_ ));
      return bilinearMass;
    }

    void FluidMechPDE::DefineSolveStep()
    {
      solveStep_ = new SolveStepFluidMech(*this);

//      solveStep_ = new StdSolveStep(*this);
    }




    // ======================================================
    // COUPLING SECTION
    // ======================================================


    void FluidMechPDE::InitCoupling(PDECoupling * Coupling)
    {
        isIterCoupled_ = true;
        ptCoupling_   = Coupling;

        const UInt numOutputCouplings = ptCoupling_->GetNumOutputCouplings();

        StdVector<StdVector<UInt> > elemNodeToCouplingNode_tmp;
        elemNodeToCouplingNode_.Resize(numOutputCouplings);
        elemNodeToCouplingNode_.Init();

        // input couplings
        for (UInt i=0; i<ptCoupling_->GetNumInputCouplings(); i++) {
//          // check for input mechanic velocity
//          if (ptCoupling_->GetInputQuantity(i) == MECH_VELOCITY) {
//            numCouplingBcs_ += (dim_ * ptCoupling_->GetInputNumNodes(i));
//          }
           // check for input smooth velocity
          if (ptCoupling_->GetInputQuantity(i) == SMOOTH_VELOCITY) {
            numCouplingBcs_ += (dim_ * ptCoupling_->GetInputNumNodes(i));
          }
        }

        for (UInt i = 0; i < numOutputCouplings; i++) {
          if (ptCoupling_->GetOutputQuantity(i) == FLUIDMECH_FORCE) {
            // Intialize the memory of the coupling values
            ptCoupling_->CreateCouplingVector(i,isComplex_);

            // now since we need a incremental formulation, 
            //  initialize some necessary vectors
            isIncrFormulation_ = true;
          }

          else if (ptCoupling_->GetOutputQuantity(i) == ACOU_RHSVAL) {
            // Intialize the memory of the coupling values
            ptCoupling_->CreateCouplingVector(i,isComplex_); 
            
            //now since we need a incremental formulation, initialize some necessary vectors
            isIncrFormulation_ = true;
          }
          else if (ptCoupling_->GetOutputQuantity(i) == ACOUSURF_RHSVAL) {
            // Intialize the memory of the coupling values
            ptCoupling_->CreateCouplingVector(i,isComplex_); 
            
            //now since we need a incremental formulation, initialize some necessary vectors
            isIncrFormulation_ = true;
          }
          else {
            EXCEPTION("unknown coupling quantitiy");
          }
          
        }
    }


  void FluidMechPDE::CalcInputCoupling() {

    std::string errMsg;
    StdVector<UInt> * nodes;
    CFSVector * val;
    Integer eqnNr;
    UInt couplingDof;
    // Determine maximal allowed equation number for algebraic system
    Integer maxAllowedEqn = (Integer)algsys_->GetDimension();

    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == false)
      return;

    // Reset counter for boundary conditions
    couplingBCsCounter_ = 0;
 

    // Outer loop over all INPUT coupling terms
    for (UInt i=0; i<ptCoupling_->GetNumInputCouplings(); i++) {

      //    ptCoupling_ = &ptCoupling_[i];
      ptCoupling_->GetInputValues(i, val);
      couplingDof = ptCoupling_->GetInputDof(i);
    
      // Up to now, Coupling is only possible with
      // Real valued solutions
      Vector<Double> const & help = dynamic_cast<Vector<Double>&>(*val);

      switch(ptCoupling_->GetInputType(i)) {
        
        // -------------------
        // COORDINATE COUPLING
        // -------------------
      case COORD: 

        ptCoupling_->GetInputNodes(i, nodes);
//         LOG_DBG(fluidMechPDE) << "Input SMOOTH_DISPLACEMENT\n"
//                               << "couplingnodes:\n" << (*nodes)
//                               << "\n values:\n" << help << std::endl;
        ptgrid_->SetNodeOffset( *nodes , help );
      
        break;

        // -------------------
        // RHS COUPLING
        // -------------------
      case RHS:
        ptCoupling_->GetInputNodes(i, nodes);

//         LOG_DBG(fluidMechPDE) << "Input FLUIDMECH RHS????\n"
//                               << "couplingnodes:\n" << (*nodes)
//                               << "\n values:\n" << help << std::endl;
       
        //for (UInt dof=0; dof<ptCoupling_->GetInputDof(i); dof++)
        for ( UInt dof = 0; dof < couplingDof; dof++ ) {
          for ( UInt j = 0; j < nodes->GetSize(); j++ ) {
            eqnNr = eqnMap_->GetNodeEqn((*nodes)[j],dof+1);
            if ( eqnNr != 0 && eqnNr <= maxAllowedEqn ) {
              algsys_->SetNodeRHS( help[ dof + couplingDof * j ], pdeId_,
                                   eqnNr);
            }

#ifdef DEBUG
            else if ( eqnNr > maxAllowedEqn ) {
              (*debug) << "FluidMechPDE::CalcInputCoupling: "
                       << "(" << pdename_ << ") "
                       << "Refused to pass "
                       << "eqnNr = " << eqnNr << " to SetNodeRHS(), since "
                       << "it execeeds numLastFreeDof = " << maxAllowedEqn
                       << std::endl;
            }
            else if ( eqnNr == 0 ) {
              (*debug) << "FluidMechPDE::CalcInputCoupling: "
                       << "(" << pdename_ << ") "
                       << "Refused to pass "
                       << "eqnNr = " << eqnNr << " to SetNodeRHS(), since "
                       << "it is fixed by hom. Dirichlet BC" << std::endl;
            }
#endif

          }
        }
        break;

        // -----------------------
        // InhomDirichlet COUPLING
        // -----------------------
      case ID_BC:
        // How do we want to treat inhomogeneous Dirichlet boundary
        // conditions?
        {
          //check whether penalty formulation is used
          // with elimination this does not work currentlly
          bool usePenalty = true;
          std::string aux;
          StdVector<std::string> keyVec;
          StdVector<std::string> attrVec;
          StdVector<std::string> valVec;
//           keyVec  = "linearSystems", "system", "setup", "idbcHandling";
//           attrVec = "", "name", "";
//           valVec  = "", pdename_, "";
//           params->Get( keyVec, attrVec, valVec, aux );
//           usePenalty = aux == "penalty" ? true : false;
          usePenalty = true;
          Info->PrintF( pdename_, "Treating IDBCs using '%s' approach\n",
                        aux.c_str() );
          if ( usePenalty == false ) {
            EXCEPTION("Cannot use inhom. Dirichlet coupling together with "
                     << "IDBC elimination, since the equation numbering "
                     << "object does currently not have the information "
                     << "required to put those values at the end of the "
                     << "equation number interval! Please set idbcHandling="
                     << '"' << "penalty" << '"' << " in your xml-file");
          }
        }

        // Set flag that the boundary conditions have to be incorporated
        updateCouplingBCs_ = true;

        ptCoupling_->GetInputNodes(i, nodes);

//         LOG_DBG(fluidMechPDE) << "Input FLUIDMECH_VELOCITY\n"
//                               << "couplingnodes:\n" << (*nodes)
//                               << "values:\n" << help << std::endl;
	
        for ( UInt dof = 0; dof < dim_; dof++ ) {
          for ( UInt j = 0; j < nodes->GetSize();
                j++, couplingBCsCounter_++) {

            //eqnNr = eqnMap_->GetNodeEqn((*nodes)[j],dof+1 );
	    eqnNr = eqnMap_->GetNodeEqn(*results_[0],(*nodes)[j],dof+1);

            if (eqnNr==0) {
              EXCEPTION( "The specified coupling node has no equation number");
            }

            algsys_->SetDirichlet( pdeId_, eqnNr,
                                   help[dof+j*dim_] );
          }
        }
        break;
          
      case MAT:
        Error( "Not implemented yet", __FILE__, __LINE__ );
        break;

      case GRID_VEL:
        ptCoupling_->GetInputNodes(i, nodes);

//         LOG_DBG(fluidMechPDE) << "Input GRID_VELOCITY\n"
//                               << "couplingnodes:\n" << (*nodes)
//                               << "values:\n" << help << std::endl;

        //gridSol_->CouplingToNodeSolution(help,(*nodes));
        for ( UInt dof = 0; dof < dim_; dof++ ) {
          for ( UInt j = 0; j < nodes->GetSize();j++) {
            eqnNr = eqnMap_->GetNodeEqn(*results_[0],(*nodes)[j],dof+1);
            //gridSol_->data_[abs(eqnNr-1)]=help[dof+j*dim_];
            gridSol_->Set(eqnNr,help[dof+j*dim_]);
          }
        }
        break;

      }
      //       for (UInt ii=0; ii<help.GetSize(); ii++)
      //         Info->PrintF( pdename_, "help[%d]=%e\n", ii, help[ii]);

    }
  }
  
    void FluidMechPDE::CalcOutputCoupling()
    {
        UInt dof;
        SolutionType quantity;
        StdVector<Elem*> * couplingElems = NULL;
        StdVector<UInt> * couplingNodes = NULL;
        CFSVector * temp_values = NULL;
        CFSVector * temp_oldValues = NULL;
        //UInt regionCount = 0;

//        Info->PrintF( "FACTORIZATION", "The fluid mechanical forces will be factorized by forceFactor = %e\n",forceFac_);
      
        // at first, check if this PDE is iterative coupled
        if (isIterCoupled_ == false)
          return;

        // loop over all output coupling quantities
        for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {

          //Welche Größe ist die Koppelgröße??
          //gibt einen Solution Type zurück
          quantity = ptCoupling_->GetOutputQuantity(i);
          ptCoupling_->GetOutputValues(i, temp_values);
          ptCoupling_->GetOutputOldValues(i, temp_oldValues);

          // hard coded cast, since coupling is only possible with
          // real valued entries
          Vector<Double> * values = dynamic_cast<Vector<Double>*>(temp_values);
          Vector<Double> * oldValues = dynamic_cast<Vector<Double>*>(temp_oldValues);
          Vector<Double> auxValues;
          auxValues.Resize(values->GetSize());
          //Info->PrintF( pdename_, "Values_Size=%d\n", values->GetSize());
          //for(UInt i=0; i< values->GetSize(); i++)
          //Info->PrintF( pdename_, "Values[%d]=%d\n", i, (*values)[i]);


          switch(ptCoupling_->GetOutputType(i)) {
          case NODE:
            
            ptCoupling_->GetOutputNodes(i, couplingNodes);
    	//for (UInt ii=0; ii<(*couplingNodes).GetSize(); ii ++)
    	//Info->PrintF( pdename_, "couplingNodes[%d]=%d\n", ii, (*couplingNodes)[ii]);

            if (quantity == FLUIDMECH_FORCE) {
              ptCoupling_->GetOutputElements(i, couplingElems);
              dof = ptCoupling_->GetOutputDof(i);
              CalcMechCouplingRHS(couplingElems, *couplingNodes, auxValues, dof);
              //(*values)=(auxValues*forceFac_) + ((*oldValues)*(1.0-forceFac_));
              (*values)=(auxValues*forceFac_);
              //std::cout << "Force:\n" << (*values) << std::endl;
              Info->PrintF( pdename_, "Coupling Output Type=Force\n");
//              LOG_DBG(fluidMechPDE) << "Output FLUIDMECH_FORCE\n"
//                                    << "couplingnodes:\n" << (*couplingNodes)
//                                    << "values:\n" << (*values) << std::endl;

            }
            if (quantity == ACOU_RHSVAL) {
                if(saveAcouSrc_ == true)
                  {
                    acou_src_.NodeSolutionToCoupling(*values, *couplingNodes, 1);
                    Info->PrintF( pdename_, "Coupling Output Type=acouRHSval\n");
//                    LOG_DBG(fluidMechPDE) << "Output ACOU_RHSVAL\n"
//                                          << "couplingNodes:\n" << (*couplingNodes)
//                                          << "values:\n" << (*values) << std::endl;
                  }else
                  Error("Sorry, but acouRHSval is not defined as a NodeResult in your XML, please do so", __FILE__,__LINE__);
                  
            } 
            if (quantity == ACOUSURF_RHSVAL) {
              ptCoupling_->GetOutputElements(i, couplingElems);
              dof = ptCoupling_->GetOutputDof(i);
              CalcAcouSurfSourceCouplingRHS(couplingElems, *couplingNodes, (*values), dof);

              Info->PrintF( pdename_, "Coupling Output Type=acouSurfRHSVal\n");
//              LOG_DBG(fluidMechPDE) << "Output FLUIDMECH_FORCE\n"
//                                    << "couplingnodes:\n" << (*couplingNodes)
//                                    << "values:\n" << (*values) << std::endl;
            } 
            break;

          case ELEM:
            EXCEPTION("No Element coupling output");
          }
        }

    }

    void FluidMechPDE::GetPresSolVecOfElement( Vector<Double>& elemSol,
                                               StdVector<UInt>& connecth ) {

      // fluidMech pressure of element nodes
      elemSol.Resize(1 * connecth.GetSize());
      elemSol.Init(0);
      Integer eqnNr; 
      //UInt presDof;

      NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
      Vector<Double> & sol = solhelp->GetAlgSysVector();
  
      for(UInt actNode=0; actNode<connecth.GetSize(); actNode++) {
        //eqnNr = eqnMap_->GetNodeEqn(connecth[actNode],presDof);
        eqnNr = eqnMap_->GetNodeEqn(*results_[1],connecth[actNode],1);
        if (eqnNr!= 0) {
          elemSol[actNode] = sol[(abs(eqnNr-1))];
          //elemSol[actNode] = presSol[(abs(eqnNr-1))];
        }
        else {
          elemSol[actNode] = 0.0;
        }
      }
    }

    template <class TYPE>
      void FluidMechPDE::GetVeloSolMatOfElement( Matrix<TYPE>& elemSol,
                                                 StdVector<UInt>& connecth ) {

      // fluidMech velocities of element nodes
      elemSol.Resize(dim_,connecth.GetSize());
      elemSol.Init(0);
      Integer eqnNr; 
      //UInt veloEndDof;

      // Determine maximal allowed equation number for algebraic system
      Integer maxAllowedEqn = (Integer)algsys_->GetDimension();

      NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
      Vector<Double> & sol = solhelp->GetAlgSysVector();
  
      for(UInt actNode=0; actNode<connecth.GetSize(); actNode++) {
        for(UInt actDof=0; actDof<dim_; actDof++) {

          //eqnNr = eqnMap_->GetNodeEqn(connecth[actNode],actDof+1);
          eqnNr = eqnMap_->GetNodeEqn(*results_[0],connecth[actNode],actDof+1);
          if ( eqnNr != 0 ) {
            elemSol[actDof][actNode] = sol[eqnNr-1];
            //elemSol[actDof][actNode] = veloSol[eqnNr-1];
          }
          else if ( eqnNr <= maxAllowedEqn ) {
            elemSol[actDof][actNode] = 0.0;
          }
          else {
            EXCEPTION( "wrong equation Number");
          }
        }
      }
    }

    void FluidMechPDE::GetVeloDeriv1SolVecOfElement( Vector<Double>& sol,
                                                     const EntityIterator& it) {
      
     
      StdVector<Integer> eqns;
      eqnMap_->GetEqns( eqns, *results_[0], it );
      sol.Resize( eqns.GetSize() );
      sol.Init( 0.0 ); 
      
      if (  analysistype_ == TRANSIENT) {
        const Vector<Double> & sol_der1 = getS1();
        
        for( UInt i = 0; i < eqns.GetSize(); i++ ) {
          if ( eqns[i] != 0 ) {
            sol[i] = sol_der1[abs(eqns[i])-1];
          } else {
            sol[i] = 0.0;
          }
        }
      }
    }

    void FluidMechPDE::AcouSourceCalc()
    {
      /******************************************************************************************
       * Hier die Akustik-Flow Kopplung
       ******************************************************************************************/
      //Hier die Funktion CalcElemVector4Quad aufrufen
      
      if(saveAcouSrc_==true){

        UInt ii, elsize = 0;
        Integer eqnNr;
        
        StdVector<Elem*> elemssd;
        StdVector<UInt> connecth;
        StdVector<Integer> connect_PDE;
        Vector<Double>  Result;
        
        Matrix<Double> ptCoordNodes;
        Matrix<Double> flowdata;
        
        BaseFE * ptEl;
      
        //Reset the values in nodestresol, due to the usage of the function "Add"
        acou_src_.Init(0.0);
        
        // loop over all subdomains
        for (UInt isd=0; isd<subdoms_.GetSize(); isd++)
          {
            //Get the density of the current subdomain
            Double density;
            BaseMaterial * actMat = materials_[subdoms_[isd]];
            actMat->GetScalar(density,DENSITY,REAL);
            
            ElemList actSDList(ptgrid_ );
            actSDList.SetRegion( subdoms_[isd] );
            EntityIterator it = actSDList.GetIterator();
            //           if(ptgrid_->RegionIdToName(subdoms_[isd]) == "fluid"){
            // loop over all elements of subdomain
            for ( it.Begin(); !it.IsEnd(); it++ ) {
              //Pointer to the reference element
              const Elem * actCoupleElem = it.GetElem();
              ptEl = actCoupleElem->ptElem;
              
              //generate the object LinearFlowNoiseInt which
              //is defined in linearForm.cc
              LinearFlowNoiseInt * linear_load = new LinearFlowNoiseInt(ptEl);
                
              //Number of Nodes of the Element
              elsize=(actCoupleElem->connect).GetSize();
              //the length of the connecth is set to the number of nodes
              //because it stores the nodesnumbers of the element
              connecth.Resize(elsize);
              for (ii=0; ii<elsize; ii++)
                connecth[ii]=(actCoupleElem->connect)[ii];
              
              
              //Schreibt die Koordinaten der Knoten des Elements in ptCoordNodes
              Matrix<Double> ptCoordNodes;
              ptgrid_->GetElemNodesCoord(ptCoordNodes, connecth);      
                
              //as sad as it is, but we have to bring the fluid calculation results
              //in Matrix notation here :-(
              GetVeloSolMatOfElement(flowdata, connecth);
              //in case, that we really use the sorting funtionality, for the results, in CalcElemVector4Quad
              //we have to add one row for the pressure!
              //elem_data.AddRow(tmp,0);
              //flowdata = elem_data;
              
              //Now we have all the stuff needed to calculate the quadrupole source term,
              //CalcElemVector4Quad, gets the coordinates of the Nodes, the Node Numbers, the flowdata
              //and returns Result
//              LOG_DBG(fluidMechPDE) << "The current Element is:\n" << actCoupleElem->elemNum << std::endl;
              linear_load->CalcElemVec_withCFSVel(ptCoordNodes, flowdata, Result, density);
              // Now the results are written out as a gid result
              for(ii = 0; ii< ptEl->GetNumNodes(); ii++)
                { 
                  eqnNr = eqnMap_->GetNodeEqn(*results_[1],connecth[ii],1);
                  acou_src_.Add(eqnNr,Result[ii]);
                }
              delete linear_load;
            }
            //          }
          }    
      }
    }

    
    void FluidMechPDE::CalcMechCouplingRHS( StdVector<Elem*> * couplingElems, 
                                            StdVector<UInt> & couplingNodes,
                                            Vector<Double>& elemCouplingSols,
                                            UInt couplingdof ) {
      Double density = 0.0;
      Double dynamicViscosity = 0.0;
      Double kinematicViscosity = 0.0;
      Double sign = 0.0;
      Double auxPres;
      Integer matIndex = -1;
      Elem * ptVolElem = NULL;
      Matrix<Double> ptCoord, ptCoordVol, elemMat;
      Vector<Double> presSol, pressureForceOnNode, normal;
      Matrix<Double> veloSol;
      Vector<Double> Eps_xx_ForceOnNode, Eps_yy_ForceOnNode, Eps_zz_ForceOnNode;
      Vector<Double> Eps_xy_ForceOnNode, Eps_yz_ForceOnNode, Eps_xz_ForceOnNode;
      Vector<Double> Eps_xx_Sol, Eps_yy_Sol, Eps_zz_Sol, Eps_xy_Sol, Eps_yz_Sol, Eps_xz_Sol;
            
      elemCouplingSols.Init(0.0);
      for (UInt actElem=0; actElem<couplingElems->GetSize(); actElem++) {
        
        // Perform cast from volume element to surface element, since
        // fluidMech-Mech coupling makes only sense on surface elements
        SurfElem * actCoupleElem = 
          dynamic_cast<SurfElem*> ((*couplingElems)[actElem]);

        if (actCoupleElem == NULL) {
          EXCEPTION( "No elements found for coupling!");
        }
        
        BaseFE * ptElem = actCoupleElem->ptElem;

        StdVector<UInt> & connecth = actCoupleElem->connect;
        ptgrid_->GetElemNodesCoord( ptCoord, connecth, true );

        Double surfArea;
        surfArea=ptElem->CalcVolume(ptCoord, isaxi_);

        // Try to find according region for first neighbouring volume
        // element of the surface element
        matIndex = subdoms_.Find(actCoupleElem->ptVolElem1->regionId);
        
        // If first volume element does not belong to fluidMech PDE, try the
        // second one
        if ( matIndex == -1 ) {
          matIndex = subdoms_.Find(actCoupleElem->ptVolElem2->regionId);
          ptVolElem = actCoupleElem->ptVolElem2;
          sign = actCoupleElem->normalSign;

        } 
        else {
          ptVolElem = actCoupleElem->ptVolElem1;
          sign = -1.0 * actCoupleElem->normalSign;
        }
        
        if ( matIndex == -1) {
          EXCEPTION("FluidMechPDE::CalcMechCouplingRHS: The two volume element neighbours of surface element Nr.actCoupleElem->elemNum do not belong to my regions!");
        }

        // Assign correct density, dynamic viscosity and kinematic viscosity
        materials_[subdoms_[matIndex]]->GetScalar(density,DENSITY,REAL);
        materials_[subdoms_[matIndex]]->GetScalar(dynamicViscosity,DYNAMIC_VISCOSITY,REAL);
        materials_[subdoms_[matIndex]]->GetScalar(kinematicViscosity,KINEMATIC_VISCOSITY,REAL);

        std::auto_ptr<BaseForm> bilinear_mass = std::auto_ptr<BaseForm>(new MassInt(1.0,1,isaxi_,coordUp_ ));

        ElemList temp( ptgrid_ );
        temp.SetElement( actCoupleElem );
        EntityIterator it = temp.GetIterator();
        bilinear_mass->CalcElementMatrix(elemMat, it, it);
        //delete bilinear_mass;     

        //fluidMech      
        GetPresSolVecOfElement(presSol, connecth);
//        LOG_DBG2(fluidMechPDE) << "Pressure values for all nodes of the act. surface element\n" 
//                               << "Pressure:\n"<< presSol << std::endl
//                               << "Nodes:\n" << connecth << std::endl;

        pressureForceOnNode = elemMat * presSol;
        
        // force has to be added on RHS with negative sign
        pressureForceOnNode *= (-1.0 * density);
//        LOG_DBG2(fluidMechPDE) << "Pressure force for all nodes of the act. surface element\n" 
//                               << "pressureForceOnNode:\n"<< pressureForceOnNode << std::endl;
        
        // the normal vector points outwards of the MECHANICAL domain
        // (see. Kaltenbacher, "Num. Sim. of Mechatr. Act. & Sens." chapter 8.2)
        ptgrid_->CalcSurfNormal(normal, *actCoupleElem);
        normal *= sign;

        //***********************************************************************************************//
//        LOG_DBG2(fluidMechPDE) << "surfaceElemNum=" << actCoupleElem->elemNum << ", surfArea=" << surfArea << std::endl
//                               << "volNeughElemNum=" << ptVolElem->elemNum 
//                               << ", neighRegionId=" << ptVolElem->regionId << std::endl
//                               << "normal:\n" << normal << std::endl;

        //transform the type
        SubTensorType type;
        String2Enum(subType_,type);
        if (type == PLANE)
        	type = PLANE_STRAIN;
        
        std::auto_ptr<FluidMechShearStress<Double> > shearStress =  
          std::auto_ptr<FluidMechShearStress<Double> >(new FluidMechShearStress<Double>(materials_[subdoms_[matIndex]],type) );

        BaseFE * ptVolFE = ptVolElem->ptElem;

        ElemList temp2( ptgrid_ );
        temp2.SetElement( ptVolElem );
        EntityIterator it2 = temp2.GetIterator();

        StdVector<UInt> & connecthVol = ptVolElem->connect;
        ptgrid_->GetElemNodesCoord( ptCoordVol, connecthVol, true );

        GetVeloSolMatOfElement(veloSol, connecthVol);
//        LOG_DBG2(fluidMechPDE) << "Velocity solution for all nodes of the act. volume element\n" 
//                               << "veloSol:\n"<< veloSol << std::endl
//                               << "Nodes:\n" << connecthVol << std::endl;

        shearStress->SetActElemSol(veloSol);

        Vector<Double> intPoint;
        ptVolFE -> GetCoordMidPoint(intPoint);
        shearStress->SetIntPoint(intPoint);

        Vector<Double> elemStress;
        elemStress.Resize(stressDim_);
        elemStress.Init(0);
        
        //calculates the stress
        shearStress->CalcShearStressVec(elemStress,1,it2 );
//        LOG_DBG2(fluidMechPDE) << "Shear stress within the neighbouring volume element\n" 
//                               << "elemShearStress:\n" << elemStress << std::endl;

        Matrix<Double> elemShearForce;
        elemShearForce.Resize(connecth.GetSize(),dim_);
        if(dim_==2){
          Eps_xx_Sol.Resize(connecth.GetSize());
          Eps_xx_Sol.Init(elemStress[0]);
          Eps_xx_ForceOnNode = elemMat * Eps_xx_Sol;
//          LOG_DBG2(fluidMechPDE) << "Eps_xx_So=\n"<< Eps_xx_Sol << std::endl
//   			       << "Eps_xx_ForceOnNode=\n"<< Eps_xx_ForceOnNode;

          Eps_yy_Sol.Resize(connecth.GetSize());
          Eps_yy_Sol.Init(elemStress[1]);
          Eps_yy_ForceOnNode = elemMat * Eps_yy_Sol;
//          LOG_DBG2(fluidMechPDE) << "Eps_yy_So=\n"<< Eps_yy_Sol << std::endl
//   			       << "Eps_yy_ForceOnNode=\n"<< Eps_yy_ForceOnNode;

          Eps_xy_Sol.Resize(connecth.GetSize());
          Eps_xy_Sol.Init(elemStress[2]);
          Eps_xy_ForceOnNode = elemMat * Eps_xy_Sol;
//          LOG_DBG2(fluidMechPDE) << "Eps_xy_So=\n"<< Eps_xy_Sol << std::endl
//   			       << "Eps_xy_ForceOnNode=\n"<< Eps_xy_ForceOnNode;

          for(UInt n=0; n<connecth.GetSize(); n++){
            elemShearForce[n][0]=Eps_xx_ForceOnNode[n]*normal[0] + Eps_xy_ForceOnNode[n]*normal[1];
            elemShearForce[n][1]=Eps_xy_ForceOnNode[n]*normal[0] + Eps_yy_ForceOnNode[n]*normal[1];
          }
        }
        else if(dim_==3){
          Eps_xx_Sol.Resize(connecth.GetSize());
          Eps_xx_Sol.Init(elemStress[0]);
          Eps_xx_ForceOnNode = elemMat * Eps_xx_Sol;

          Eps_yy_Sol.Resize(connecth.GetSize());
          Eps_yy_Sol.Init(elemStress[1]);
          Eps_yy_ForceOnNode = elemMat * Eps_yy_Sol;

          Eps_zz_Sol.Resize(connecth.GetSize());
          Eps_zz_Sol.Init(elemStress[2]);
          Eps_zz_ForceOnNode = elemMat * Eps_zz_Sol;

          Eps_yz_Sol.Resize(connecth.GetSize());
          Eps_yz_Sol.Init(elemStress[3]);
          Eps_yz_ForceOnNode = elemMat * Eps_yz_Sol;

          Eps_xz_Sol.Resize(connecth.GetSize());
          Eps_xz_Sol.Init(elemStress[4]);
          Eps_xz_ForceOnNode = elemMat * Eps_xz_Sol;

          Eps_xy_Sol.Resize(connecth.GetSize());
          Eps_xy_Sol.Init(elemStress[5]);
          Eps_xy_ForceOnNode = elemMat * Eps_xy_Sol;

          for(UInt n=0; n<connecth.GetSize(); n++){
            elemShearForce[n][0]=Eps_xx_ForceOnNode[n]*normal[0] + Eps_xy_ForceOnNode[n]*normal[1] + Eps_xz_ForceOnNode[n]*normal[2];
            elemShearForce[n][1]=Eps_xy_ForceOnNode[n]*normal[0] + Eps_yy_ForceOnNode[n]*normal[1] + Eps_yz_ForceOnNode[n]*normal[2];
            elemShearForce[n][2]=Eps_xz_ForceOnNode[n]*normal[0] + Eps_yz_ForceOnNode[n]*normal[1] + Eps_zz_ForceOnNode[n]*normal[2];
          }

        }
//        LOG_DBG2(fluidMechPDE) << "Shear stress force for all nodes of the act. surface element\n" 
//                               << "elemShearForce:\n" << elemShearForce << std::endl;

        //***********************************************************************************************//
       
        for (UInt actNode=0; actNode<connecth.GetSize(); actNode++) {
        	UInt nodePos = 0;
        	while(connecth[actNode] != couplingNodes[nodePos] && nodePos < couplingNodes.GetSize()) 
        		nodePos++;
        	for (UInt actDof=0; actDof < couplingdof ; actDof++) 
        		elemCouplingSols[nodePos*couplingdof +actDof] += (pressureForceOnNode[actNode] * normal[actDof] + elemShearForce[actNode][actDof]);
        }

//        LOG_DBG2(fluidMechPDE) << "Complete Force Vector of all coupling elements\n" 
//                               << "elemCouplingSols:\n" << elemCouplingSols << std::endl;

      }//end loop over elem
    }

    void FluidMechPDE::CalcAcouSurfSourceCouplingRHS( StdVector<Elem*> * couplingElems, 
                                                      StdVector<UInt> & couplingNodes,
                                                      Vector<Double>& elemCouplingSols,
                                                      UInt couplingdof ) {
      Double density = 0.0;
      Double sign = 0.0;
      Integer matIndex = -1;
      Elem * ptVolElem = NULL;
      Matrix<Double> ptCoord, elemMat;
      Vector<Double> normal, actElemAcouRHS;
      Vector<Double> veloDerivSol, nveloDerivSol;
            
      elemCouplingSols.Init(0.0);
      for (UInt actElem=0; actElem<couplingElems->GetSize(); actElem++) {
        
        // Perform cast from volume element to surface element, since
        // acouSurfSource coupling makes only sense on surface elements
        SurfElem * actCoupleElem = 
          dynamic_cast<SurfElem*> ((*couplingElems)[actElem]);

        if (actCoupleElem == NULL) {
          EXCEPTION( "No elements found for coupling!");
        }
        
        BaseFE * ptElem = actCoupleElem->ptElem;

        StdVector<UInt> & connecth = actCoupleElem->connect;
        ptgrid_->GetElemNodesCoord( ptCoord, connecth, true );

        //Double surfArea;
        //surfArea=ptElem->CalcVolume(ptCoord, isaxi_);

        // Try to find according region for first neighbouring volume
        // element of the surface element
        matIndex = subdoms_.Find(actCoupleElem->ptVolElem1->regionId);
        
        // If first volume element does not belong to fluidMech PDE, try the
        // second one
        if ( matIndex == -1 ) {
          matIndex = subdoms_.Find(actCoupleElem->ptVolElem2->regionId);
          ptVolElem = actCoupleElem->ptVolElem2;
          sign = -1.0*actCoupleElem->normalSign;

        } 
        else {
          ptVolElem = actCoupleElem->ptVolElem1;
          sign = actCoupleElem->normalSign;
        }
        
        if ( matIndex == -1) {
        	EXCEPTION("FluidMechPDE::CalcAcouSurfSourceCouplingRHS: The two volume element neighbours of surface element do not belong to my regions!");
        }

        materials_[subdoms_[matIndex]]->GetScalar(density,DENSITY,REAL);
        
        std::auto_ptr<BaseForm> bilinear_mass = std::auto_ptr<BaseForm>(new MassInt(1.0,isaxi_,coordUp_ ));

        ElemList temp( ptgrid_ );
        temp.SetElement( actCoupleElem );
        EntityIterator it = temp.GetIterator();
        bilinear_mass->CalcElementMatrix(elemMat, it, it);
        
        GetVeloDeriv1SolVecOfElement(veloDerivSol, it);
        nveloDerivSol.Resize(connecth.GetSize());

        // the normal vector points outwards of the MECHANICAL domain
        // (see. Kaltenbacher, "Num. Sim. of Mechatr. Act. & Sens." chapter 8.2)
        ptgrid_->CalcSurfNormal(normal, *actCoupleElem);
        normal *= sign;
        for (UInt actNode=0; actNode < connecth.GetSize(); actNode++)
          for (UInt actDof=0; actDof<dim_; actDof++)
            nveloDerivSol[actNode] += veloDerivSol[actDof + actNode*dim_] * normal[actDof];


        actElemAcouRHS = elemMat * nveloDerivSol;
        
        // force has to be added on RHS with negative sign
        actElemAcouRHS *= (-1.0 * density);
        


        for (UInt actNode=0; actNode<connecth.GetSize(); actNode++) {
          UInt nodePos = 0;
            
          while(connecth[actNode] != couplingNodes[nodePos] &&
                nodePos < couplingNodes.GetSize()) {
            nodePos++;
          }
          elemCouplingSols[nodePos] += actElemAcouRHS[actNode];
        }
      }//end loop over elem
    }

    bool FluidMechPDE::HasOutput(SolutionType output)
    {

      if (output == FLUIDMECH_VELOCITY || output == FLUIDMECH_PRESSURE || output == FLUIDMECH_FORCE || output == ACOU_RHSVAL || output == ACOUSURF_RHSVAL)
        return true;

      return false;
    }



    // ======================================================
    // TIME STEPPING SECTION
    // ======================================================


    void FluidMechPDE::InitTimeStepping()
    {
    	// timestepping formulation
    	TS_alg_ = new Bdf2( algsys_);
//    	TS_alg_ = new Trapezoidal( algsys_);
    }

    // ***********************************************************************
    //   Obtain information on desired output quantities from parameter file
    // ***********************************************************************
    void FluidMechPDE::ReadSpecialResults() {

        //gridSol_ = NULL;     
        gridSol_->SetNumSolutions(1);
        gridSol_->SetNumNodes(numPDENodes_);
        gridSol_->SetSolutionTypeName(FLUIDMECH_VELOCITY,GRID_VELOCITY);
        gridSol_->SetNumDofs(dim_);
        gridSol_->SetResults( results_ );
        gridSol_->SetPtrEQNData( eqnMap_.get(), ptgrid_);
        gridSol_->SetRegions( subdoms_ );
        gridSol_->Init();

        if(saveAcouSrc_ == true){
          acou_src_.SetNumSolutions(1);
          acou_src_.SetSolutionTypeName(FLUIDMECH_PRESSURE,ACOU_RHSVAL);
          acou_src_.SetNumNodes(numPDENodes_);
          acou_src_.SetNumDofs(1);   //dim_);
          acou_src_.SetResults(results_);
          acou_src_.SetPtrEQNData( eqnMap_.get(), ptgrid_);
          acou_src_.SetRegions( subdoms_ );
          acou_src_.Init();
        }

    }

    void FluidMechPDE::DefineAvailResults() {

      //the size of the matrix with the stabilization parameters is currently
      //Num2DElements/Num3DElements long. 
      //Bec currently i have no idea how to access local elem numbers in the integrator
      stabilParams_.Resize(ptgrid_->GetNumElems()+1,4);
      stabilParams_.Init(-1.0);    	    	
    	
      // Check for subType
      StdVector<std::string> dispDofNames, stressDofNames;
    
      if( subType_ == "3d" ) {
        dispDofNames = "x", "y", "z";
        stressDofNames = "xx", "yy", "zz", "yz", "xz", "xy";

      } else if( subType_ == "plane" ) {
        dispDofNames = "x", "y";
        stressDofNames = "xx", "yy", "xy";

      } else if( subType_ == "axi" ) {
        dispDofNames = "r", "z";
        stressDofNames = "rr", "zz", "rz", "phiphi";

      }

      // === FLUIDMECHANIC VELOCITY ===
      shared_ptr<ResultInfo> vel(new ResultInfo);
      vel->resultType = FLUIDMECH_VELOCITY;
      vel->dofNames = dispDofNames;
      vel->unit = "m/s";
      vel->entryType = ResultInfo::VECTOR;

      // check if problem is lagrange, legendre or taylorHood
      if ( approxType_ == "lagrange" ) {	
    	  shared_ptr<AnsatzFct> fct(new LagrangeFct);
    	  vel->fctType = fct;
    	  vel->definedOn = ResultInfo::NODE;
      } else if ( approxType_ == "legendre" ) {
    	  shared_ptr<LegendreFct> fct(new LegendreFct);
    	  UInt order = myParam_->Get("order")->AsUInt();
    	  fct->SetIsoOrder( order );
    	  vel->fctType = fct;
    	  vel->definedOn = ResultInfo::PFEM;
      } else if (approxType_ == "taylorHood"){
    	  if (stabilizationType_=="none") {
    		  // define Legendre type
    		  shared_ptr<LegendreFct> fct(new LegendreFct);
    		  fct->SetIsoOrder( 2 );
    		  vel->fctType = fct;
    		  vel->definedOn = ResultInfo::PFEM;
    	  } else 
    		  Error( "stabilization combined with taylorHood needs to be implemented", __FILE__, __LINE__ );
      }
      else
    	  Error( "wrong approximation type", __FILE__, __LINE__ );	
      
      results_.Push_back( vel );
      availResults_.insert( vel);

      // === FLUIDMECHANIC PRESSURE ===
      shared_ptr<ResultInfo> pres(new ResultInfo);
      pres->resultType = FLUIDMECH_PRESSURE;
      pres->dofNames = "";
      pres->unit = "N/m^2";
      pres->entryType = ResultInfo::SCALAR;

      // check if problem is lagrange or legendre
      if ( approxType_ == "lagrange" ) {
    	  shared_ptr<AnsatzFct> fct(new LagrangeFct);
    	  pres->fctType = fct;
    	  pres->definedOn = ResultInfo::NODE;
      } else if ( approxType_ == "legendre" ){
    	  // define Legendre type
    	  shared_ptr<LegendreFct> fct(new LegendreFct);
    	  UInt order = myParam_->Get("order")->AsUInt();
    	  fct->SetIsoOrder( order );
    	  pres->fctType = fct;
    	  pres->definedOn = ResultInfo::PFEM;
      } else if (approxType_ == "taylorHood"){
    	  if (stabilizationType_=="none") {
    		  shared_ptr<LegendreFct> fct(new LegendreFct);
    		  fct->SetIsoOrder( 1 );
    		  pres->fctType = fct;
    		  pres->definedOn = ResultInfo::PFEM;
    	  } else 
    		  Error( "stabilization combined with taylorHood needs to be implemented", __FILE__, __LINE__ );
      }
      else
    	  Error( "wrong approximation type", __FILE__, __LINE__ );	

      
      results_.Push_back( pres );
      availResults_.insert( pres);

//      // === FLUIDMECHANIC STRESS ===
//      shared_ptr<ResultInfo> stress(new ResultInfo);    
//      stress->resultType = FLUIDMECH_STRESS;
//      stress->dofNames = stressDofNames;
//      stress->unit =  "N/m^2";
//      stress->entryType = ResultInfo::TENSOR;
//      stress->definedOn = ResultInfo::ELEMENT;
//      stress->fctType = shared_ptr<ConstFct>(new ConstFct() );
//      availResults_.insert( stress );
//
//      // === FLUIDMECHANIC STRAIN RATE===
//      shared_ptr<ResultInfo> strainRate(new ResultInfo);    
//      strainRate->resultType = FLUIDMECH_STRAINRATE;
//      strainRate->dofNames = stressDofNames;
//      strainRate->unit =  "";
//      strainRate->entryType = ResultInfo::TENSOR;
//      strainRate->definedOn = ResultInfo::ELEMENT;
//      strainRate->fctType = shared_ptr<ConstFct>(new ConstFct() );
//      availResults_.insert( strainRate );
//
//      // === FLUIDMECHANIC ENERGY ===
//      shared_ptr<ResultInfo> energy(new ResultInfo);    
//      energy->resultType = FLUIDMECH_ENERGY;
//      energy->dofNames = "";
//      energy->unit = "Ws";
//      energy->entryType = ResultInfo::SCALAR;
//      energy->definedOn = ResultInfo::REGION;
//      energy->fctType = shared_ptr<ConstFct>(new ConstFct() );
//      availResults_.insert( energy );
//    
//      // === ACOUSTIC SOURCES (Lighthill) ===
//      shared_ptr<ResultInfo> acouSrc(new ResultInfo);    
//      acouSrc->resultType = ACOU_RHS_LOAD;
//      acouSrc->dofNames = "";
//      acouSrc->unit = "";
//      acouSrc->entryType = ResultInfo::SCALAR;
//      acouSrc->definedOn = ResultInfo::NODE;
//      acouSrc->fctType = shared_ptr<ConstFct>(new ConstFct() );
//      availResults_.insert( acouSrc );
//
//
//      // === STABILIZATION PARAMETERS ===
//      shared_ptr<ResultInfo> stabilParam(new ResultInfo);    
//      stabilParam->resultType = FLUIDMECH_STABILPARAM;
//      stabilParam->dofNames = "";
//      stabilParam->unit = "";
//      stabilParam->entryType = ResultInfo::VECTOR;
//      stabilParam->definedOn = ResultInfo::ELEMENT;
//      stabilParam->fctType = shared_ptr<ConstFct>(new ConstFct() );
//      availResults_.insert( stabilParam );
    
      // === OPT_RESULT_1/2/3 ===
      // this is added via the optimization stuff in DesignSpace.
    }

    void FluidMechPDE::CalcResults( shared_ptr<BaseResult> result ) {
    
      switch (result->GetResultInfo()->resultType ) {
      
      case FLUIDMECH_VELOCITY:
        if( isComplex_ ) {
          ExtractResult<Complex>( result, sol_ );
        } else {
          ExtractResult<Double>( result, sol_ );
        }
        break;
      
      case FLUIDMECH_PRESSURE:
        if( isComplex_ ) {
          ExtractResult<Complex>( result, sol_ );
        } else {
          ExtractResult<Double>( result, sol_ );
        }
        break;
      
      case ACOU_RHS_LOAD:
        if( isComplex_ ) {
          ExtractRhsResult<Complex>( result, results_[0] );
        } else {
          ExtractRhsResult<Double>( result, results_[0] );
        }
        break;

      case FLUIDMECH_STRESS:
        if( isComplex_ ) {
          CalcStresses<Complex>( result );
        } else {
          CalcStresses<Double>( result );
        }
        break;

      case FLUIDMECH_STRAINRATE:
        if( isComplex_ ) {
          CalcStrainRates<Complex>( result );
        } else {
          CalcStrainRates<Double>( result );
        }
        break;

      case FLUIDMECH_ENERGY:
        if( isComplex_ ) {
          CalcEnergy<Complex>( result );
        } else {
          CalcEnergy<Double>( result );
        }
        break;

        // the actual case is given in the result info in result
      case OPT_RESULT_1:
      case OPT_RESULT_2:
      case OPT_RESULT_3:
        // design should work, this is checked in AvailabeResults()
        domain->GetErsatzMaterial()->ExtractResults(result, isComplex_);
        break;


      default:
        Warning( "Resulttype not computable by mechanic PDE",
                 __FILE__, __LINE__ );
      }
    }


    template <class TYPE>
      void FluidMechPDE::CalcStresses( shared_ptr<BaseResult> res ) {
    
      //get the correct bilinear form
      Vector<Double> intPoint;
  
      //transform the type
      SubTensorType type;
      String2Enum(subType_,type);

      Vector<TYPE> elemStress;
      elemStress.Resize(stressDim_);
      elemStress.Init(0);
    
      Result<TYPE> &  actRes = 
        dynamic_cast<Result<TYPE>&>(*res);
      EntityIterator it = actRes.GetEntityList()->GetIterator();
    
      Vector<TYPE> & actVal = actRes.GetVector();
      actVal.Resize( actRes.GetEntityList()->GetSize() * stressDim_ );
    
      // Fetch material: As we assume, that all elements belong to
      // one and the same region, we simply take the subdomain of the first 
      // element
      it.Begin();
      BaseMaterial* actSDMat = materials_[it.GetElem()->regionId];
      MechStressStrain<TYPE> *stress = new MechStressStrain<TYPE>(actSDMat,type);

      // loop over elements
      for ( it.Begin(); !it.IsEnd(); it++ ) {
        it.GetElem()->ptElem->GetCoordMidPoint(intPoint);
      
        //set element solution        
        Matrix<TYPE> elSol;
        sol_->GetElemSolutionAsMatrix(elSol, it);
        stress->SetActElemSol(elSol);

        //calculates the stress
        stress->SetIntPoint(intPoint);
        stress->CalcStressVec(elemStress,1,it);
        stress->UnsetIntPoint();
        for( UInt iDof = 0; iDof < stressDim_; iDof++ ) {
          actVal[it.GetPos()*stressDim_ + iDof] = elemStress[iDof];
        }
      }
    
      delete stress;
    }

    template <class TYPE>
      void FluidMechPDE::CalcStrainRates( shared_ptr<BaseResult> res ) {
      //transform the subType of the pde
      SubTensorType type;
      String2Enum(subType_,type);

      Vector<Double> intPoint;
      Vector<TYPE> elemStrainRate;
      elemStrainRate.Resize(stressDim_);
      elemStrainRate.Init(0);
    
      Result<TYPE> &  actRes = 
        dynamic_cast<Result<TYPE>&>(*res);
      EntityIterator it = actRes.GetEntityList()->GetIterator();
    
      Vector<TYPE> & actVal = actRes.GetVector();
      actVal.Resize( actRes.GetEntityList()->GetSize() * stressDim_ );
    
      // Fetch material: As we assume, that all elements belong to
      // one and the same region, we simply take the subdomain of the first 
      // element
      it.Begin();
      BaseMaterial* actSDMat = materials_[it.GetElem()->regionId];
      MechStressStrain<TYPE> *strainRate = new MechStressStrain<TYPE>(actSDMat,type);

      // loop over elements
      for ( it.Begin(); !it.IsEnd(); it++ ) {
        it.GetElem()->ptElem->GetCoordMidPoint(intPoint);
      
        StdVector<UInt> connecth = it.GetElem()->connect;

        //set element solution        
        Matrix<TYPE> elSol;
        GetVeloSolMatOfElement(elSol, connecth);

        strainRate->SetActElemSol(elSol);

        //calculates the element strainRate
        strainRate->SetIntPoint(intPoint);
        strainRate->CalcStrainVec(elemStrainRate,1,it);
        strainRate->UnsetIntPoint();
        for( UInt iDof = 0; iDof < stressDim_; iDof++ ) {
          actVal[it.GetPos()*stressDim_ + iDof] = elemStrainRate[iDof];
        }
      }
    
      delete strainRate;
    }




    // ======================================================
    // POSTPROCESSING  
    // ======================================================
    template<class TYPE>
      void FluidMechPDE::ComputeVolDefSurf( shared_ptr<BaseResult> vals ) {

    
      // convert result and get region iterator
      Result<TYPE> &  actSol = 
        dynamic_cast<Result<TYPE>&>(*vals);      
      EntityIterator regionIt = actSol.GetEntityList()->GetIterator();
    
      // resize vector
      Vector<TYPE> & actVal = actSol.GetVector();
      actVal.Resize( actSol.GetEntityList()->GetSize() );
    
      // get dof for direction
      std::string dir = volAboveDefSurfDir_[actSol.GetEntityList()];
      Integer dof = results_[0]->dofNames.Find( dir );
      if( dof < 0 ) {
        EXCEPTION( "ComputeVolDefSurf: dof = '" << dir << "' is not "
                   << "allowed! Must be one of 'x', 'y' or 'z'" );
      }
    
      Vector<TYPE> elemDisp, elemDispDof;
      TYPE actVolume;
    
      // Loop over regions
      for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {
        SurfElemList actSDList(ptgrid_ );
        actSDList.SetRegion( regionIt.GetRegion() );
        EntityIterator elemIt = actSDList.GetIterator();
      
        actVolume = 0.0;
        // Loop over elements
        for( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
        
          BaseFE * ptSurfEl = elemIt.GetSurfElem()->ptElem;
          StdVector<UInt> connecth = elemIt.GetSurfElem()->connect;
        
          Matrix<Double> ptSurfCoord;
          ptgrid_->GetElemNodesCoord( ptSurfCoord, connecth, false );
        
          GetSolVecOfElement( elemDisp, elemIt, results_[0] );
          elemDispDof.Resize(connecth.GetSize() );
          UInt actEntry = dof;
          for( UInt i = 0; i < connecth.GetSize(); i++ ) {
            elemDispDof[i] = elemDisp[actEntry];
            actEntry += dim_;
          }
          actVolume += ComputeVolElem<TYPE>( ptSurfEl, 
                                             ptSurfCoord, elemDispDof );
        
        }
        actVal[regionIt.GetPos()] = actVolume;
      }
    
    }

    template<class TYPE>
      TYPE FluidMechPDE::ComputeVolElem(BaseFE * surfEl, Matrix<Double>& surfCoord, 
                                        Vector<TYPE> disp) {


      TYPE elemVol, averageDis;
      UInt nrSurfNodes = surfEl->GetNumNodes();


      //compute average displacedment
      averageDis = 0;
      for (UInt i=0; i<nrSurfNodes; i++) {
        averageDis += disp[i]; 
      }
      averageDis /= (Double)nrSurfNodes;

      //compute the deformed volume    
      elemVol = averageDis * surfEl->CalcVolume(surfCoord,isaxi_);
      
      return elemVol;
    }


   
    void FluidMechPDE::ReadPressureLoads() {

      // try to get bcsAndLoads node
      ParamNode * bcNode = myParam_->Get("bcsAndLoads", false);
      if( !bcNode )
        return;
      StdVector<ParamNode*> presNodes = bcNode->GetList("pressure");
    
      // iterate over all pressure definitions
      std::string name, value, phase;
      for( UInt i = 0; i < presNodes.GetSize(); i++ ) {
      
        presNodes[i]->Get("name", name );
        presNodes[i]->Get("value", value );
        presNodes[i]->Get("phase", phase );
      
        pressSurf_.Push_back(ptgrid_->GetEntityList( EntityList::SURF_ELEM_LIST,
        		name, EntityList::REGION ) );
        pressVals_.Push_back( value );
        pressPhase_.Push_back( phase );

      }
    }


    template <class TYPE>
      void FluidMechPDE::CalcEnergy( shared_ptr<BaseResult> vals )
      {

        Matrix<Double> elemmat;  
        Vector<TYPE> help, elvel;
        TYPE energy ;

        Result<TYPE> &  actSol = 
          dynamic_cast<Result<TYPE>&>(*vals);      

        // resize vector
        Vector<TYPE> & actVal = actSol.GetVector();
        actVal.Resize( actSol.GetEntityList()->GetSize() );

        // loop over regions
        EntityIterator regionIt = actSol.GetEntityList()->GetIterator();
        for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {

          //get material
          BaseMaterial* actSDMat = materials_[regionIt.GetRegion()];

          // get bilinear stiffness integrator
          BaseForm * bilinear_stiff = GetStiffIntegrator(actSDMat, 
                                                         regionIt.GetRegion() );
          ElemList actSDList(ptgrid_ );
          actSDList.SetRegion( regionIt.GetRegion() );
          EntityIterator it = actSDList.GetIterator();
          energy = 0.0;
      
          // Loop over elements
          for ( it.Begin(); !it.IsEnd(); it++ ) {
            bilinear_stiff->SetAnsatzFct( results_[0]->fctType );
            bilinear_stiff->CalcElementMatrix(elemmat,it,it);
            sol_->GetElemSolution(elvel, it);
            help = elemmat * elvel;
            energy += ( help * elvel) * 0.5;
          } 
          actVal[regionIt.GetPos()] = energy;
          delete bilinear_stiff;    
        }
      }

    void FluidMechPDE::InitStabParams(){

    	if( stabilizationType_ != "none"){
    		//the size of the matrix with the stabilization parameters is currently
    		//Num2DElements/Num3DElements long. 
    		//Bec currently i have no idea how to access local elem numbers in the integrator
    		//std::cerr << "init\n"; 
    		for(Integer i=0; i<ptgrid_->GetNumElems(); i++) 
    			stabilParams_[i][3]=(-1.0);
    	}
    }

    void FluidMechPDE::PrintStabParams(){

    	if( stabilizationType_ != "none"){
    		LOG_DBG(fluidmechpde) << "Stabilization parameters\n"
    		<< stabilParams_ << std::endl;
    	}
    }

  
    void FluidMechPDE::RegisterFunctions() {
    
      typedef FctPointer<FluidMechPDE> FCPT;
      StdVector<ArgList> a;
      StdVector<FCPT*> pt;
      StdVector<std::string> name;

      // --- ReadRegionLoad ---
      a.Push_back();
      a.Last().RegisterParam( "name", ArgList::STRING );
      a.Last().RegisterParam( "value", ArgList::STRING );
      a.Last().RegisterParam( "dof", ArgList::STRING );
      a.Last().RegisterParam( "coordSysId", ArgList::STRING );
      a.Last().RegisterParam( "type", ArgList::STRING );
      pt.Push_back( new FCPT ( this, &FluidMechPDE::ReadRegionLoads ) );
      name.Push_back( "setRegionLoad" );

    
      // Now register all functions with scripting 
      for (UInt i = 0; i < pt.GetSize(); i++ ) {
        Script_RegisterFct(name[i], pt[i], a[i] );
      }          

    }

  
  } // end namespace CoupledField
