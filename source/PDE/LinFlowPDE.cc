// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <set>

#include "LinFlowPDE.hh"

#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionStabParams.hh"
#include "Domain/CoefFunction/CoefFunctionMeanFlowConvection.hh"
#include "Domain/CoefFunction/CoefFunctionMeanFlowConvection.hh"
#include "Domain/CoefFunction/CoefFunctionComplexToReal.hh"
#include "Domain/CoefFunction/CoefFunctionConversion.hh"
#include "Domain/CoefFunction/CoefFunctionDiagTensorFromScalar.hh"
#include "Forms/Operators/IdentityOperatorInVector.hh"

#include "Utils/StdVector.hh"
#include "Driver/SolveSteps/SolveStepElec.hh"
#include "Driver/Assemble.hh"
#include "Utils/ApproxData.hh"
#include "Utils/SmoothSpline.hh"
#include "Materials/Models/Hysteresis.hh"
#include "Materials/Models/Preisach.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/FeFunctions.hh"
#include "Domain/ElemMapping/ElemShapeMap.hh"

// used by Mortar coupling
#include "Domain/Mesh/NcInterfaces/MortarInterface.hh"
#include "FeBasis/HCurl/HCurlElems.hh" // needed for the surface operators
#include "Forms/Operators/SurfaceOperators.hh"
#include "Forms/Operators/SurfaceNormalStressOperator.hh"

// new integrator concept
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/BiLinForms/ADBInt.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/BUInt.hh"

#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/DivOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorProjected.hh"
#include "Forms/Operators/MultiIdOp.hh"
#include "Forms/Operators/LaplOp.hh"
#include "Forms/Operators/ConvectiveOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"
#include "Forms/Operators/IdentityOperatorNormalTrans.hh"
#include "Forms/Operators/StrainOperator.hh"
#include "Forms/Operators/SurfaceOperators.hh"
#include "Forms/Operators/SurfaceNormalStressOperator.hh"
#include "Domain/CoefFunction/CoefXpr.hh"


#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
//#include "Domain/Results/ExternalFieldFunctors.hh"

namespace CoupledField {

  DEFINE_LOG(linfluidmechpde, "pde.linfluidmech")

  // ***************
  //   Constructor
  // ***************
  LinFlowPDE::LinFlowPDE( Grid* grid, PtrParamNode paramNode,
                                      PtrParamNode infoNode,
                                      shared_ptr<SimState> simState, 
                                      Domain* domain )
  :SinglePDE( grid, paramNode, infoNode,simState, domain )
  {
    // ======================================================================
    // set solution information
    // ======================================================================
    pdename_          = "fluidMechLin";
    pdematerialclass_ = FLOW;
 
    nonLin_    = false;
    nonLinMaterial_ = false;
    isAlwaysStatic_ = false;
    isHeatPDECoupled_  = false;
    isCouplingFormulationSymmetric_ = false;

    //! Always use total Lagrangian formulation 
    updatedGeo_        = true;
 
    //! check, if a compressible formulation is specified
    isCompressible_ = false;
    std::string formulation = myParam_->Get("formulation")->As<std::string>();
    if ( formulation == "compressible")
        isCompressible_ = true;


    // Check the subtype of the problem
    paramNode->GetValue("subType", subType_);

    if(ptGrid_->GetDim() == 3) {
      subType_ = "3d";
    } else if(isaxi_) {
      subType_ = "axi";
    }
    if( !(subType_ == "3d" || subType_ == "plane" || subType_ == "axi") ) {
      EXCEPTION( "Subtype '" << subType_ << "' unknown for linearized fluid-mechanic physic" );
    }

    //type of geometry
    isaxi_ = ptGrid_->IsAxi();
    
      
    //get oder of FE basis functions
    presPolyId_ = myParam_->Get("presPolyId")->As<std::string>();
    velPolyId_  = myParam_->Get("velPolyId")->As<std::string>();
    lagrangeMultPolyId_ = myParam_->Get("presPolyId")->As<std::string>();

    //get order of integration
    presIntegId_ = myParam_->Get("presIntegId")->As<std::string>();
    velIntegId_  = myParam_->Get("velIntegId")->As<std::string>();
    lagrangeMultIntegId_ = myParam_->Get("presIntegId")->As<std::string>();

    // Obtain pressure surface regions
    PtrParamNode presSurfaceNode = myParam_->Get("presSurfaceList",
                                                 ParamNode::PASS);
    ParamNodeList regionNodes;
    if(presSurfaceNode) {
      regionNodes = presSurfaceNode->GetList("presSurf");
    }

    // output and set regions_
    for( UInt i = 0; i < regionNodes.GetSize(); i++ )
    {
      std::string regionName = regionNodes[i]->Get("name")->As<std::string>();
      RegionIdType actRegionId = ptGrid_->GetRegion().Parse(regionName);

      // Check, if region was already defined an issue a warning otherwise
      if( std::find(presSurfaces_.Begin(), presSurfaces_.End(), actRegionId )
          != presSurfaces_.End() )  {
        WARN( "The pressure surface '" << regionName
              << "' was already defined for PDE '" << pdename_
              << "'. Please remove duplicate entries." );
      }

      presSurfaces_.Push_back( actRegionId );
    }

    enableC2_ = false;
    enableC2_ = myParam_->Get("enableC2")->As<bool>();
    if ( enableC2_ ) WARN("ENABLING SECOND CONVECTIVE TERM (C2)");

    enableC3_ = myParam_->Get("enableC3")->As<bool>();
    if ( enableC3_ ) WARN("ENABLING Third CONVECTIVE TERM (C3 = V0.Grad(V0)");

    factorC1_ = 2.0;
    factorC1_ = myParam_->Get("factorC1")->As<Double>();

    enableGridVelC1_ = false;
    enableGridVelC1_ = myParam_->Get("enableGridVelC1")->As<bool>();
    if ( enableGridVelC1_ )
          std::cout << "\n ENABLING FIRST CONVECTIVE TERM CAUSED BY THE GRID VELOCITY (GridVelC1, ALE)\n" << std::endl;

    enableGridVelC2_ = false;
    enableGridVelC2_ = myParam_->Get("enableGridVelC2")->As<bool>();
    if ( enableGridVelC2_ )
          std::cout << "\n ENABLING SECOND CONVECTIVE TERM CAUSED BY THE GRID VELOCITY (GridVelC2, ALE)\n" << std::endl;

    // grid velocity coefFunction
    gridVelCoef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, ptGrid_->GetDim(), 1, isComplex_, true));


    // Check if the lagrange multiplier is needed
    PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    ParamNodeList velocityConstraintNodes = bcNode->GetList( "velocityConstraint" );
    useLagrangeMultVec_ = false;
    if( !velocityConstraintNodes.IsEmpty() ) {
      // enables usage of LAGRANGE_MULT
      useLagrangeMultVec_ = true;
    }

     // Check if the lagrange multiplier 1 is needed
    ParamNodeList scalarVelocityConstraintNodes = bcNode->GetList( "scalarVelocityConstraint" );
    useLagrangeMultScal_ = false;
    if( !scalarVelocityConstraintNodes.IsEmpty() ) {
      // enables usage of LAGRANGE_MULT_1
      useLagrangeMultScal_ = true;
    }
  }

  void LinFlowPDE::SetHeatPDECouplingFlags(bool useSymmtericForm) {
  	// Set flag for coupling
    isHeatPDECoupled_ = true;
    // Set flag whether to use symmetric form or not
    isCouplingFormulationSymmetric_= useSymmtericForm;
  }
  
  double LinFlowPDE::GetBalanceOfMomentumSign() const { return balanceOfMomentumSign_; }
  
  void LinFlowPDE::InitNonLin() {
    SinglePDE::InitNonLin();
  }

  void LinFlowPDE::DefineIntegrators() {

    RegionIdType actRegion;

    // Get FESpace and FeFunction of fluid velocity
    shared_ptr<BaseFeFunction> velFct = feFunctions_[FLUIDMECH_VELOCITY];
    shared_ptr<FeSpace> velSpace = velFct->GetFeSpace();

    // Get FESpace and FeFunction of fluid pressure
    shared_ptr<BaseFeFunction> presFct = feFunctions_[FLUIDMECH_PRESSURE];
    shared_ptr<FeSpace> presSpace = presFct->GetFeSpace();

    shared_ptr<BaseFeFunction> meanVelFct = meanFlowFeFct_;
    shared_ptr<FeSpace> meanVelSpace = meanVelFct->GetFeSpace();

    // Create coefficient functions for all fluid densities
    std::map< RegionIdType, PtrCoefFct > oneFuncs;
    std::set< RegionIdType > flowRegions;

    shared_ptr<BaseFeFunction> lagrangeMultFct = NULL;
    shared_ptr<FeSpace> lagrangeMultSpace = NULL;

    shared_ptr<BaseFeFunction> lagrangeMultFct_1 = NULL;
    shared_ptr<FeSpace> lagrangeMultSpace_1 = NULL;

    if ( useLagrangeMultVec_ ) {
      lagrangeMultFct = feFunctions_[LAGRANGE_MULT];
      lagrangeMultSpace = lagrangeMultFct->GetFeSpace();
    }
    if ( useLagrangeMultScal_ ) {
      lagrangeMultFct_1 = feFunctions_[LAGRANGE_MULT_1];
      lagrangeMultSpace_1 = lagrangeMultFct_1->GetFeSpace();
    }

    //  Loop over all regions
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion++){

      // Set current region and material
      actRegion = regions_[iRegion];

      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      // Get ParamNode for current region
      PtrParamNode curRegNode = myParam_->Get("regionList")
          ->GetByVal("region","name",regionName.c_str());

      // create new entity list and add it fefunction
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      velFct->AddEntityList( actSDList );
      presFct->AddEntityList( actSDList );
      meanVelFct->AddEntityList( actSDList );


      // Lagrange multiplier
      if ( useLagrangeMultVec_ ) {
        lagrangeMultFct->AddEntityList( actSDList );
        lagrangeMultSpace->SetRegionApproximation(actRegion, lagrangeMultPolyId_, lagrangeMultIntegId_);
      }
      if ( useLagrangeMultScal_ ) {
        lagrangeMultFct_1->AddEntityList( actSDList );
        lagrangeMultSpace_1->SetRegionApproximation(actRegion, lagrangeMultPolyId_, lagrangeMultIntegId_);
      }


      // --------------------------
      //  Set region approximation
      // --------------------------
      // We hardcode the Taylor-Hood spaces for the moment
      velSpace->SetRegionApproximation(actRegion, velPolyId_, velIntegId_);
      meanVelSpace->SetRegionApproximation(actRegion, velPolyId_, velIntegId_);
      presSpace->SetRegionApproximation(actRegion, presPolyId_, presIntegId_);

      PtrCoefFct density =
              materials_[actRegion]->GetScalCoefFnc(DENSITY, Global::REAL);

      // Create set of flow regions and map of density functions for surface integrators.
      flowRegions.insert(actRegion);
      oneFuncs[actRegion] = 
        CoefFunction::Generate( mp_, Global::REAL, lexical_cast<std::string>(1.0) );

      PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");

      // ====================================================================
      // stiffness integrators: conservation of mass
      //
      //  K_PV Integrator (lower off-diagonal integrator):
      //  \int_{\Omega_f} phi div(v') d\Omega
      // ===================================================================
      // This integrator is not created if the LinFlowPDE is coupled to the HeatPDE in symmetric form as only one of the
      // two symmteric integrators is needed and the counter part in the balance of momentum will be created instead of
      // this one
      if(!isHeatPDECoupled_) {

        BiLinearForm * stiffIntPV = NULL;
        if( dim_ == 2 ) {
          if( isaxi_ ) {
            stiffIntPV = new ABInt<>( new MultiIdOp<FeH1,2>(), 
                                    new DivOperatorAxi<FeH1>(), constOne, 1.0, updatedGeo_ );
          } else if(subType_ == "plane") {
            stiffIntPV = new ABInt<>( new MultiIdOp<FeH1,2>(), 
                                    new DivOperator<FeH1,2>(), constOne, 1.0, updatedGeo_ );
          }
        } else {
          stiffIntPV = new ABInt<>( new MultiIdOp<FeH1,3>(),
                                    new DivOperator<FeH1,3>(), constOne, 1.0, updatedGeo_ );
        }
        stiffIntPV->SetName("LinFlowStiffIntPV");
        BiLinFormContext *stiffContPV = NULL;
        stiffContPV = new BiLinFormContext(stiffIntPV, STIFFNESS );

        stiffContPV->SetEntities( actSDList, actSDList );
        stiffContPV->SetFeFunctions( presFct, velFct);
        assemble_->AddBiLinearForm( stiffContPV );
      }

      PtrCoefFct adiabaticExp = materials_[actRegion]->GetScalCoefFnc(
          FLUID_ADIABATIC_EXPONENT, Global::REAL);
      PtrCoefFct compressionModulus = materials_[actRegion]->GetScalCoefFnc( FLUID_BULK_MODULUS, Global::REAL );
      if ( isCompressible_ ) {
        // ====================================================================
        // C_PP: damping integrator, conservation of mass
        // add time derivative of density expressed by pressure according to
          // thermodynamic relation
        // In general form : adiabaticExp \(\rho * c^2)
        // ====================================================================
          PtrCoefFct fnc;
          if ( isHeatPDECoupled_ ) {
              fnc = CoefFunction::Generate( mp_, Global::REAL,
                              CoefXprBinOp(mp_, adiabaticExp, compressionModulus, CoefXpr::OP_DIV ) );
          }
          else {
              fnc = CoefFunction::Generate( mp_, Global::REAL,
                  CoefXprBinOp(mp_, constOne, compressionModulus, CoefXpr::OP_DIV ) );
          }

          BiLinearForm *dampIntpp = NULL;
          if( dim_ == 2 ) {
              dampIntpp = new BBInt<Double>(new IdentityOperator<FeH1,2,1,
                          Double>, fnc, 1.0, updatedGeo_ );
          } else {
              dampIntpp = new BBInt<Double>(new IdentityOperator<FeH1,3,1,
                          Double>, fnc, 1.0, updatedGeo_ );
          }
          dampIntpp->SetName("FlowDampIntPP");

          BiLinFormContext *dampContextpp = NULL;
          dampContextpp = new BiLinFormContext(dampIntpp, DAMPING );

          dampContextpp->SetEntities( actSDList, actSDList );
          dampContextpp->SetFeFunctions( presFct, presFct );
          assemble_->AddBiLinearForm( dampContextpp );
      }

      // ====================================================================
      // M_VV mass integrator, conservation of momentum
      // j omega rho v'.v (=inertia term)
      // ====================================================================
      BiLinearForm *dampIntvv = NULL;
      if( dim_ == 2 ) {
        dampIntvv = new BBInt<>(new IdentityOperator<FeH1,2,2>(),
                                density, balanceOfMomentumSign_, updatedGeo_ );
      } else {
        dampIntvv = new BBInt<>(new IdentityOperator<FeH1,3,3>(),
                                density, balanceOfMomentumSign_, updatedGeo_ );
      }
      dampIntvv->SetName("FlowDampIntVV");

      BiLinFormContext *dampContextvv = NULL;
      dampContextvv = new BiLinFormContext(dampIntvv, DAMPING );

      dampContextvv->SetEntities( actSDList, actSDList );
      dampContextvv->SetFeFunctions( velFct, velFct );
      assemble_->AddBiLinearForm( dampContextvv );

      // ====================================================================
      //  K_VP Integrator: conservation of momentum
      //  Div(v') p = Grad(v'):(pI) ... the pressure term of the stress tensor
      // TODO: write as scalar DivOp = more efficient
      // ====================================================================
      PtrCoefFct coeffKVP = CoefFunction::Generate(mp_,Global::REAL, "1.0");
      BiLinearForm * stiffIntVP = NULL;
      if( dim_ == 2 ) {
        if( isaxi_ ) {
          stiffIntVP = new ABInt<>( new DivOperatorAxi<FeH1>(),
                                      new MultiIdOp<FeH1,2>(), coeffKVP, -balanceOfMomentumSign_, updatedGeo_ );
        } else if(subType_ == "plane") {
          stiffIntVP = new ABInt<>( new DivOperator<FeH1,2>(),
                                    new MultiIdOp<FeH1,2>(), coeffKVP, -balanceOfMomentumSign_, updatedGeo_ );
        }
      } else {
        stiffIntVP = new ABInt<>( new DivOperator<FeH1,3>(),
                                  new MultiIdOp<FeH1,3>(), coeffKVP, -balanceOfMomentumSign_, updatedGeo_ );
      }
      stiffIntVP->SetName("LinFlowStiffIntVP");
      BiLinFormContext *stiffContVP = NULL;
      stiffContVP = new BiLinFormContext(stiffIntVP, STIFFNESS );

      stiffContVP->SetEntities( actSDList, actSDList );
      stiffContVP->SetFeFunctions( velFct, presFct );
      // In case the LinFLowPDE is coupled to the HeatPDE in a symmetric form the counterpart needs to be set, to also
      // create the other symmteric integrator in the balance of mass
      stiffContVP->SetCounterPart(isHeatPDECoupled_);
      assemble_->AddBiLinearForm( stiffContVP );


      // ====================================================================
      //  K_VV Integrator: conservation of momentum
      //  diagonal integrator - partially integrated
      //  2*mu*B(v'):B(v) ... the total strain-rate term
      // ====================================================================
      PtrCoefFct shearViscosity = materials_[actRegion]->GetScalCoefFnc(FLUID_DYNAMIC_VISCOSITY, Global::REAL); // mu
      PtrCoefFct shearViscosityDouble = CoefFunction::Generate( mp_,  Global::REAL,
          CoefXprBinOp(mp_, shearViscosity, CoefFunction::Generate( mp_, Global::REAL, "2"), CoefXpr::OP_MULT));
      PtrCoefFct coefZero = CoefFunction::Generate( mp_, Global::REAL, "0");
      BaseBDBInt * stiffIntLaplace = NULL;
      StdVector<PtrCoefFct> tensorComponents(dim_ == 2 ? (isaxi_ == true ? 16 : 9) : 36);
      tensorComponents.Init(coefZero);
      if( dim_ == 2 ) {
        if( isaxi_ ) {
          tensorComponents[0] = shearViscosityDouble;
          tensorComponents[5] = shearViscosityDouble;
          tensorComponents[10] = shearViscosity;
          tensorComponents[15] = shearViscosityDouble;
          PtrCoefFct coefBB = CoefFunction::Generate(mp_,Global::REAL,4,4,tensorComponents);
          stiffIntLaplace = new BDBInt<>( new StrainOperatorAxi<FeH1>(), coefBB, balanceOfMomentumSign_, updatedGeo_ );
        } else if( subType_ == "plane" ) {
          tensorComponents[0] = shearViscosityDouble;
          tensorComponents[4] = shearViscosityDouble;
          tensorComponents[8] = shearViscosity;
          PtrCoefFct coefBB = CoefFunction::Generate(mp_,Global::REAL,3,3,tensorComponents);
          stiffIntLaplace = new BDBInt<>( new StrainOperator2D<FeH1>(), coefBB, balanceOfMomentumSign_, updatedGeo_ );
        }
      } else {
        tensorComponents[0] = shearViscosityDouble;
        tensorComponents[7] = shearViscosityDouble;
        tensorComponents[14] = shearViscosityDouble;
        tensorComponents[21] = shearViscosity;
        tensorComponents[28] = shearViscosity;
        tensorComponents[35] = shearViscosity;
        PtrCoefFct coefBB = CoefFunction::Generate(mp_,Global::REAL,6,6,tensorComponents);
        stiffIntLaplace = new BDBInt<>( new StrainOperator3D<FeH1>(), coefBB, balanceOfMomentumSign_, updatedGeo_ );
      }
      stiffIntLaplace->SetName("LinFlowStiffIntViscous");
      BiLinFormContext *stiffContLaplace;
      stiffContLaplace = new BiLinFormContext(stiffIntLaplace, STIFFNESS );

      stiffContLaplace->SetEntities( actSDList, actSDList );
      stiffContLaplace->SetFeFunctions( velFct, velFct );
      assemble_->AddBiLinearForm( stiffContLaplace );

      // Important: Add bdb-integrator to global list, as we need them later
      // for calculation of postprocessing results
      // We can't add two integrators to a coefFunction, hence, by adding both to bdbInts_ we need to add name requests for the coefFunction
      bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffIntLaplace) );
      LOG_TRACE(linfluidmechpde) << "Add Lin BDB (strain part)" << std::endl;

      if ( isCompressible_ ) { // we need to subtract 2*mu/3 Div(v)I and add the bulk part lambda*Div(v)I
        PtrCoefFct bulkViscosity = materials_[actRegion]->GetScalCoefFnc(FLUID_BULK_VISCOSITY, Global::REAL);
        PtrCoefFct coefDivDiv = CoefFunction::Generate( mp_,  Global::REAL,
            CoefXprBinOp(mp_,
                bulkViscosity,
                CoefXprBinOp(mp_,shearViscosityDouble,CoefFunction::Generate( mp_, Global::REAL, "3"),CoefXpr::OP_DIV),
            CoefXpr::OP_SUB ));

        BaseBDBInt * stiffIntDivDiv = NULL;
        if( dim_ == 2 ) {
          if( isaxi_ ) {
            stiffIntDivDiv = new BBInt<>(new ScalarDivergenceOperatorAxi<FeH1>(), coefDivDiv, balanceOfMomentumSign_, updatedGeo_ );
          } else if(subType_ == "plane") {
            stiffIntDivDiv = new BBInt<>(new ScalarDivergenceOperator<FeH1,2,Double>(), coefDivDiv, balanceOfMomentumSign_, updatedGeo_);
          }         
        } else {
          stiffIntDivDiv = new BBInt<>(new ScalarDivergenceOperator<FeH1,3,Double>(), coefDivDiv, balanceOfMomentumSign_, updatedGeo_);
        }
        stiffIntDivDiv->SetName("LinFlowStiffIntBulkViscous");
        BiLinFormContext *stiffContDivDiv;
        stiffContDivDiv = new BiLinFormContext(stiffIntDivDiv, STIFFNESS );

        stiffContDivDiv->SetEntities( actSDList, actSDList );
        stiffContDivDiv->SetFeFunctions( velFct, velFct );
        assemble_->AddBiLinearForm( stiffContDivDiv );

        // Important: Add bdb-integrator to global list, as we need them later
        // for calculation of postprocessing results
        // We can't add two integrators to a coefFunction, hence, by adding both to bdbInts_ we need to add name requests for the coefFunction
        bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffIntDivDiv) );
        LOG_TRACE(linfluidmechpde) << "Add Lin BDB (compressible part)" << std::endl;
      }


      //======================================================================
      // CHECK FOR BASE BACKGROUND FLOW
      //======================================================================
      std::string flowId = curRegNode->Get("flowId")->As<std::string>();
       
      if(flowId == "") {
        
          // Check if no meanFluidMechVelocity is defined if no background flow exists 
          if(myParam_->Get("storeResults")->HasByVal("nodeResult","type","meanFluidMechVelocity")) {
            EXCEPTION("Result value meanFluidMechVelocity defined but no background flow is defined!");
          }
      }

      // This makes sure that the meanFluidMechVelocity is only defined if a background flow is defined in the region 
      else {

        if( isaxi_ ) {
          EXCEPTION("Axi formulation not checked for background flow, please check and verify!");
        }

        //Check if result value meanFluidMechVelocity is defined
        if(myParam_->Get("storeResults")->HasByVal("nodeResult","type","meanFluidMechVelocity")) {
          
          ParamNodeList regionListAll = myParam_->Get("regionList")->GetList("region");

          //Check for allRegions in result value meanFluidMechVelocity
          if(myParam_->Get("storeResults")->GetByVal("nodeResult","type","meanFluidMechVelocity")->Has("allRegions")) { 
            
            //iterate over allRegions 
            for (boost::shared_ptr<ParamNode> region : regionListAll) {
              
              //Check if a background flow is defined
              if(!region->Has("flowId")) {
                 EXCEPTION("Result value meanFluidMechVelocity defined an a region that has no background flow defined!");
              }
            }
          }

          //Check for specific defined regions
          else if(myParam_->Get("storeResults")->GetByVal("nodeResult","type","meanFluidMechVelocity")->Has("regionList")) { 
            
            ParamNodeList regionListMeanFlow = myParam_->Get("storeResults")->GetByVal("nodeResult","type","meanFluidMechVelocity")->Get("regionList")->GetList("region");
            
            //iterate over all regions meanFluidMechVelocity is calculated on
            for (boost::shared_ptr<ParamNode> regionMean : regionListMeanFlow) {

              //Check if the meanFluidMechVelocity region is defined in the region list
              if(myParam_->Get("regionList")->HasByVal("region","name",regionMean->Get("name")->As<std::string>())) {
                
                //Check if the meanFluidMechVelocity region has a flowId defined in the region list
                if(!myParam_->Get("regionList")->GetByVal("region","name",regionMean->Get("name")->As<std::string>())->Has("flowId")) {
                  EXCEPTION("Result value meanFluidMechVelocity defined on region that has no background flow defined!");
                }
              } else {
                EXCEPTION("Unknown region defined result value meanFluidMechVelocity");
              }
            }
          }
        } 


        // Get result info object for flow
        shared_ptr<ResultInfo> velocityInfo = GetResultInfo(MEAN_FLUIDMECH_VELOCITY);
        // Read coefficient flow coefficient function for this region
        bool coefUpdateGeo = false;
        PtrCoefFct velocityCoef;
        std::set<UInt> definedDofs;

        //read from XML scheme and add the region information
        PtrParamNode velocityNode = myParam_->Get("flowList")->GetByVal("velocity","name", flowId.c_str());
        ReadUserFieldValues( actSDList, velocityNode, velocityInfo->dofNames, velocityInfo->entryType, isComplex_, velocityCoef, definedDofs, coefUpdateGeo );
        meanVelocityCoef_->AddRegion( actRegion, velocityCoef );

        // add stuff to the mean flow FE function
        meanVelFct->AddEntityList( actSDList );
        meanVelFct->AddExternalDataSource( velocityCoef, actSDList );

        // first convective term in balance of momentum
        // v' . ( rho0 v0.Grad(v) )
        BiLinearForm *convectiveVv = NULL;
        if( dim_ == 2 ) {
          if(isComplex_) 
          {
            convectiveVv = new ABInt<Complex>( new IdentityOperator<FeH1,2,2>(), new ConvectiveOperator<FeH1,2,2,Complex>(), density, balanceOfMomentumSign_, coefUpdateGeo );
          }
          else
          {
            convectiveVv = new ABInt<Double>( new IdentityOperator<FeH1,2,2>(), new ConvectiveOperator<FeH1,2,2>(), density, balanceOfMomentumSign_, coefUpdateGeo );
          }
        } else {
          if(isComplex_) 
          {
            convectiveVv = new ABInt<Complex>( new IdentityOperator<FeH1,3,3>(), new ConvectiveOperator<FeH1,3,3,Complex>(), density, balanceOfMomentumSign_, coefUpdateGeo );
          }
          else
          {
            convectiveVv = new ABInt<Double>( new IdentityOperator<FeH1,3,3>(), new ConvectiveOperator<FeH1,3,3>(), density, balanceOfMomentumSign_, coefUpdateGeo );
          }
        }
        // here we set the velocity vector used in the ConvectiveOperator above
        convectiveVv->SetBCoefFunctionOpB( meanVelocityCoef_ );
        convectiveVv->SetName("LinFlowStiffIntConvectiveVv");

        BiLinFormContext *convectiveContextVv = NULL;
        convectiveContextVv = new BiLinFormContext(convectiveVv, STIFFNESS );

        convectiveContextVv->SetEntities( actSDList, actSDList );
        convectiveContextVv->SetFeFunctions( velFct, velFct );
        assemble_->AddBiLinearForm( convectiveContextVv );

        PtrCoefFct rho0GradV0,coeffV0inGradV0;
        BaseBOperator* gradOp;
        if(isComplex_) {
          if (dim_==2){
            gradOp = new GradientOperator<FeH1,2, 1, Complex>();
            rho0GradV0.reset( new CoefFunctionMeanFlowConvection<Complex,2>( density, gradOp, meanVelFct ) );
            coeffV0inGradV0.reset( new CoefFunctionMeanFlowConvection<Complex,2>(gradOp, meanVelFct) );
          } else {
            gradOp = new GradientOperator<FeH1,3, 1, Complex>();
            rho0GradV0.reset( new CoefFunctionMeanFlowConvection<Complex,3>( density, gradOp, meanVelFct ) );
            coeffV0inGradV0.reset( new CoefFunctionMeanFlowConvection<Complex,3>(gradOp, meanVelFct) );
          }
        } else {
          if (dim_==2){
            gradOp = new GradientOperator<FeH1,2, 1>();
            rho0GradV0.reset( new CoefFunctionMeanFlowConvection<Double,2>( density, gradOp, meanVelFct ) );
            coeffV0inGradV0.reset( new CoefFunctionMeanFlowConvection<Double,2>(gradOp, meanVelFct) );
          } else {
            gradOp = new GradientOperator<FeH1,3, 1>();
            rho0GradV0.reset( new CoefFunctionMeanFlowConvection<Double,3>( density, gradOp, meanVelFct ) );
            coeffV0inGradV0.reset( new CoefFunctionMeanFlowConvection<Double,3>(gradOp, meanVelFct) );
          }
        }

        if(enableC2_)
        {
          // Second convective term. Derivative tensor of mean flow field is a
          // factor computed by a CoefFunction.
          // v' . ( rho0 Grad(v0) . v )
          BaseBOperator* bOpId;
          if( dim_ == 2 ) {
            bOpId = new IdentityOperator<FeH1,2,2>();
          }
          else {
            bOpId = new IdentityOperator<FeH1,3,3>();
          }
          //now create the integrators
          BiLinearForm *convectivevV = NULL;
          PtrCoefFct coeffConvec;
          if(isComplex_) {
            convectivevV = new BDBInt<Complex,Complex>( bOpId, rho0GradV0, balanceOfMomentumSign_ );
          } else {
            convectivevV = new BDBInt<Double,Double>( bOpId, rho0GradV0, balanceOfMomentumSign_ );
          }
          convectivevV->SetName("LinFlowStiffIntConvectivevV");
          // assign to context
          BiLinFormContext *convectiveContextvV = new BiLinFormContext(convectivevV, STIFFNESS );
          convectiveContextvV->SetEntities( actSDList, actSDList );
          convectiveContextvV->SetFeFunctions( velFct, velFct );
          assemble_->AddBiLinearForm( convectiveContextvV );
        }
        if ( isCompressible_ ){
          // convective term for constant base density and pressure
          // p' v0/K . Grad(p)
          PtrCoefFct coefPp = CoefFunction::Generate( mp_,  Global::REAL,
                      CoefXprBinOp(mp_,
                          constOne,
                          compressionModulus,
                      CoefXpr::OP_DIV ));          
          BiLinearForm *convectivePp = NULL;
          if(isComplex_) {
            if( dim_ == 2 ) { // 2D
              convectivePp = new ABInt<Complex>( new IdentityOperator<FeH1,2>(), new ConvectiveOperator<FeH1,2,1>(), coefPp, 1.0, coefUpdateGeo );
            } else { // 3D
              convectivePp = new ABInt<Complex>( new IdentityOperator<FeH1,3>(), new ConvectiveOperator<FeH1,3,1>(), coefPp, 1.0, coefUpdateGeo );
            }
          }
          else {
             if (dim_==2) {
               convectivePp = new ABInt<Double>( new IdentityOperator<FeH1,2>(), new ConvectiveOperator<FeH1,2,1,Complex>(), coefPp, 1.0, coefUpdateGeo );
             } else {
               convectivePp = new ABInt<Double>( new IdentityOperator<FeH1,3>(), new ConvectiveOperator<FeH1,3,1,Complex>(), coefPp, 1.0, coefUpdateGeo );
             }
          }
          convectivePp->SetName("LinFlowStiffIntConvectivePp");
          //SetCoefFunction : meanVelocity for ConvectiveOperator
          convectivePp->SetBCoefFunctionOpB(meanVelocityCoef_);

          BiLinFormContext *convectiveContextPp = new BiLinFormContext(convectivePp, STIFFNESS );
          convectiveContextPp->SetEntities( actSDList, actSDList );
          convectiveContextPp->SetFeFunctions( presFct, presFct );
          assemble_->AddBiLinearForm( convectiveContextPp );

          // 3rd convective term for momentum equation
          // v' . rho0/(gamma p0) v0.Grad(v0) p
          if(enableC3_){
          // scalar factor for the term: rho0/K
          PtrCoefFct rho0OverGammaP0 = CoefFunction::Generate( mp_,  Global::REAL,
                                          CoefXprBinOp(mp_,
                                              density,
                                              compressionModulus,
                                          CoefXpr::OP_DIV ));
          // multiply vector v0 . Grad(v0) with scalar factor
          PtrCoefFct rho0OverGammaP0V0inGradV0 = CoefFunction::Generate( mp_,  Global::COMPLEX, CoefXprVecScalOp(mp_,coeffV0inGradV0,rho0OverGammaP0,CoefXpr::OP_MULT));

          // finally convert to a D-tensor for the ADB integrator
          PtrCoefFct coefPv;
          if (isComplex_) {
            coefPv.reset( new CoefFunctionVectorToTensor<Complex>(rho0OverGammaP0V0inGradV0,true) );
          } else {
            coefPv.reset( new CoefFunctionVectorToTensor<Double>(rho0OverGammaP0V0inGradV0,true) );
          }

          // This would be the fast way, but fails with index error - no idea why
          //PtrCoefFct coefPv = CoefFunction::Generate( mp_,  Global::COMPLEX, CoefXprTensScalOp(mp_,tensor,rho0OverKappaP0,CoefXpr::OP_MULT));

          BiLinearForm *convectiveVp = NULL;
          if(isComplex_) {
            if(dim_==2){
              convectiveVp = new ADBInt<Complex>( new MultiIdOp<FeH1,1,2,Complex>(), new MultiIdOp<FeH1,2,1,Complex>(), coefPv, balanceOfMomentumSign_, coefUpdateGeo );
            } else {
              convectiveVp = new ADBInt<Complex>( new MultiIdOp<FeH1,1,3,Complex>(), new MultiIdOp<FeH1,3,1,Complex>(), coefPv, balanceOfMomentumSign_, coefUpdateGeo );
            }
          } else {
            if(dim_==2){
              convectiveVp = new ADBInt<Double>( new MultiIdOp<FeH1,1,3,Double>(), new MultiIdOp<FeH1,3,1,Double>(), coefPv, balanceOfMomentumSign_, coefUpdateGeo );
            } else {
              convectiveVp = new ADBInt<Double>( new MultiIdOp<FeH1,1,3,Double>(), new MultiIdOp<FeH1,3,1,Double>(), coefPv, balanceOfMomentumSign_, coefUpdateGeo );
            }
          }
          convectiveVp->SetName("LinFlowStiffIntConvectiveVp");
          BiLinFormContext *convectiveContextVp = new BiLinFormContext(convectiveVp, STIFFNESS );
          convectiveContextVp->SetEntities( actSDList, actSDList );
          convectiveContextVp->SetFeFunctions( velFct, presFct );
          assemble_->AddBiLinearForm( convectiveContextVp );
          }
        }
      } // is flow

      //======================================================================
      // ALE-section (check for grid-velocity)
      //======================================================================
      
      std::string movingMeshId = curRegNode->Get("movingMeshId")->As<std::string>();
      if(movingMeshId != ""){
        if(flowId != ""){
          EXCEPTION("Moving mesh with background flow is currently not supported!")
        }

        //if( useMovingMeshTerms_ ){
        // Manually set this bool since the Mech-PDE sets the updatetGeo_ bool to false
        // by default and this term only makes sense when the geometry is updated
        bool updateGeoLF;
        updateGeoLF = true; // this is just a dummy: if the mechPDE is not computed on the updated geometry, this will get set to false

        // Get result info object for flow
        shared_ptr<ResultInfo> movingMeshInfo;
        movingMeshInfo = GetResultInfo(FLUIDMECH_MESH_VELOCITY);
        // Read ALE coefficient function for this region
        PtrCoefFct regionMovingMesh;
        std::set<UInt> definedDofs;

        //Add the region information
        PtrParamNode movingMeshNode =
          myParam_->Get("movingMeshList")->GetByVal("movingMesh",
                                              "name",
                                              movingMeshId.c_str());

        //TODO Test ReadRHSExcitation and let entList be defined by it (instead of actSDList which are elements)
        StdVector<shared_ptr<EntityList> > ent;
        StdVector<PtrCoefFct > tCoef;
        //ReadRhsExcitation("velocity", movingMeshInfo->dofNames, ResultInfo::VECTOR, isComplex_, ent, tCoef, updateGeoLF, movingMeshNode);
        ReadUserFieldValues( actSDList, movingMeshNode, movingMeshInfo->dofNames,
                              movingMeshInfo->entryType, isComplex_, regionMovingMesh,
                              definedDofs, updateGeoLF );
        gridVelCoef_->AddRegion( actRegion, regionMovingMesh );


        if ( enableGridVelC1_ ) {
          if ( isCompressible_ ) {
            //now create the integrators
            // ====================================================================
            // K_PP (MG): stiffness integrator, moving grid (grid velocity v_g)
            // \int_{Omega_f} v_g/K \phi \grad{p_a} d\Omega
            // ====================================================================
            BiLinearForm *stiffIntMovingGridPP = NULL;

            if (isHeatPDECoupled_) {
              EXCEPTION("Heat coupling with moving grid is not implemented yet")
            }
            else {
              PtrCoefFct fnc;
              fnc = CoefFunction::Generate( mp_, Global::REAL,
                                CoefXprBinOp(mp_, constOne, compressionModulus, CoefXpr::OP_DIV ) );

              if( dim_ == 2 ) {
                if(isComplex_)
                {
                  stiffIntMovingGridPP = new ABInt<Complex>( new MultiIdOp<FeH1,2>(),
                                                      new ScaledGradientOperator<FeH1,2,Complex>(),
                                                      fnc, -1.0, updatedGeo_ );
                }
                else
                {
                  stiffIntMovingGridPP = new ABInt<Double>( new MultiIdOp<FeH1,2>(),
                                                    new ScaledGradientOperator<FeH1,2>(),
                                                    fnc, -1.0, updatedGeo_ );
                }
              } else {
                if(isComplex_)
                {
                  stiffIntMovingGridPP = new ABInt<Complex>( new MultiIdOp<FeH1,3>(),
                                                      new ScaledGradientOperator<FeH1,3,Complex>(),
                                                      fnc, -1.0, updatedGeo_ );
                }
                else
                {
                  stiffIntMovingGridPP = new ABInt<Double>( new MultiIdOp<FeH1,3>(),
                                                    new ScaledGradientOperator<FeH1,3>(),
                                                    fnc, -1.0, updatedGeo_ );
                }
              }
            }

            stiffIntMovingGridPP->SetBCoefFunctionOpB(gridVelCoef_);
            stiffIntMovingGridPP->SetSolDependent(true); // TODO: does this trigger the update of the coefFunction?

            stiffIntMovingGridPP->SetName("LinFlowStiffIntMovingGridPP");

            BiLinFormContext *stiffIntMovingGridContextPP = NULL;
            stiffIntMovingGridContextPP = new BiLinFormContext(stiffIntMovingGridPP, STIFFNESS );

            stiffIntMovingGridContextPP->SetEntities( actSDList, actSDList );
            stiffIntMovingGridContextPP->SetFeFunctions( presFct, presFct );
            assemble_->AddBiLinearForm( stiffIntMovingGridContextPP );
          } else {
            EXCEPTION("GridVel C1 term activated but formulation is incompressible!")
          }
        }

        if ( enableGridVelC2_ ) {
          // ====================================================================
          // K_VV (MG): stiffness integrator, moving grid (grid velocity v_g)
          // \int_{Omega_f} \rho v' (v_g \cdot \grad{v}) d\Omega
          // ====================================================================
          BiLinearForm *stiffIntMovingGridVV = NULL;

          if( dim_ == 2 ) {
            // Axi-symmetric version is identical to plane when setting v_g,\phi = 0 and assuming no changes in \phi direction
            if(isComplex_)
            {
              stiffIntMovingGridVV = new ABInt<Complex>( new IdentityOperator<FeH1,2,2>(),
                                                  new ConvectiveOperator<FeH1,2,2,Complex>(),
                                                  density, -balanceOfMomentumSign_, updatedGeo_ );
            }
            else
            {
              stiffIntMovingGridVV = new ABInt<Double>( new IdentityOperator<FeH1,2,2>(),
                                                new ConvectiveOperator<FeH1,2,2>(),
                                                density, -balanceOfMomentumSign_, updatedGeo_ );
            }
          } else {
            if(isComplex_)
            {
              stiffIntMovingGridVV = new ABInt<Complex>( new IdentityOperator<FeH1,3,3>(),
                                                  new ConvectiveOperator<FeH1,3,3,Complex>(),
                                                  density, -balanceOfMomentumSign_, updatedGeo_ );
            }
            else
            {
              stiffIntMovingGridVV = new ABInt<Double>( new IdentityOperator<FeH1,3,3>(),
                                                new ConvectiveOperator<FeH1,3,3>(),
                                                density, -balanceOfMomentumSign_, updatedGeo_ );
            }
          }

          stiffIntMovingGridVV->SetBCoefFunctionOpB(gridVelCoef_);
          stiffIntMovingGridVV->SetSolDependent(true);

          stiffIntMovingGridVV->SetName("LinFlowStiffIntMovingGridVV");

          BiLinFormContext *stiffIntMovingGridContextVV = NULL;
          stiffIntMovingGridContextVV = new BiLinFormContext(stiffIntMovingGridVV, STIFFNESS );

          stiffIntMovingGridContextVV->SetEntities( actSDList, actSDList );
          stiffIntMovingGridContextVV->SetFeFunctions( velFct, velFct );
          assemble_->AddBiLinearForm( stiffIntMovingGridContextVV );
          //}
        }
      }
    }
    
    // ====================================================================
    // When pressure is prescribed, we need to add the corresponding
    // surface integrator, since we have integrated by parts the gradient
    // of the pressure in the momentum conservation equation
    // ====================================================================
    PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    if( bcNode ) {
        ParamNodeList surfPresNodes = bcNode->GetList( "pressure" );
        //std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );

        for( UInt i = 0; i < surfPresNodes.GetSize(); i++ ) {
            std::string regionName = surfPresNodes[i]->Get("name")->As<std::string>();
            shared_ptr<EntityList> actSDList = ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );

            velFct->AddEntityList( actSDList );
            // --------------------------------------------------------------------
            //  VERSION 2: K_VP Integrator
            //  (upper off-diagonal integrators - partially integrated, surface)
            // --------------------------------------------------------------------
            BiLinearForm * stiffIntVPSurf = NULL;
            if( dim_ == 2 ) {
                stiffIntVPSurf = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                     new IdentityOperatorNormal<FeH1,2>(),
                                     oneFuncs, balanceOfMomentumSign_, flowRegions, updatedGeo_);
            } else {
                stiffIntVPSurf = new SurfaceABInt<>(new IdentityOperator<FeH1,3,3>(),
                                     new IdentityOperatorNormal<FeH1,3>(),
                                     oneFuncs, balanceOfMomentumSign_, flowRegions, updatedGeo_);
            }
            stiffIntVPSurf->SetName("LinFlowStiffIntVPSurf");
            BiLinFormContext *stiffContVP = NULL;
            stiffContVP = new BiLinFormContext(stiffIntVPSurf, STIFFNESS );

            stiffContVP->SetEntities( actSDList, actSDList );
            stiffContVP->SetFeFunctions( velFct, presFct );
            // In case the LinFLowPDE is coupled to the HeatPDE in a symmetric form the counterpart needs to be set, to also
            // create symmetric counterpart
            stiffContVP->SetCounterPart(isHeatPDECoupled_);
            assemble_->AddBiLinearForm( stiffContVP );
        }
    }

    // Since we manage the mean flow FeFunction, we also have to finalize here.
    // c.f. SinglePDE::InitStage3
    meanVelSpace->Finalize();
    meanVelSpace->PreCalcShapeFncs();
    meanVelFct->Finalize();
  }
  
  void LinFlowPDE::DefineRhsLoadIntegrators(PtrParamNode input)
  {
    // Get FESpace and FeFunction of fluid velocity
    shared_ptr<BaseFeFunction> velFct = feFunctions_[FLUIDMECH_VELOCITY];
    shared_ptr<FeSpace> velSpace = velFct->GetFeSpace();

    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    LinearForm * lin = NULL;
    StdVector<std::string> dispDofNames = velFct->GetResultInfo()->dofNames;
    bool coefUpdateGeo = false;
    // ==================
    //  SURFACE TRACTION
    // ==================
    //      LOG_DBG(mechpde) << "Reading surface tractions";

    ReadRhsExcitation("traction", dispDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo, input);
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Surface traction must be defined on elements")
      }
      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim != (dim_-1) ) {
        EXCEPTION("Surface traction can only be defined on surface elements");
      }
      if( dim_ == 2) {
        if(isComplex_) {
          lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,2,2>(),
              Complex(balanceOfMomentumSign_), coef[i], coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,2,2>(),
              balanceOfMomentumSign_, coef[i],coefUpdateGeo);
        }
      } else  {
        if(isComplex_) {
          lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,3,3>(),
              Complex(balanceOfMomentumSign_), coef[i], coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,3,3>(),
              balanceOfMomentumSign_, coef[i], coefUpdateGeo);
        }
      }
      lin->SetName("TractionIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(velFct);
      assemble_->AddLinearForm(ctx);
      velFct->AddEntityList(ent[i]);
    }

    // ===============
    //  Normal Traction
    // ===============
    //    LOG_DBG(mechpde) << "Reading mechanical pressure";
    StdVector<std::string> empty;
    ReadRhsExcitation("normalTraction", empty, ResultInfo::SCALAR, isComplex_, ent, coef, coefUpdateGeo, input);
    std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );

    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Mechanical pressure must be defined on elements")
      }

      // to have correct normal direction
      const Double tracFac = -1.0;

      if( dim_ == 2) {
        if(isComplex_) {
          lin = new BUIntegrator<Complex, true> ( new IdentityOperatorNormalTrans<FeH1,2>(),
              Complex(balanceOfMomentumSign_*tracFac), coef[i],
              volRegions, coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double,true> ( new IdentityOperatorNormalTrans<FeH1,2>(),
              balanceOfMomentumSign_*tracFac, coef[i], volRegions,
              coefUpdateGeo);
        }
      } else  {
        if(isComplex_) {
          lin = new BUIntegrator<Complex, true> ( new IdentityOperatorNormalTrans<FeH1,3>(),
              Complex(balanceOfMomentumSign_*tracFac), coef[i],
              volRegions, coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double, true> ( new IdentityOperatorNormalTrans<FeH1,3>(),
              balanceOfMomentumSign_*tracFac, coef[i],
              volRegions, coefUpdateGeo);
        }
      }
      lin->SetName("PressureInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(velFct);
      assemble_->AddLinearForm(ctx);
      velFct->AddEntityList(ent[i]);
    } // for


    //========================================================================================
    // Lagrange multiplier based weak Dirichlet BCs
    //========================================================================================

    // Get FESpace and FeFunction of the lagrange multiplier
    shared_ptr<BaseFeFunction> lagrangeMultFct = NULL;
    shared_ptr<FeSpace> lagrangeMultSpace = NULL;
    shared_ptr<BaseFeFunction> lagrangeMultFct_1 = NULL;
    shared_ptr<FeSpace> lagrangeMultSpace_1 = NULL;

    // we only assemble the integrator if there are some BCs specified and if mesh movement is present
    PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    if( bcNode ) {
      StdVector<shared_ptr<EntityList> > ent;
      StdVector<PtrCoefFct > kCoef;
      StdVector<std::string> volumeRegions;

      // Define FE space etc.
      if ( useLagrangeMultScal_ ) {
        lagrangeMultFct_1 = feFunctions_[LAGRANGE_MULT_1];
        lagrangeMultSpace_1 = lagrangeMultFct_1->GetFeSpace();
      }

      ReadRhsExcitation("scalarVelocityConstraint", feFunctions_[FLUIDMECH_VELOCITY]->GetResultInfo()->dofNames,
          ResultInfo::SCALAR, isComplex_, ent, kCoef, updatedGeo_, volumeRegions);

      for( UInt i = 0; i < ent.GetSize(); ++i ) {
        // get the volume region for defining the correct normal direction
        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volumeRegions[i]);
        std::set<RegionIdType> volRegion;
        volRegion.insert(aRegion);
        // check type of entitylist
        if (ent[i]->GetType() != EntityList::SURF_ELEM_LIST) {
          EXCEPTION("Velocity constraint BC must be defined on (surface) elements")
        }
        std::string regionName = ent[i]->GetName();
        shared_ptr<EntityList> actSDList = ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );

        // check for grid velocity, if the coefFunction is not present, skip the definition of the integrator
        std::string volRegionName = ptGrid_->GetRegion().ToString(aRegion);
        PtrParamNode regionNode = myParam_->Get("regionList")->GetByVal("region","name",volRegionName.c_str());
        std::string movingMeshID;
        regionNode->GetValue( "movingMeshId", movingMeshID, ParamNode::PASS );

        if( !movingMeshID.empty() ) {
          // Get polynomial and integration order for Lagrange multiplier
          lagrangeMultSpace_1->SetRegionApproximation(aRegion, lagrangeMultPolyId_, lagrangeMultIntegId_);

          // check for noPenetration BC
          PtrParamNode scalVelConstNode = bcNode->GetByVal("scalarVelocityConstraint","name", regionName);
          ParamNodeList scalVelConstNodeList = scalVelConstNode->GetChildren();

          for (UInt ii = 0; ii < scalVelConstNodeList.GetSize(); ii++){

            std::string nodeName = scalVelConstNodeList[ii]->GetName();
        
            if( nodeName=="noPenetration" ) {
              // build the RHS
              // constraint equation: \int_\Gamma_D \lambda^\star n \cdot  v_g
              // since this is dependent on the grid velocity we will always have a geometry update
              PtrCoefFct gridVelCoefRegion = gridVelCoef_->GetRegionCoef(aRegion);

              if( dim_ == 2) {
                if(isComplex_) {
                  lin = new BUIntegrator<Complex, true> ( new IdentityOperatorNormal<FeH1,2,1>(),
                      Complex(balanceOfMomentumSign_,0), gridVelCoefRegion, volRegion, true, true);
                } else {
                  lin = new BUIntegrator<Double,true> ( new IdentityOperatorNormal<FeH1,2,1>(),
                      balanceOfMomentumSign_, gridVelCoefRegion, volRegion, true, true);
                }
              } else  {
                if(isComplex_) {
                  lin = new BUIntegrator<Complex, true> ( new IdentityOperatorNormal<FeH1,3,1>(),
                      Complex(balanceOfMomentumSign_,0), gridVelCoefRegion, volRegion, true, true);
                } else {
                  lin = new BUIntegrator<Double, true> ( new IdentityOperatorNormal<FeH1,3,1>(),
                      balanceOfMomentumSign_, gridVelCoefRegion, volRegion, true, true);
                }
              }
              lin->SetName("noPenetrationInt");
              LinearFormContext *ctx = new LinearFormContext( lin );
              ctx->SetEntities( ent[i] );
              ctx->SetFeFunction(lagrangeMultFct_1);
              assemble_->AddLinearForm(ctx);
              lagrangeMultFct_1->AddEntityList(ent[i]);
            }
          }
        }
      }

      // Define FE space etc.
      if ( useLagrangeMultVec_ ) {
        lagrangeMultFct = feFunctions_[LAGRANGE_MULT];
        lagrangeMultSpace = lagrangeMultFct->GetFeSpace();
      }
      
      ReadRhsExcitation("velocityConstraint", feFunctions_[FLUIDMECH_VELOCITY]->GetResultInfo()->dofNames,
          ResultInfo::SCALAR, isComplex_, ent, kCoef, updatedGeo_, volumeRegions);

      for( UInt i = 0; i < ent.GetSize(); ++i ) {
        // get the volume region for defining the correct normal direction
        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volumeRegions[i]);
        std::set<RegionIdType> volRegion;
        volRegion.insert(aRegion);
        // check type of entitylist
        if (ent[i]->GetType() != EntityList::SURF_ELEM_LIST) {
          EXCEPTION("Velocity constraint BC must be defined on (surface) elements")
        }
        std::string regionName = ent[i]->GetName();
        shared_ptr<EntityList> actSDList = ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );

        // check for grid velocity, if the coefFunction is not present, skip the definition of the integrator
        std::string volRegionName = ptGrid_->GetRegion().ToString(aRegion);
        PtrParamNode regionNode = myParam_->Get("regionList")->GetByVal("region","name",volRegionName.c_str());
        std::string movingMeshID;
        regionNode->GetValue( "movingMeshId", movingMeshID, ParamNode::PASS );

        if( !movingMeshID.empty() ) {
          // Get polynomial and integration order for Lagrange multiplier
          lagrangeMultSpace->SetRegionApproximation(aRegion, lagrangeMultPolyId_, lagrangeMultIntegId_);

          // check for noSlip and Maxwell slip BC
          PtrParamNode velConstNode = bcNode->GetByVal("velocityConstraint","name", regionName);
          ParamNodeList velConstNodeList = velConstNode->GetChildren();

          for (UInt ii = 0; ii < velConstNodeList.GetSize(); ii++){

            std::string nodeName = velConstNodeList[ii]->GetName();
        
            if( nodeName=="MaxwellFirstOrderSlip" || nodeName=="noSlip" ) {
              // build the RHS
              // constraint equation: \int_\Gamma_D \lambda^\star \cdot  v_g
              // since this is dependent on the grid velocity we will always have a geometry update
              PtrCoefFct gridVelCoefRegion = gridVelCoef_->GetRegionCoef(aRegion);

              if( dim_ == 2) {
                if(isComplex_) {
                  lin = new BUIntegrator<Complex, true> ( new IdentityOperator<FeH1,2,2>(),
                      Complex(balanceOfMomentumSign_,0), gridVelCoefRegion, volRegion, true, true);
                } else {
                  lin = new BUIntegrator<Double,true> ( new IdentityOperator<FeH1,2,2>(),
                      balanceOfMomentumSign_, gridVelCoefRegion, volRegion, true, true);
                }
              } else  {
                if(isComplex_) {
                  lin = new BUIntegrator<Complex, true> ( new IdentityOperator<FeH1,3,3>(),
                      Complex(balanceOfMomentumSign_,0), gridVelCoefRegion, volRegion, true, true);
                } else {
                  lin = new BUIntegrator<Double, true> ( new IdentityOperator<FeH1,3,3>(),
                      balanceOfMomentumSign_, gridVelCoefRegion, volRegion, true, true);
                }
              }
              lin->SetName("SlipInt");
              LinearFormContext *ctx = new LinearFormContext( lin );
              ctx->SetEntities( ent[i] );
              ctx->SetFeFunction(lagrangeMultFct);
              assemble_->AddLinearForm(ctx);
              lagrangeMultFct->AddEntityList(ent[i]);
            }
          }
        }
      }
    } // LM based weak Dirichlet BC

  }

  void LinFlowPDE::DefineNcIntegrators() {
//	if ( complexFluidFormulation_ )
//		EXCEPTION("Complex fluid and NC-interfaces currently not allowed");

    StdVector< NcInterfaceInfo >::iterator ncIt = ncInterfaces_.Begin(),
                                           endIt = ncInterfaces_.End();
    for ( ; ncIt != endIt; ++ncIt ) {
      switch (ncIt->type) {
      case NC_MORTAR:
        if (dim_ == 2)
          EXCEPTION("Mortar ncInterface not implemented")
        else
          EXCEPTION("Mortar ncInterface not implemented")
        break;
      case NC_NITSCHE:
        if (dim_ == 2)
          DefineNitscheCoupling<2,2,1>(*ncIt );
        else
          DefineNitscheCoupling<3,3,1>(*ncIt );
        break;
      default:
        EXCEPTION("Unknown type of ncInterface");
        break;
      }
    }
  }

  template<UInt DIM, UInt D_DOF, UInt P_DOF>
  void LinFlowPDE::DefineNitscheCoupling( NcInterfaceInfo &iface,
                                         shared_ptr<CoefFunctionMulti> additionalCoef )
  {
    // TODO check and remove this if everything is ok
    if( isMaterialComplex_ ) {
      EXCEPTION("LinFlow ncInterfaces not tested for complex values. Please verify if its working and remove this exception.")
    }

    shared_ptr<BaseNcInterface> ncIf = ptGrid_->GetNcInterface(iface.interfaceId);
    MortarInterface * nitscheIf = dynamic_cast<MortarInterface*>(ncIf.get());
    assert(nitscheIf);
    
    //in case of Nitsche coupling edge/face information is required
    this->ptGrid_->MapEdges();
    this->ptGrid_->MapFaces();

    // currently we have a moving formulation only for acoustics
    updatedGeo_ = updatedGeo_ || ncIf->NeedsUpdate(); // TODO jens: isn't is this too late?
    bool isMoving = updatedGeo_;

    // create new entity list
    shared_ptr<ElemList> actSDList = ncIf->GetElemList();

    // we set here the penalty factor
    Double beta = iface.nitscheFactor;

    // material parameter and adaption of penalty term
    // check if both interfaces use the same material of type CONSTANT, othewise throw an exception
    PtrCoefFct shearViscosityMaster = materials_[nitscheIf->GetMasterVolRegion()]->GetScalCoefFnc(FLUID_DYNAMIC_VISCOSITY, Global::REAL); // mu
    PtrCoefFct shearViscositySlave = materials_[nitscheIf->GetSlaveVolRegion()]->GetScalCoefFnc(FLUID_DYNAMIC_VISCOSITY, Global::REAL); // mu
    
    PtrCoefFct bulkViscosityMaster = materials_[nitscheIf->GetMasterVolRegion()]->GetScalCoefFnc(FLUID_BULK_VISCOSITY, Global::REAL);
    PtrCoefFct bulkViscositySlave = materials_[nitscheIf->GetSlaveVolRegion()]->GetScalCoefFnc(FLUID_BULK_VISCOSITY, Global::REAL);
  
    if( shearViscosityMaster->GetDependency() == CoefFunction::CONSTANT || shearViscositySlave->GetDependency() == CoefFunction::CONSTANT ) {
      double shearViscosityMasterVal, shearViscositySlaveVal, bulkViscosityMasterVal, bulkViscositySlaveVal;
      
      materials_[nitscheIf->GetMasterVolRegion()]->GetScalar(shearViscosityMasterVal,FLUID_DYNAMIC_VISCOSITY, Global::REAL);
      materials_[nitscheIf->GetSlaveVolRegion()]->GetScalar(shearViscositySlaveVal,FLUID_DYNAMIC_VISCOSITY, Global::REAL);
      
      materials_[nitscheIf->GetMasterVolRegion()]->GetScalar(bulkViscosityMasterVal,FLUID_BULK_VISCOSITY, Global::REAL);
      materials_[nitscheIf->GetSlaveVolRegion()]->GetScalar(bulkViscositySlaveVal,FLUID_BULK_VISCOSITY, Global::REAL);

      if( abs(shearViscosityMasterVal-shearViscositySlaveVal)>1e-14 || abs(bulkViscosityMasterVal-bulkViscositySlaveVal)>1e-14 ) {
        EXCEPTION("LinFlow ncInterfaces do not support different materials for master/slave region at the moment!")
      }
    } else {
      EXCEPTION("LinFlow ncInterfaces only support constant material parameters at the moment!")
    }

    PtrCoefFct shearViscosity = shearViscosityMaster;
    PtrCoefFct bulkViscosity = bulkViscosityMaster;
    PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
    PtrCoefFct constTwo = CoefFunction::Generate( mp_, Global::REAL, "2.0");

    PtrCoefFct shearViscosityDouble = CoefFunction::Generate( mp_,  Global::REAL,
        CoefXprBinOp(mp_, shearViscosity, CoefFunction::Generate( mp_, Global::REAL, "2"), CoefXpr::OP_MULT));
    PtrCoefFct coefZero = CoefFunction::Generate( mp_, Global::REAL, "0");
    StdVector<PtrCoefFct> tensorComponents(dim_ == 2 ? 9 : 36);
    tensorComponents.Init(coefZero);
    PtrCoefFct coefBB;
    if( dim_ == 2 ) {
      tensorComponents[0] = shearViscosityDouble;
      tensorComponents[4] = shearViscosityDouble;
      tensorComponents[8] = shearViscosity;
      coefBB = CoefFunction::Generate(mp_,Global::REAL,3,3,tensorComponents);
    } else {
      tensorComponents[0] = shearViscosityDouble;
      tensorComponents[7] = shearViscosityDouble;
      tensorComponents[14] = shearViscosityDouble;
      tensorComponents[21] = shearViscosity;
      tensorComponents[28] = shearViscosity;
      tensorComponents[35] = shearViscosity;
      coefBB = CoefFunction::Generate(mp_,Global::REAL,6,6,tensorComponents);
    }

    PtrCoefFct coefDivDiv = CoefFunction::Generate( mp_,  Global::REAL,
      CoefXprBinOp(mp_,
          bulkViscosity,
          CoefXprBinOp(mp_,shearViscosityDouble,CoefFunction::Generate( mp_, Global::REAL, "3"),CoefXpr::OP_DIV),
      CoefXpr::OP_SUB ));

    // get feFunctions
    shared_ptr<BaseFeFunction> velFct = feFunctions_[FLUIDMECH_VELOCITY];
    shared_ptr<BaseFeFunction> presFct = feFunctions_[FLUIDMECH_PRESSURE];

    //notation> assume the test function is called v, the velocity u and the pressure p
    //the notation differs slightly to the one used in the SinglePDE since the integrators use the test-function for the first entry and the fe-function for the second (enhanced readability)
    BiLinearForm *penalty_v1_u1 = NULL;
    BiLinearForm *penalty_v1_u2 = NULL;
    BiLinearForm *penalty_v2_u2 = NULL;
    //now bilinear forms related to the normal derivatives
    //du1 refers to the normal derivative directing from 1 to 2
    //for the flux part we have 3 individual terms which we have to take into account (main difference to standard stuff going on in the SinglePDE)
    // terms caused by -p I
    BiLinearForm *flux1_v1_p1 = NULL; // we have this term twice in the formulation (equivalent to flux1_p1_v1)
    BiLinearForm *flux1_v1_p2 = NULL;
    // terms caused by mu ( grad(v)+grad(v)^T )
    BiLinearForm *flux2_dv1_u1 = NULL;
    BiLinearForm *flux2_dv1_u2 = NULL;
    BiLinearForm *flux2_v1_du1 = NULL;
    // terms caused by (lambda - 2/3 mu) div(v) I
    BiLinearForm *flux3_dv1_u1 = NULL;
    BiLinearForm *flux3_dv1_u2 = NULL;
    BiLinearForm *flux3_v1_du1 = NULL;

    BiLinearForm::CouplingDirection curcpl;

    curcpl = BiLinearForm::MASTER_MASTER;

    // NOTE: the algebraic system sets the system matrix to
    // nonSym  if any bilinear form with the same fctID1 and fctID2 is nonSym.
    // We set here the symmetric flag to true in the constructor
    // of the SurfaceNitscheABInt even though the bilinear form itself is
    // not symmetric. Nitsche formulation is basically sym due to the
    // set counterpart directive for the context.
    
    // NOTE: Since the LinFlowPDE is non-sym per definition, the symmetrization might be unnecessary
    // TODO: Test this and if there is no difference, ommit the additional terms

    // penalty term
    if ( isMaterialComplex_) {
    	penalty_v1_u1 = new SurfaceNitscheABInt<Complex,Complex>
        	( new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
        	  new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
              shearViscosity, beta, curcpl, updatedGeo_, true, true);
    } else  {
      penalty_v1_u1 = new SurfaceNitscheABInt<Double,Double>
          ( new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
            new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
              shearViscosity, beta, curcpl, updatedGeo_, true, true);
    }

    // flux terms
    // 2 v1 p1 (due to symmetrization)
    if ( isMaterialComplex_) {
      flux1_v1_p1 = new SurfaceNitscheABInt<Complex,Complex>
          ( new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
            new SurfaceIdentityOperator<FeH1,DIM,P_DOF>(),
            constTwo, 1.0, curcpl, updatedGeo_, true);
    }
    else {
      flux1_v1_p1 = new SurfaceNitscheABInt<Double,Double>
          ( new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
            new SurfaceIdentityOperator<FeH1,DIM,P_DOF>(),
            constTwo, 1.0, curcpl, updatedGeo_, true);
    }

    // -mu ( grad(v1)+grad(v1)^T ) \cdot u1
    if ( isMaterialComplex_ ) {
      flux2_dv1_u1 = new SurfaceNitscheABInt<Complex,Complex>
          ( new SurfaceNormalStressOperator<FeH1,DIM,D_DOF>(subType_, false),
            new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
              constOne, -1.0, curcpl, updatedGeo_, true);
      flux2_dv1_u1->SetBCoefFunctionOpA(coefBB);
    } else {
      flux2_dv1_u1 = new SurfaceNitscheABInt<Double,Double>
          ( new SurfaceNormalStressOperator<FeH1,DIM,D_DOF>(subType_, false),
            new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
              constOne, -1.0, curcpl, updatedGeo_, true);
      flux2_dv1_u1->SetBCoefFunctionOpA(coefBB);
    }

    // -mu ( grad(u1)+grad(u1)^T ) \cdot v1
    if ( isMaterialComplex_ ) {
      flux2_v1_du1 = new SurfaceNitscheABInt<Complex,Complex>
          ( new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
            new SurfaceNormalStressOperator<FeH1,DIM,D_DOF>(subType_, false),
              constOne, -1.0, curcpl, updatedGeo_, true);
      flux2_v1_du1->SetBCoefFunctionOpB(coefBB);
    } else {
      flux2_v1_du1 = new SurfaceNitscheABInt<Double,Double>
          ( new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
            new SurfaceNormalStressOperator<FeH1,DIM,D_DOF>(subType_, false),
              constOne, -1.0, curcpl, updatedGeo_, true);
      flux2_v1_du1->SetBCoefFunctionOpB(coefBB);
    }

    // (lambda - 2/3 mu) div(v1) I \cdot u1
    if ( isMaterialComplex_ ) {
      flux3_dv1_u1 = new SurfaceNitscheABInt<Complex,Complex>
          ( new SurfaceNormalDivOperator<FeH1,DIM>(),
            new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
              coefDivDiv, -1.0, curcpl, updatedGeo_, true);
    } else {
      flux3_dv1_u1 = new SurfaceNitscheABInt<Double,Double>
          ( new SurfaceNormalDivOperator<FeH1,DIM>(),
            new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
              coefDivDiv, -1.0, curcpl, updatedGeo_, true);
    }

    // (lambda - 2/3 mu) div(u1) I \cdot v1
    if ( isMaterialComplex_ ) {
      flux3_v1_du1 = new SurfaceNitscheABInt<Complex,Complex>
          ( new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
            new SurfaceNormalDivOperator<FeH1,DIM>(),
              coefDivDiv, -1.0, curcpl, updatedGeo_, true);
    } else {
      flux3_v1_du1 = new SurfaceNitscheABInt<Double,Double>
          ( new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
            new SurfaceNormalDivOperator<FeH1,DIM>(),
              coefDivDiv, -1.0, curcpl, updatedGeo_, true);
    }


    curcpl = BiLinearForm::MASTER_SLAVE;

    // penalty term
    if ( isMaterialComplex_) {
    	penalty_v1_u2 = new SurfaceNitscheABInt<Complex,Complex>
        	( new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
        	  new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
              shearViscosity, beta * -1.0, curcpl, updatedGeo_, true, true);
    } else  {
      penalty_v1_u2 = new SurfaceNitscheABInt<Double,Double>
        ( new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
          new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
            shearViscosity, beta * -1.0, curcpl, updatedGeo_, true, true);
    }

    // flux terms
    // v1 p2
    if ( isMaterialComplex_) {
      flux1_v1_p2 = new SurfaceNitscheABInt<Complex,Complex>
          ( new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
            new SurfaceIdentityOperator<FeH1,DIM,P_DOF>(),
            constOne, -1.0, curcpl, updatedGeo_, true);
    }
    else {
      flux1_v1_p2 = new SurfaceNitscheABInt<Double,Double>
          ( new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
            new SurfaceIdentityOperator<FeH1,DIM,P_DOF>(),
            constOne, -1.0, curcpl, updatedGeo_, true);
    }

    // mu ( grad(v1)+grad(v1)^T ) \cdot u2
    if ( isMaterialComplex_ ) {
      flux2_dv1_u2 = new SurfaceNitscheABInt<Complex,Complex>
          ( new SurfaceNormalStressOperator<FeH1,DIM,D_DOF>(subType_, false),
            new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
              constOne, 1.0, curcpl, updatedGeo_, true);
      flux2_dv1_u2->SetBCoefFunctionOpA(coefBB);
    } else {
      flux2_dv1_u2 = new SurfaceNitscheABInt<Double,Double>
          ( new SurfaceNormalStressOperator<FeH1,DIM,D_DOF>(subType_, false),
            new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
              constOne, 1.0, curcpl, updatedGeo_, true);
      flux2_dv1_u2->SetBCoefFunctionOpA(coefBB);
    }
    
    // (lambda - 2/3 mu) div(v1) I \cdot u2
    if ( isMaterialComplex_ ) {
      flux3_dv1_u2 = new SurfaceNitscheABInt<Complex,Complex>
          ( new SurfaceNormalDivOperator<FeH1,DIM>(),
            new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
              coefDivDiv, 1.0, curcpl, updatedGeo_, true);
    } else {
      flux3_dv1_u2 = new SurfaceNitscheABInt<Double,Double>
          ( new SurfaceNormalDivOperator<FeH1,DIM>(),
            new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
              coefDivDiv, 1.0, curcpl, updatedGeo_, true);
    }


    curcpl = BiLinearForm::SLAVE_SLAVE;

    // penalty term
    if ( isMaterialComplex_) {
    	penalty_v2_u2 = new SurfaceNitscheABInt<Complex,Complex>
        	( new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
        	  new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
              shearViscosity, beta, curcpl, updatedGeo_, true, true);
    } else  {
      penalty_v2_u2 = new SurfaceNitscheABInt<Double,Double>
        ( new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
          new SurfaceIdentityOperator<FeH1,DIM,D_DOF>(),
            shearViscosity, beta, curcpl, updatedGeo_, true, true);
    }


    // penalty terms
    SurfaceBiLinFormContext *penalty_v1_u1_Context = NULL;
    SurfaceBiLinFormContext *penalty_v1_u2_Context = NULL;
    SurfaceBiLinFormContext *penalty_v2_u2_Context = NULL;
    // terms caused by -p I
    SurfaceBiLinFormContext *flux1_v1_p1_Context = NULL;
    SurfaceBiLinFormContext *flux1_v1_p2_Context = NULL;
    // terms caused by mu ( grad(v)+grad(v)^T )
    SurfaceBiLinFormContext *flux2_dv1_u1_Context = NULL;
    SurfaceBiLinFormContext *flux2_dv1_u2_Context = NULL;
    SurfaceBiLinFormContext *flux2_v1_du1_Context = NULL;
    // terms caused by (lambda - 2/3 mu) div(v) I
    SurfaceBiLinFormContext *flux3_dv1_u1_Context = NULL;
    SurfaceBiLinFormContext *flux3_dv1_u2_Context = NULL;
    SurfaceBiLinFormContext *flux3_v1_du1_Context = NULL;


    FEMatrixType targetMatrix = STIFFNESS;
    if(isMoving){
    	targetMatrix = STIFFNESS_UPDATE;
    }

    curcpl = BiLinearForm::MASTER_MASTER;
    penalty_v1_u1_Context = new SurfaceBiLinFormContext(penalty_v1_u1, targetMatrix, curcpl);
    flux1_v1_p1_Context = new SurfaceBiLinFormContext(flux1_v1_p1, targetMatrix, curcpl);
    flux2_dv1_u1_Context = new SurfaceBiLinFormContext(flux2_dv1_u1, targetMatrix, curcpl);
    flux2_v1_du1_Context = new SurfaceBiLinFormContext(flux2_v1_du1, targetMatrix, curcpl);
    flux3_dv1_u1_Context = new SurfaceBiLinFormContext(flux3_dv1_u1, targetMatrix, curcpl);
    flux3_v1_du1_Context = new SurfaceBiLinFormContext(flux3_v1_du1, targetMatrix, curcpl);

    curcpl = BiLinearForm::MASTER_SLAVE;
    penalty_v1_u2_Context = new SurfaceBiLinFormContext(penalty_v1_u2, targetMatrix, curcpl);
    flux1_v1_p2_Context = new SurfaceBiLinFormContext(flux1_v1_p2, targetMatrix, curcpl);
    flux2_dv1_u2_Context = new SurfaceBiLinFormContext(flux2_dv1_u2, targetMatrix, curcpl);
    flux3_dv1_u2_Context = new SurfaceBiLinFormContext(flux3_dv1_u2, targetMatrix, curcpl);

    curcpl = BiLinearForm::SLAVE_SLAVE;
    penalty_v2_u2_Context = new SurfaceBiLinFormContext(penalty_v2_u2, targetMatrix, curcpl);


    // for moving meshes we need to update the motion of the interface
    if (isMoving) {
      penalty_v1_u1_Context->SetMotion(true);
      penalty_v1_u2_Context->SetMotion(true);
      penalty_v2_u2_Context->SetMotion(true);
      flux1_v1_p1_Context->SetMotion(true);
      flux1_v1_p2_Context->SetMotion(true);

      flux2_dv1_u1_Context->SetMotion(true);
      flux2_dv1_u2_Context->SetMotion(true);
      flux2_v1_du1_Context->SetMotion(true);

      flux3_dv1_u1_Context->SetMotion(true);
      flux3_dv1_u2_Context->SetMotion(true);
      flux3_v1_du1_Context->SetMotion(true);
    }

    // set names
    penalty_v1_u1->SetName("penalty_v1_u1");
    penalty_v1_u2->SetName("penalty_v1_u2");
    penalty_v2_u2->SetName("penalty_v2_u2");

    flux1_v1_p1->SetName("flux1_v1_p1");
    flux1_v1_p2->SetName("flux1_v1_p2");

    flux2_dv1_u1->SetName("flux2_dv1_u1");
    flux2_dv1_u2->SetName("flux2_dv1_u2");
    flux2_v1_du1->SetName("flux2_v1_du1");

    flux3_dv1_u1->SetName("flux3_dv1_u1");
    flux3_dv1_u2->SetName("flux3_dv1_u2");
    flux3_v1_du1->SetName("flux3_v1_du1");

    // set entities
    penalty_v1_u1_Context->SetEntities(actSDList,actSDList);
    penalty_v1_u2_Context->SetEntities(actSDList,actSDList);
    penalty_v2_u2_Context->SetEntities(actSDList,actSDList);
    
    flux1_v1_p1_Context->SetEntities(actSDList,actSDList);
    flux1_v1_p2_Context->SetEntities(actSDList,actSDList);

    flux2_dv1_u1_Context->SetEntities(actSDList,actSDList);
    flux2_dv1_u2_Context->SetEntities(actSDList,actSDList);
    flux2_v1_du1_Context->SetEntities(actSDList,actSDList);

    flux3_dv1_u1_Context->SetEntities(actSDList,actSDList);
    flux3_dv1_u2_Context->SetEntities(actSDList,actSDList);
    flux3_v1_du1_Context->SetEntities(actSDList,actSDList);

    // set feFunctions
    penalty_v1_u1_Context->SetFeFunctions( velFct, velFct );
    penalty_v1_u2_Context->SetFeFunctions( velFct, velFct );
    penalty_v2_u2_Context->SetFeFunctions( velFct, velFct );

    flux1_v1_p1_Context->SetFeFunctions( velFct, presFct );
    flux1_v1_p2_Context->SetFeFunctions( velFct, presFct );

    flux2_dv1_u1_Context->SetFeFunctions( velFct, velFct );
    flux2_dv1_u2_Context->SetFeFunctions( velFct, velFct );
    flux2_v1_du1_Context->SetFeFunctions( velFct, velFct );
    
    flux3_dv1_u1_Context->SetFeFunctions( velFct, velFct );
    flux3_dv1_u2_Context->SetFeFunctions( velFct, velFct );
    flux3_v1_du1_Context->SetFeFunctions( velFct, velFct );

    // add counter part (for symmetry)
    // NOTE: Since the overall formulation is not sym, we could also ommit this, but before doing so, check if the previously defined integrators are the one for traction consistency!!
    penalty_v1_u2_Context->SetCounterPart(true);
    flux1_v1_p2_Context->SetCounterPart(true);
    flux2_dv1_u2_Context->SetCounterPart(true);
    flux3_dv1_u2_Context->SetCounterPart(true);

    // add to assemble
    assemble_->AddBiLinearForm( penalty_v1_u1_Context );
    assemble_->AddBiLinearForm( penalty_v1_u2_Context );
    assemble_->AddBiLinearForm( penalty_v2_u2_Context );

    assemble_->AddBiLinearForm( flux1_v1_p1_Context );
    assemble_->AddBiLinearForm( flux1_v1_p2_Context );

    assemble_->AddBiLinearForm( flux2_dv1_u1_Context );
    assemble_->AddBiLinearForm( flux2_dv1_u2_Context );
    assemble_->AddBiLinearForm( flux2_v1_du1_Context );

    assemble_->AddBiLinearForm( flux3_dv1_u1_Context );
    assemble_->AddBiLinearForm( flux3_dv1_u2_Context );
    assemble_->AddBiLinearForm( flux3_v1_du1_Context );


    ncIf->RegisterIntegrator( penalty_v1_u1_Context );
    ncIf->RegisterIntegrator( penalty_v1_u2_Context );
    ncIf->RegisterIntegrator( penalty_v2_u2_Context );

    ncIf->RegisterIntegrator( flux1_v1_p1_Context );
    ncIf->RegisterIntegrator( flux1_v1_p2_Context );

    ncIf->RegisterIntegrator( flux2_dv1_u1_Context );
    ncIf->RegisterIntegrator( flux2_dv1_u2_Context );
    ncIf->RegisterIntegrator( flux2_v1_du1_Context );

    ncIf->RegisterIntegrator( flux3_dv1_u1_Context );
    ncIf->RegisterIntegrator( flux3_dv1_u2_Context );
    ncIf->RegisterIntegrator( flux3_v1_du1_Context );

  }

  void LinFlowPDE::DefineSurfaceIntegrators(){
    // Get FESpace and FeFunction of fluid velocity
    shared_ptr<BaseFeFunction> velFct = feFunctions_[FLUIDMECH_VELOCITY];
    shared_ptr<FeSpace> velSpace = velFct->GetFeSpace();

    // Get FESpace and FeFunction of fluid pressure
    shared_ptr<BaseFeFunction> presFct = feFunctions_[FLUIDMECH_PRESSURE];
    shared_ptr<FeSpace> presSpace = presFct->GetFeSpace();

    // Get FESpace and FeFunction of the lagrange multiplier
    shared_ptr<BaseFeFunction> lagrangeMultFct = NULL;
    shared_ptr<FeSpace> lagrangeMultSpace = NULL;

    shared_ptr<BaseFeFunction> lagrangeMultFct_1 = NULL;
    shared_ptr<FeSpace> lagrangeMultSpace_1 = NULL;

    PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    if( bcNode ) {
      StdVector<shared_ptr<EntityList> > ent;
      StdVector<PtrCoefFct > kCoef;
      StdVector<std::string> volumeRegions;
      //========================================================================================
      // normalSurfaceMass : assumes a normal traction proprtional to the normal velocity
      //========================================================================================
      ReadRhsExcitation("normalSurfaceMass", feFunctions_[FLUIDMECH_VELOCITY]->GetResultInfo()->dofNames,
          ResultInfo::SCALAR, isComplex_, ent, kCoef, updatedGeo_, volumeRegions);
      for( UInt i = 0; i < ent.GetSize(); ++i ) {
        // get the volume region for defining the correct normal direction
        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volumeRegions[i]);
        std::set<RegionIdType> volRegion;
        volRegion.insert(aRegion);
        // check type of entitylist
        if (ent[i]->GetType() == EntityList::NODE_LIST) {
          EXCEPTION("normalSurfaceMass must be defined on (surface) elements")
        }
        // setup the integrator for: u'*t_n = u'*k u_n = u'*(k u*n n) = k u'*n u*n
        BiLinearForm * surfMassInt = NULL;
        if(isComplex_) {
          if (dim_ == 2){
            surfMassInt = new SurfaceBBInt<Complex,Complex>(new IdentityOperatorNormalTrans<FeH1,2,2,Complex>(), kCoef[i], Complex(balanceOfMomentumSign_,0), volRegion, updatedGeo_ );
          } else {
            surfMassInt = new SurfaceBBInt<Complex,Complex>(new IdentityOperatorNormalTrans<FeH1,3,3,Complex>(), kCoef[i], Complex(balanceOfMomentumSign_,0), volRegion, updatedGeo_ );
          }
        } else {
          if (dim_ == 2){
            surfMassInt = new SurfaceBBInt<>(new IdentityOperatorNormalTrans<FeH1,2,2>(), kCoef[i], balanceOfMomentumSign_, volRegion, updatedGeo_ );
          } else {
            surfMassInt = new SurfaceBBInt<>(new IdentityOperatorNormalTrans<FeH1,3,3>(), kCoef[i], balanceOfMomentumSign_, volRegion, updatedGeo_ );
          }
        }
        surfMassInt->SetName("surfMassIntegrator");
        BiLinFormContext *surfMassContext = new BiLinFormContext(surfMassInt, DAMPING );
        surfMassContext->SetEntities( ent[i], ent[i]);
        surfMassContext->SetFeFunctions( feFunctions_[FLUIDMECH_VELOCITY], feFunctions_[FLUIDMECH_VELOCITY]);
        feFunctions_[FLUIDMECH_VELOCITY]->AddEntityList( ent[i] );
        assemble_->AddBiLinearForm( surfMassContext );
      }
      //========================================================================================
      //Impedance BC : sigma.n = z_0 (v.n).n
      // high z_0 : wall(sound hard) BC
      // low z_0 : sound soft BC
      //========================================================================================
      ReadRhsExcitation("normalImpedance", feFunctions_[FLUIDMECH_VELOCITY]->GetResultInfo()->dofNames,
          ResultInfo::SCALAR, isComplex_, ent, kCoef, updatedGeo_, volumeRegions);
      for( UInt i = 0; i < ent.GetSize(); ++i ) {
        // get the volume region for defining the correct normal direction
        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volumeRegions[i]);
        std::set<RegionIdType> volRegion;
        volRegion.insert(aRegion);
        // check type of entitylist
        if (ent[i]->GetType() == EntityList::NODE_LIST) {
          EXCEPTION("normalImpedance must be defined on (surface) elements")
        }
        // setup the integrator for: u'*t_n = u'*k u_n = u'*(k u*n n) = k u'*n u*n
        BiLinearForm * impedanceInt = NULL;
        if(isComplex_) {
          if (dim_ == 2){
            impedanceInt = new SurfaceBBInt<Complex,Complex>(new IdentityOperatorNormalTrans<FeH1,2,2,Complex>(), kCoef[i], Complex(balanceOfMomentumSign_,0), volRegion, updatedGeo_ );
          } else {
            impedanceInt = new SurfaceBBInt<Complex,Complex>(new IdentityOperatorNormalTrans<FeH1,3,3,Complex>(), kCoef[i], Complex(balanceOfMomentumSign_,0), volRegion, updatedGeo_ );
          }
        } else {
          if (dim_ == 2){
            impedanceInt = new SurfaceBBInt<>(new IdentityOperatorNormalTrans<FeH1,2,2>(), kCoef[i], balanceOfMomentumSign_, volRegion, updatedGeo_ );
          } else {
            impedanceInt = new SurfaceBBInt<>(new IdentityOperatorNormalTrans<FeH1,3,3>(), kCoef[i], balanceOfMomentumSign_, volRegion, updatedGeo_ );
          }
        }
        impedanceInt->SetName("ImpedanceIntegrator");
        BiLinFormContext *impedanceContext = new BiLinFormContext(impedanceInt, STIFFNESS );
        impedanceContext->SetEntities( ent[i], ent[i]);
        impedanceContext->SetFeFunctions( feFunctions_[FLUIDMECH_VELOCITY], feFunctions_[FLUIDMECH_VELOCITY]);
        feFunctions_[FLUIDMECH_VELOCITY]->AddEntityList( ent[i] );
        assemble_->AddBiLinearForm( impedanceContext );
      }
      //========================================================================================
      // Absorbing BC : sigma.n = z_0 (v.n).n
      // z_0 is the LinFLow impedance at the boundary
      // for traveling waves:
      // z_0 = rho*c /(1+real(complexWaveNum)/complexWaveNum)
      //complexWaveNumber = sqrt((omega^2*-rho/(K+(4/3 shearVisc+ bulkVisc)*jomega)))
      //========================================================================================
      ParamNodeList abcNodes = bcNode->GetList( "absorbingBCs" );
      if ( this->analysistype_ != HARMONIC && !abcNodes.IsEmpty() )
        EXCEPTION("AbsorbingBCs are only allowed in harmonic analysis");
      // omega = 2*pi*f
      PtrCoefFct omega = CoefFunction::Generate( mp_, Global::REAL,CoefXprBinOp(mp_,
          CoefFunction::Generate( mp_,  Global::REAL, "pi*f"),
          CoefFunction::Generate( mp_, Global::REAL, "2"),CoefXpr::OP_MULT));

      PtrCoefFct jomega = CoefFunction::Generate( mp_, Global::COMPLEX,CoefXprBinOp(mp_,
          omega, CoefFunction::Generate( mp_, Global::COMPLEX, "0", "1"),
          CoefXpr::OP_MULT));

      PtrCoefFct omega2 = CoefFunction::Generate( mp_, Global::COMPLEX,CoefXprBinOp(mp_,
          omega, omega,CoefXpr::OP_MULT));
      //========================================================================================
      // complexWaveNumber = sqrt(complexWaveNumber1)
      // https://math.stackexchange.com/questions/44406/how-do-i-get-the-square-root-of-a-complex-number
      // First: complex wave number is calculated
      // Second: impedance is computed as rhoC/(1+real(complexWaveNum)/complexWaveNum)
      //========================================================================================
      for( UInt i = 0; i < abcNodes.GetSize(); i++ ) {
        std::string regionName = abcNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );
        std::string volRegName = abcNodes[i]->Get("volumeRegion")->As<std::string>();
        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);
        std::set<RegionIdType> volRegion;
        volRegion.insert(aRegion);
        PtrCoefFct density = materials_[aRegion]->GetScalCoefFnc(DENSITY, Global::REAL);
        PtrCoefFct compressionModulus = materials_[aRegion]->GetScalCoefFnc( FLUID_BULK_MODULUS, Global::REAL );
        PtrCoefFct shearViscosity = materials_[aRegion]->GetScalCoefFnc(FLUID_DYNAMIC_VISCOSITY, Global::REAL); // mu
        PtrCoefFct bulkViscosity = materials_[aRegion]->GetScalCoefFnc(FLUID_BULK_VISCOSITY, Global::REAL);
        // visc = (4/3 shearVisc+ bulkVisc)
        PtrCoefFct visc = CoefFunction::Generate( mp_,  Global::REAL,
            CoefXprBinOp(mp_,
                bulkViscosity,
                CoefXprBinOp(mp_,
                    CoefFunction::Generate( mp_,
                        Global::REAL,CoefXprBinOp(mp_,
                            shearViscosity, CoefFunction::Generate( mp_, Global::REAL, "4"),
                            CoefXpr::OP_MULT)),
                            CoefFunction::Generate( mp_, Global::REAL, "3"),
                            CoefXpr::OP_DIV),
                            CoefXpr::OP_ADD ));

        // complexWaveNum2 = K +(4/3 shearVisc+ bulkVisc)*jomega
        PtrCoefFct complexWaveNum2 = CoefFunction::Generate(mp_,  Global::COMPLEX,
            CoefXprBinOp(mp_,
                CoefFunction::Generate(mp_,
                    Global::COMPLEX,CoefXprBinOp(mp_,
                        visc, jomega,
                        CoefXpr::OP_MULT)),
                        compressionModulus,
                        CoefXpr::OP_ADD));

        //complexWaveNumber1 = -rho*omega^2/(K+(4/3 shearVisc+ bulkVisc)*jomega))
        PtrCoefFct complexWaveNum1 = CoefFunction::Generate(mp_, Global::COMPLEX,
            CoefXprBinOp(mp_,
                omega2,
                CoefFunction::Generate(mp_,
                    Global::COMPLEX,CoefXprBinOp(mp_,
                        CoefFunction::Generate( mp_,
                            Global::COMPLEX,CoefXprBinOp(mp_,
                                density, CoefFunction::Generate( mp_, Global::REAL, "-1"),
                                CoefXpr::OP_MULT)), complexWaveNum2,
                                CoefXpr::OP_DIV)),
                                CoefXpr::OP_MULT ));
        //complexWaveNumber = sqrt[-rho*omega^2/(K+(4/3 shearVisc+ bulkVisc)*jomega))]                       
        PtrCoefFct complexWaveNum = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprUnaryOp(mp_, complexWaveNum1, CoefXpr::OP_SQRT));
        PtrCoefFct complexWaveNumRe = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_, complexWaveNum, CoefXpr::OP_RE));
        // rhoC = rho*C =sqrt(rho*k)
        PtrCoefFct rhoC = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_,
            CoefFunction::Generate( mp_, Global::REAL,CoefXprBinOp(mp_, density, compressionModulus, CoefXpr::OP_MULT)),
            CoefXpr::OP_SQRT));
        //complex impedance of a viscous fluid:
        //impCoeff = rhoC/(1+real(complexWaveNum)/complexWaveNum)
        PtrCoefFct impCoeff = CoefFunction::Generate( mp_, Global::COMPLEX,CoefXprBinOp(mp_, rhoC,
            CoefFunction::Generate( mp_, Global::COMPLEX,CoefXprBinOp(mp_, CoefFunction::Generate( mp_, Global::REAL, "1"),
                CoefFunction::Generate( mp_, Global::COMPLEX,CoefXprBinOp(mp_, complexWaveNumRe, complexWaveNum, CoefXpr::OP_DIV)),
                CoefXpr::OP_ADD)),
                CoefXpr::OP_DIV));

        // setup the integrator for: u'*t_n = u'*k u_n = u'*(k u*n n) = z_0 u'*n u*n
        BiLinearForm * abcInt = NULL;
        if(isComplex_) {
          if (dim_ == 2){
            abcInt = new SurfaceBBInt<Complex,Complex>(new IdentityOperatorNormalTrans<FeH1,2,2,Complex>(), impCoeff, Complex(balanceOfMomentumSign_,0), volRegion, updatedGeo_ );
          } else {
            abcInt = new SurfaceBBInt<Complex,Complex>(new IdentityOperatorNormalTrans<FeH1,3,3,Complex>(), impCoeff, Complex(balanceOfMomentumSign_,0), volRegion, updatedGeo_ );
          }
        } else {
          if (dim_ == 2){
            abcInt = new SurfaceBBInt<>(new IdentityOperatorNormalTrans<FeH1,2,2>(), impCoeff, balanceOfMomentumSign_, volRegion, updatedGeo_ );
          } else {
            abcInt = new SurfaceBBInt<>(new IdentityOperatorNormalTrans<FeH1,3,3>(), impCoeff, balanceOfMomentumSign_, volRegion, updatedGeo_ );
          }
        }
        abcInt->SetName("abcIntegrator");
        BiLinFormContext *abcContext = new BiLinFormContext(abcInt, STIFFNESS );
        abcContext->SetEntities( actSDList, actSDList);
        abcContext->SetFeFunctions( feFunctions_[FLUIDMECH_VELOCITY], feFunctions_[FLUIDMECH_VELOCITY]);
        feFunctions_[FLUIDMECH_VELOCITY]->AddEntityList( actSDList );
        assemble_->AddBiLinearForm( abcContext );
      }
      
      
      //========================================================================================
      // Do nothing BC
      // Implements the classical do nothing BC which is also required for LM based BCs
      // Although we could provide the doNothing BC in the xml scheme, we restrict it to be an
      // "inner" function. The reason for this is that the formulation seems unstable when
      // inflow occours, hence, the usage besides the weakly enforced Maxwell BC is questionable
      //========================================================================================

      ReadRhsExcitation("doNothing", feFunctions_[FLUIDMECH_VELOCITY]->GetResultInfo()->dofNames,
          ResultInfo::SCALAR, isComplex_, ent, kCoef, updatedGeo_, volumeRegions);

      for( UInt i = 0; i < ent.GetSize(); ++i ) {
        // Axi-symmetric formulation is not finished yet, throw an exception of requested
        // Note: The strain operator is already updated to handle axi-symmetric stuff, only the surface div-operator is missing
        if( isaxi_ ) {
          EXCEPTION("Axi-symmtric formulation not implemented for the doNothing BC (missing surface div-operator in axi-formulation, please implement me)");
        }
        
        AddSurfaceIntegratorFromPartialIntegration(i, volumeRegions, ent, "all");
        
      }


      //========================================================================================
      // velocity constraint
      //========================================================================================

      ReadRhsExcitation("velocityConstraint", feFunctions_[FLUIDMECH_VELOCITY]->GetResultInfo()->dofNames,
          ResultInfo::SCALAR, isComplex_, ent, kCoef, updatedGeo_, volumeRegions);

      // Define FE space etc.
      if ( useLagrangeMultVec_ ) {
        lagrangeMultFct = feFunctions_[LAGRANGE_MULT];
        lagrangeMultSpace = lagrangeMultFct->GetFeSpace();
      }

      for( UInt i = 0; i < ent.GetSize(); ++i ) {
        // get the volume region for defining the correct normal direction
        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volumeRegions[i]);
        std::set<RegionIdType> volRegion;
        volRegion.insert(aRegion);
        // check type of entitylist
        if (ent[i]->GetType() != EntityList::SURF_ELEM_LIST) {
          EXCEPTION("Velocity constraint BC must be defined on (surface) elements")
        }
        std::string regionName = ent[i]->GetName();
        shared_ptr<EntityList> actSDList = ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );

        // Get polynomial and integration order for Lagrange multiplier
        lagrangeMultSpace->SetRegionApproximation(aRegion, lagrangeMultPolyId_, lagrangeMultIntegId_);


        PtrParamNode velConstNode = bcNode->GetByVal("velocityConstraint","name", regionName);
        ParamNodeList velConstNodeList = velConstNode->GetChildren();

        bool isNoSlipBC = false;
        bool isMaxwellBC = false;
        for (UInt ii = 0; ii < velConstNodeList.GetSize(); ii++){
          std::string nodeName = velConstNodeList[ii]->GetName();
          if( nodeName=="noSlip" ){
            isNoSlipBC = true;
          } else if( nodeName=="MaxwellFirstOrderSlip" ) {
            isMaxwellBC = true;
          }
        }

        // setup the integrators


        if( isNoSlipBC || isMaxwellBC ) {
          // since the velocity is enforced weakly we have to add the missing surface term stemming from partial integration
          // similar to the noPenetration BC we only add the pressure term, otherwise stuff gets unstable for the noSlip BC
          AddSurfaceIntegratorFromPartialIntegration(i, volumeRegions, ent, "P1");


          // conservation of momentum: \int_\Gamma_D v' \lambda
          BiLinearForm * surfVelocityConstraintIntVL = NULL;
          PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
          if(isComplex_) {
            if (dim_ == 2){
              surfVelocityConstraintIntVL = new SurfaceABInt<Complex,Complex>(new IdentityOperator<FeH1,2,2>(),
                                                                              new IdentityOperator<FeH1,2,2>(),
                                                                              constOne, Complex(balanceOfMomentumSign_,0), volRegion, updatedGeo_ );
            } else {
              surfVelocityConstraintIntVL = new SurfaceABInt<Complex,Complex>(new IdentityOperator<FeH1,3,3>(),
                                                                              new IdentityOperator<FeH1,3,3>(),
                                                                              constOne, Complex(balanceOfMomentumSign_,0), volRegion, updatedGeo_ );
            }
          } else {
            if (dim_ == 2){
              surfVelocityConstraintIntVL = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                                                new IdentityOperator<FeH1,2,2>(),
                                                                constOne, balanceOfMomentumSign_, volRegion, updatedGeo_ );
            } else {
              surfVelocityConstraintIntVL = new SurfaceABInt<>(new IdentityOperator<FeH1,3,3>(),
                                                                new IdentityOperator<FeH1,3,3>(),
                                                                constOne, balanceOfMomentumSign_, volRegion, updatedGeo_ );
            }
          }
          surfVelocityConstraintIntVL->SetName("velocityConstraintIntVL");
          BiLinFormContext *surfVelocityConstraintContextVL = new BiLinFormContext(surfVelocityConstraintIntVL, STIFFNESS );
          surfVelocityConstraintContextVL->SetEntities( actSDList, actSDList);
          surfVelocityConstraintContextVL->SetFeFunctions( feFunctions_[FLUIDMECH_VELOCITY], feFunctions_[LAGRANGE_MULT]);
          feFunctions_[FLUIDMECH_VELOCITY]->AddEntityList( actSDList );
          feFunctions_[LAGRANGE_MULT]->AddEntityList( actSDList );
          assemble_->AddBiLinearForm( surfVelocityConstraintContextVL );


          // constraint equation: \int_\Gamma_D \lambda' \cdot  v 
          BiLinearForm * surfVelocityConstraintIntLV = NULL;
          if(isComplex_) {
            if (dim_ == 2){
              surfVelocityConstraintIntLV = new SurfaceABInt<Complex,Complex>(new IdentityOperator<FeH1,2,2>(),
                                                                              new IdentityOperator<FeH1,2,2>(),
                                                                              constOne, Complex(balanceOfMomentumSign_,0), volRegion, updatedGeo_ );
            } else {
              surfVelocityConstraintIntLV = new SurfaceABInt<Complex,Complex>(new IdentityOperator<FeH1,3,3>(),
                                                                              new IdentityOperator<FeH1,3,3>(),
                                                                              constOne, Complex(balanceOfMomentumSign_,0), volRegion, updatedGeo_ );
            }
          } else {
            if (dim_ == 2){
              surfVelocityConstraintIntLV = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                                                new IdentityOperator<FeH1,2,2>(),
                                                                constOne, balanceOfMomentumSign_, volRegion, updatedGeo_ );
            } else {
              surfVelocityConstraintIntLV = new SurfaceABInt<>(new IdentityOperator<FeH1,3,3>(),
                                                                new IdentityOperator<FeH1,3,3>(),
                                                                constOne, balanceOfMomentumSign_, volRegion, updatedGeo_ );
            }
          }
          surfVelocityConstraintIntLV->SetName("velocityConstraintIntLV");
          BiLinFormContext *surfVelocityConstraintContextLV = new BiLinFormContext(surfVelocityConstraintIntLV, STIFFNESS );
          surfVelocityConstraintContextLV->SetEntities( actSDList, actSDList);
          surfVelocityConstraintContextLV->SetFeFunctions( feFunctions_[LAGRANGE_MULT], feFunctions_[FLUIDMECH_VELOCITY]);
          feFunctions_[FLUIDMECH_VELOCITY]->AddEntityList( actSDList );
          feFunctions_[LAGRANGE_MULT]->AddEntityList( actSDList );
          assemble_->AddBiLinearForm( surfVelocityConstraintContextLV );
        }


        if( isMaxwellBC ){
          std::string maxwellMeanFreePathString;
          std::string maxwellC1String;
          std::string maxwellFormulationString;
          velConstNode->Get("MaxwellFirstOrderSlip")->GetValue( "meanFreePath", maxwellMeanFreePathString );
          velConstNode->Get("MaxwellFirstOrderSlip")->GetValue( "C1", maxwellC1String );
          velConstNode->Get("MaxwellFirstOrderSlip")->GetValue( "formulation", maxwellFormulationString );

          //double meanFreePath = std::stod(maxwellMeanFreePathString);
          //double maxwellC1 = std::stod(maxwellC1String);
          
          PtrCoefFct coefMaxwellBC = CoefFunction::Generate( mp_,  Global::REAL,
            CoefXprBinOp(mp_, CoefFunction::Generate( mp_, Global::REAL, maxwellC1String),
            CoefFunction::Generate( mp_, Global::REAL, maxwellMeanFreePathString), CoefXpr::OP_MULT));

          
          // for the Maxwell BC we also add the rest of the terms stemming from partial integration
          // this is necessary in order to not get a mesh dependent BC, which would be the case if the missing terms are neglected
          // for the Maxwell BC adding these terms does not seem to introduce instability when the normal vector is parallel to an axis
          AddSurfaceIntegratorFromPartialIntegration(i, volumeRegions, ent, "P2");
          AddSurfaceIntegratorFromPartialIntegration(i, volumeRegions, ent, "P3");
          
          
          // constraint equation: \int_\Gamma_D \lambda^\star \cdot  M(u)
          // here M() denotes the operator used for the calculation of the slip velocity based on the normal derivative of the velocity 
          // to be consistent with the derivation we also set the counterpart for this integrator
          // since we used the weak noSlip BC as a basis and the Maxwell part only acts in the tangential plane, we also enforce a noPenetration BC
          BiLinearForm * surfVelocityConstraintMaxwellIntLV = NULL;
          if(isComplex_) {
            if (dim_ == 2){
              if (maxwellFormulationString == "hollowIncompressibleStress") {
                surfVelocityConstraintMaxwellIntLV = new SurfaceABInt<Complex,Complex>(new IdentityOperator<FeH1,2,2>(),
                                                                                        new SurfaceTangentialHollowIncompressibleStrainOperator2D<FeH1,2,2>(),
                                                                                        coefMaxwellBC, Complex(balanceOfMomentumSign_,0), volRegion, updatedGeo_ );
              } else if (maxwellFormulationString == "incompressibleStress") {
                surfVelocityConstraintMaxwellIntLV = new SurfaceABInt<Complex,Complex>(new IdentityOperator<FeH1,2,2>(),
                                                                                        new SurfaceTangentialIncompressibleStrainOperator2D<FeH1,2,2>(),
                                                                                        coefMaxwellBC, Complex(balanceOfMomentumSign_,0), volRegion, updatedGeo_ );
              } else if (maxwellFormulationString == "compressibleStress") {
                surfVelocityConstraintMaxwellIntLV = new SurfaceABInt<Complex,Complex>(new IdentityOperator<FeH1,2,2>(),
                                                                                        new SurfaceTangentialCompressibleStrainOperator2D<FeH1,2,2>(),
                                                                                        coefMaxwellBC, Complex(balanceOfMomentumSign_,0), volRegion, updatedGeo_ );
              } else {
                EXCEPTION("Unknown formulation for the Maxwell slip BC!");
              }
            } else {
              EXCEPTION("3D weakly enforced Dirichlet BCs not available for LinFlow, please implement me!");
            }
          } else {
            if (dim_ == 2) {
              if (maxwellFormulationString == "hollowIncompressibleStress") {
                surfVelocityConstraintMaxwellIntLV = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                                                        new SurfaceTangentialHollowIncompressibleStrainOperator2D<FeH1,2,2>(),
                                                                        coefMaxwellBC, balanceOfMomentumSign_, volRegion, updatedGeo_ );
              } else if (maxwellFormulationString == "incompressibleStress") {
                surfVelocityConstraintMaxwellIntLV = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                                                        new SurfaceTangentialIncompressibleStrainOperator2D<FeH1,2,2>(),
                                                                        coefMaxwellBC, balanceOfMomentumSign_, volRegion, updatedGeo_ );
              } else if (maxwellFormulationString == "compressibleStress") {
                surfVelocityConstraintMaxwellIntLV = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                                                        new SurfaceTangentialCompressibleStrainOperator2D<FeH1,2,2>(),
                                                                        coefMaxwellBC, balanceOfMomentumSign_, volRegion, updatedGeo_ );
              }
              else {
                EXCEPTION("Unknown formulation for the Maxwell slip BC!");
              }
            } else {
              EXCEPTION("3D weakly enforced Dirichlet BCs not available for LinFlow, please implement me!");
            }
          }

          //surfWeakVelocityMaxwellIntLV->SetRequestVolElem();
          surfVelocityConstraintMaxwellIntLV->SetName("velocityConstraintMaxwellIntLV");
          // we set the volume evaluation for the BiLinForm since the context will automatically use the value of the BiLinForm
          // we have to do this before the initialization of the context!
          surfVelocityConstraintMaxwellIntLV->SetUseVolEqnB( true ); 

          //Adding the coefficient lambda/mu to the compressible Term 
          PtrCoefFct shearViscosity = materials_[aRegion]->GetScalCoefFnc(FLUID_DYNAMIC_VISCOSITY, Global::REAL); // mu
          PtrCoefFct bulkViscosity =  materials_[aRegion]->GetScalCoefFnc(FLUID_BULK_VISCOSITY, Global::REAL); //lambda
          //PtrCoefFct bulkViscosity = CoefFunction::Generate( mp_, Global::REAL, "0.0"); // Setting lambda to 0 for the incompressible flow forumation
          // if ( isCompressible_ ) {
          //  bulkViscosity = materials_[aRegion]->GetScalCoefFnc(FLUID_BULK_VISCOSITY, Global::REAL);
          //}
          
          //Generating CoefFunktion: lambda/mu
          PtrCoefFct coefLambdaOverMu = CoefFunction::Generate( mp_,  Global::REAL, CoefXprBinOp(mp_, bulkViscosity, shearViscosity, CoefXpr::OP_DIV ));
          surfVelocityConstraintMaxwellIntLV->SetBCoefFunctionOpB(coefLambdaOverMu);

          BiLinFormContext *surfVelocityConstraintMaxwellContextLV = new BiLinFormContext(surfVelocityConstraintMaxwellIntLV, STIFFNESS );
          surfVelocityConstraintMaxwellContextLV->SetEntities( actSDList, actSDList);
          surfVelocityConstraintMaxwellContextLV->SetFeFunctions( feFunctions_[LAGRANGE_MULT], feFunctions_[FLUIDMECH_VELOCITY]);
          feFunctions_[LAGRANGE_MULT]->AddEntityList( actSDList );
          feFunctions_[FLUIDMECH_VELOCITY]->AddEntityList( actSDList );
          surfVelocityConstraintMaxwellContextLV->SetCounterPart(true);
          assemble_->AddBiLinearForm( surfVelocityConstraintMaxwellContextLV );
        }
      }


      //========================================================================================
      // scalar velocity constraint
      //========================================================================================

      ReadRhsExcitation("scalarVelocityConstraint", feFunctions_[FLUIDMECH_VELOCITY]->GetResultInfo()->dofNames,
          ResultInfo::SCALAR, isComplex_, ent, kCoef, updatedGeo_, volumeRegions);

      // Define FE space etc.
      if ( useLagrangeMultScal_ ) {
        lagrangeMultFct_1 = feFunctions_[LAGRANGE_MULT_1];
        lagrangeMultSpace_1 = lagrangeMultFct_1->GetFeSpace();
      }

      for( UInt i = 0; i < ent.GetSize(); ++i ) {
        // get the volume region for defining the correct normal direction
        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volumeRegions[i]);
        std::set<RegionIdType> volRegion;
        volRegion.insert(aRegion);

        // check type of entitylist
        if (ent[i]->GetType() != EntityList::SURF_ELEM_LIST) {
          EXCEPTION("Velocity constraint BC must be defined on (surface) elements")
        }

        std::string regionName = ent[i]->GetName();
        shared_ptr<EntityList> actSDList = ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );

        // Get polynomial and integration order for Lagrange multiplier
        lagrangeMultSpace_1->SetRegionApproximation(aRegion, lagrangeMultPolyId_, lagrangeMultIntegId_);


        PtrParamNode scalVelConstNode = bcNode->GetByVal("scalarVelocityConstraint","name", regionName);
        ParamNodeList scalVelConstNodeList = scalVelConstNode->GetChildren();

        bool isNoPenetrationBC = false;
        for (UInt i = 0; i < scalVelConstNodeList.GetSize(); i++){
          std::string nodeName = scalVelConstNodeList[i]->GetName();
          if( nodeName=="noPenetration" ){
            isNoPenetrationBC = true;
          }
        }

        // setup the integrators

        // adding the additional velocity based surface integrators from integration by parts here causes some trouble, hence, we skip them
        // it is ok for this case since the contribution would be similar (in normal direction) and will be just "taken over" by the Lagrange multiplier
        AddSurfaceIntegratorFromPartialIntegration(i, volumeRegions, ent, "P1");

        if( isNoPenetrationBC ) {
          // conservation of momentum: \int_\Gamma_D v' \cdot n \lambda
          BiLinearForm * surfScalarVelocityConstraintIntVL = NULL;
          PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
          if(isComplex_) {
            if (dim_ == 2){
              surfScalarVelocityConstraintIntVL = new SurfaceABInt<Complex,Complex>(new IdentityOperator<FeH1,2,2>(),
                                                                                    new IdentityOperatorNormal<FeH1,2,1>(),
                                                                                    constOne, Complex(balanceOfMomentumSign_,0), volRegion, updatedGeo_ );
            } else {
              surfScalarVelocityConstraintIntVL = new SurfaceABInt<Complex,Complex>(new IdentityOperator<FeH1,3,3>(),
                                                                                    new IdentityOperatorNormal<FeH1,3,1>(),
                                                                                    constOne, Complex(balanceOfMomentumSign_,0), volRegion, updatedGeo_ );
            }
          } else {
            if (dim_ == 2){
              surfScalarVelocityConstraintIntVL = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                                                      new IdentityOperatorNormal<FeH1,2,1>(),
                                                                      constOne, balanceOfMomentumSign_, volRegion, updatedGeo_ );
            } else {
              surfScalarVelocityConstraintIntVL = new SurfaceABInt<>(new IdentityOperator<FeH1,3,3>(),
                                                                      new IdentityOperatorNormal<FeH1,3,1>(),
                                                                      constOne, balanceOfMomentumSign_, volRegion, updatedGeo_ );
            }
          }
          surfScalarVelocityConstraintIntVL->SetName("scalarVelocityConstraintIntVL");
          BiLinFormContext *surfScalarVelocityConstraintContextVL = new BiLinFormContext(surfScalarVelocityConstraintIntVL, STIFFNESS );
          surfScalarVelocityConstraintContextVL->SetEntities( actSDList, actSDList);
          surfScalarVelocityConstraintContextVL->SetFeFunctions( feFunctions_[FLUIDMECH_VELOCITY], feFunctions_[LAGRANGE_MULT_1]);
          feFunctions_[FLUIDMECH_VELOCITY]->AddEntityList( actSDList );
          feFunctions_[LAGRANGE_MULT_1]->AddEntityList( actSDList );
          assemble_->AddBiLinearForm( surfScalarVelocityConstraintContextVL );

          // constraint equation: \int_\Gamma_D \lambda' n \cdot v
          BiLinearForm * surfScalarVelocityConstraintIntLV = NULL;
          if(isComplex_) {
            if (dim_ == 2){
              surfScalarVelocityConstraintIntLV = new SurfaceABInt<Complex,Complex>(new IdentityOperatorNormal<FeH1,2,1>(),
                                                                                    new IdentityOperator<FeH1,2,2>(),
                                                                                    constOne, Complex(balanceOfMomentumSign_,0), volRegion, updatedGeo_ );
            } else {
              surfScalarVelocityConstraintIntLV = new SurfaceABInt<Complex,Complex>(new IdentityOperatorNormal<FeH1,3,1>(),
                                                                                    new IdentityOperator<FeH1,3,3>(),
                                                                                    constOne, Complex(balanceOfMomentumSign_,0), volRegion, updatedGeo_ );
            }
          } else {
            if (dim_ == 2){
              surfScalarVelocityConstraintIntLV = new SurfaceABInt<>(new IdentityOperatorNormal<FeH1,2,1>(),
                                                                      new IdentityOperator<FeH1,2,2>(),
                                                                      constOne, balanceOfMomentumSign_, volRegion, updatedGeo_ );
            } else {
              surfScalarVelocityConstraintIntLV = new SurfaceABInt<>(new IdentityOperatorNormal<FeH1,3,1>(),
                                                                      new IdentityOperator<FeH1,3,3>(),
                                                                      constOne, balanceOfMomentumSign_, volRegion, updatedGeo_ );
            }
          }
          surfScalarVelocityConstraintIntLV->SetName("surfScalarVelocityConstraintIntLV");
          BiLinFormContext *surfScalarVelocityConstraintContextLV = new BiLinFormContext(surfScalarVelocityConstraintIntLV, STIFFNESS );
          surfScalarVelocityConstraintContextLV->SetEntities( actSDList, actSDList);
          surfScalarVelocityConstraintContextLV->SetFeFunctions( feFunctions_[LAGRANGE_MULT_1], feFunctions_[FLUIDMECH_VELOCITY]);
          feFunctions_[FLUIDMECH_VELOCITY]->AddEntityList( actSDList );
          feFunctions_[LAGRANGE_MULT_1]->AddEntityList( actSDList );
          assemble_->AddBiLinearForm( surfScalarVelocityConstraintContextLV );
        }
      }
    }
  }

  void LinFlowPDE::DefineSolveStep() {
    solveStep_ = new StdSolveStep(*this);
  }

  void LinFlowPDE::InitTimeStepping() {

    GLMScheme * schemeV = new Trapezoidal(0.5);
    GLMScheme * schemeP = new Trapezoidal(0.5);
    // Important: Create a new time scheme just for the Lagrange multiplier unknowns, as otherwise the
	  // size of the vectors does not match!
    GLMScheme * schemeL = NULL;
    GLMScheme * schemeL_1 = NULL;
    if ( useLagrangeMultVec_ ) {
	    schemeL = new Trapezoidal(0.5);
    }
    if ( useLagrangeMultScal_ ) {
      schemeL_1 = new Trapezoidal(0.5);
    }

    shared_ptr<BaseTimeScheme> mySchemeV(new TimeSchemeGLM(schemeV, 0) );
    shared_ptr<BaseTimeScheme> mySchemeP(new TimeSchemeGLM(schemeP, 0) );
    feFunctions_[FLUIDMECH_VELOCITY]->SetTimeScheme(mySchemeV);
    feFunctions_[FLUIDMECH_PRESSURE]->SetTimeScheme(mySchemeP);
    if ( useLagrangeMultVec_ ) {
      shared_ptr<BaseTimeScheme> mySchemeL(new TimeSchemeGLM(schemeL, 0) );
      feFunctions_[LAGRANGE_MULT]->SetTimeScheme(mySchemeL);
    }
    if ( useLagrangeMultScal_ ) {
      shared_ptr<BaseTimeScheme> mySchemeL_1(new TimeSchemeGLM(schemeL_1, 0) );
      feFunctions_[LAGRANGE_MULT_1]->SetTimeScheme(mySchemeL_1);
    }
  }
  

  void LinFlowPDE::DefinePrimaryResults() {
    shared_ptr<BaseFeFunction> feFct = feFunctions_[FLUIDMECH_VELOCITY];

    // Check for subType
    StdVector<std::string> velDofNames;
    if ( ptGrid_->GetDim() == 3 ) {
      velDofNames = "x", "y", "z";
    } else { 
      if( ptGrid_->IsAxi() ) {
        velDofNames = "r", "z";
      } else {
        velDofNames = "x", "y";
      }
    }

    // === Primary result according to definition ===
    // PRESSURE
    shared_ptr<ResultInfo> pressure( new ResultInfo);
    pressure->resultType = FLUIDMECH_PRESSURE;
    pressure->dofNames = "";
    pressure->unit = MapSolTypeToUnit(pressure->resultType);

    pressure->definedOn = ResultInfo::NODE;
    pressure->entryType = ResultInfo::SCALAR;
    feFunctions_[FLUIDMECH_PRESSURE]->SetResultInfo(pressure);
    results_.Push_back( pressure );
    availResults_.insert( pressure );

    pressure->SetFeFunction(feFunctions_[FLUIDMECH_PRESSURE]);
    DefineFieldResult( feFunctions_[FLUIDMECH_PRESSURE], pressure);

    // VELOCITY
    shared_ptr<ResultInfo> velocity( new ResultInfo);
    velocity->resultType = FLUIDMECH_VELOCITY;
    velocity->dofNames = velDofNames;
    velocity->unit = MapSolTypeToUnit(velocity->resultType);

    velocity->definedOn = ResultInfo::NODE;
    velocity->entryType = ResultInfo::VECTOR;
    feFunctions_[FLUIDMECH_VELOCITY]->SetResultInfo(velocity);
    results_.Push_back( velocity );
    availResults_.insert( velocity );

    velocity->SetFeFunction(feFunctions_[FLUIDMECH_VELOCITY]);
    DefineFieldResult( feFunctions_[FLUIDMECH_VELOCITY], velocity );

    // MEAN VELOCITY
    CreateMeanFlowFunction(velDofNames);


    if ( useLagrangeMultVec_ ) {
      // === LAGRANGE MULTIPLIER ===
      // used for weakly enforced Dirichlet BC (e.g. Maxwell Slip BC)
      shared_ptr<BaseFeFunction> lagrangeMultFct = feFunctions_[LAGRANGE_MULT];
      shared_ptr<ResultInfo> lagrangeMult(new ResultInfo);
      lagrangeMult->resultType = LAGRANGE_MULT;
      lagrangeMult->dofNames = velDofNames;
      lagrangeMult->unit = "";
      lagrangeMult->definedOn = ResultInfo::NODE;
      lagrangeMult->entryType = ResultInfo::VECTOR;
      results_.Push_back( lagrangeMult );
      availResults_.insert( lagrangeMult );
      lagrangeMultFct->SetResultInfo( lagrangeMult );
      DefineFieldResult( lagrangeMultFct, lagrangeMult );
    }
      
    if ( useLagrangeMultScal_ ) {
      // === LAGRANGE MULTIPLIER 1 ===
      // used for weakly enforced Dirichlet BC in only one direction (e.g. noPenetration BC)
      shared_ptr<BaseFeFunction> lagrangeMultFct_1 = feFunctions_[LAGRANGE_MULT_1];
      shared_ptr<ResultInfo> lagrangeMult_1(new ResultInfo);
      lagrangeMult_1->resultType = LAGRANGE_MULT_1;
      lagrangeMult_1->dofNames = "";
      lagrangeMult_1->unit = "";
      lagrangeMult_1->definedOn = ResultInfo::NODE;
      lagrangeMult_1->entryType = ResultInfo::SCALAR;
      results_.Push_back( lagrangeMult_1 );
      availResults_.insert( lagrangeMult_1 );
      lagrangeMultFct_1->SetResultInfo( lagrangeMult_1 );
      DefineFieldResult( lagrangeMultFct_1, lagrangeMult_1 );
    }

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    idbcSolNameMap_.insert( std::pair<SolutionType, std::string>(FLUIDMECH_PRESSURE,"pressure") );
    idbcSolNameMap_.insert( std::pair<SolutionType, std::string>(FLUIDMECH_VELOCITY,"velocity") );
    if ( useLagrangeMultVec_ ) {
      idbcSolNameMap_.insert( std::pair<SolutionType, std::string>(LAGRANGE_MULT,"lagrangeMultiplier") );
    }
    if ( useLagrangeMultScal_ ) {
      idbcSolNameMap_.insert( std::pair<SolutionType, std::string>(LAGRANGE_MULT_1,"lagrangeMultiplier1") );
    }
  
    hdbcSolNameMap_[FLUIDMECH_PRESSURE] = "noPressure";
    hdbcSolNameMap_[FLUIDMECH_VELOCITY] = "noSlip";
  }
  
  
  void LinFlowPDE::DefinePostProcResults() {

    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    shared_ptr<BaseFeFunction> feFct = feFunctions_[FLUIDMECH_PRESSURE];

      StdVector<std::string> stressComponents;
      if( subType_ == "3d" ) {
        stressComponents = "xx", "yy", "zz", "yz", "xz", "xy";
      } else if( subType_ == "plane" ) {
        stressComponents = "xx", "yy", "xy";
      } else if( subType_ == "axi" ) {
        stressComponents = "rr", "zz", "rz", "phiphi";
      }
      StdVector<std::string > dispDofNames;
      dispDofNames = feFunctions_[FLUIDMECH_VELOCITY]->GetResultInfo()->dofNames;
      shared_ptr<BaseFeFunction> velFeFct = feFunctions_[FLUIDMECH_VELOCITY];

      // === FLUID-MECHANIC STRESS ===
      // mu ( grad(v)+grad(v)^T )
      shared_ptr<ResultInfo> stress(new ResultInfo);
      stress->resultType = FLUIDMECH_STRESS;
      stress->dofNames = stressComponents;
      stress->unit = MapSolTypeToUnit(FLUIDMECH_STRESS);
      stress->entryType = ResultInfo::TENSOR;
      stress->definedOn = ResultInfo::ELEMENT;
      stress->SetFeFunction(feFunctions_[FLUIDMECH_VELOCITY]);
      availResults_.insert( stress );
      shared_ptr<CoefFunctionFormBased> sigmaFunc;
      if( isComplex_ ) {
        sigmaFunc.reset(new CoefFunctionFlux<Complex>(velFeFct, stress));
      } else {
        sigmaFunc.reset(new CoefFunctionFlux<Double>(velFeFct, stress));
      }
      // we add a specific integrator name to ensure that only this integrator gets passed to the coefFunction during the assignment in the SinglePDE during FinalizePostProcResults
      sigmaFunc->SetIntegratorName("LinFlowStiffIntViscous");
      DefineFieldResult( sigmaFunc, stress );
      stiffFormCoefs_.insert(sigmaFunc);

      // === FLUID-MECHANIC ZERO PRESSURE ===
      // This is a dummy result in order to get the iterative coupling to recognize
      // both PDEs as iterative coupled although it's just a forward coupling
      shared_ptr<ResultInfo> pressureZero(new ResultInfo);
      pressureZero->resultType = FLUIDMECH_ZERO_PRESSURE;
      pressureZero->dofNames = "";
      pressureZero->unit = MapSolTypeToUnit(FLUIDMECH_ZERO_PRESSURE);
      pressureZero->entryType = ResultInfo::SCALAR;
      pressureZero->definedOn = ResultInfo::NODE;
      pressureZero->SetFeFunction(feFunctions_[FLUIDMECH_PRESSURE]);
      availResults_.insert( pressureZero );
      PtrCoefFct constZero = CoefFunction::Generate( mp_, Global::REAL, "0.0");

      DefineFieldResult( constZero, pressureZero );

      // === FLUID-MECHANIC STRAINRATE ===
      shared_ptr<ResultInfo> strain(new ResultInfo);
      strain->resultType = FLUIDMECH_STRAINRATE;
      strain->dofNames = stressComponents;
      strain->unit =  MapSolTypeToUnit(FLUIDMECH_STRAINRATE);;
      strain->entryType = ResultInfo::TENSOR;
      strain->definedOn = ResultInfo::ELEMENT;
      strain->SetFeFunction(feFunctions_[FLUIDMECH_VELOCITY]);
      availResults_.insert( strain );
      shared_ptr<CoefFunctionFormBased> strainFunc; // Coef function for result evaluation
      if( isComplex_ ) {
        strainFunc.reset(new CoefFunctionBOp<Complex>(velFeFct, strain));
      } else {
        strainFunc.reset(new CoefFunctionBOp<Double>(velFeFct, strain));
      }
      // we add a specific integrator name to ensure that only this integrator gets passed to the coefFunction during the assignment in the SinglePDE during FinalizePostProcResults
      strainFunc->SetIntegratorName("LinFlowStiffIntViscous"); // coefFunctionBOp acts on strain part
      DefineFieldResult( strainFunc, strain );
      // we use the same stiffness bdbInt as for the stress, but evaluating it with CoefFunctionBOp seems to ignore the material and return the strain
      stiffFormCoefs_.insert(strainFunc);


      // check for grid velocity
      bool hasGridVel = false;
      PtrParamNode bcsNode = myParam_->Get("movingMeshList", ParamNode::PASS );
      if( bcsNode ) {
        if( bcsNode->Has("movingMesh") ) {
          hasGridVel = true;
        }
      }

      // check if we have a grid velocity and define it
      // note: at the moment we use element results since nodal results are extremely slow
      // TODO: fix the nodal results and change the results to nodal ones
      if ( hasGridVel )  {
        // === FLUID-MECHANIC MESH VELOCITY ELEM ===
        shared_ptr<ResultInfo> meshVelocityElem(new ResultInfo);
        meshVelocityElem->resultType = FLUIDMECH_MESH_VELOCITY_ELEM;
        meshVelocityElem->dofNames = dispDofNames;
        meshVelocityElem->unit =  MapSolTypeToUnit(FLUIDMECH_MESH_VELOCITY_ELEM);;
        meshVelocityElem->entryType = ResultInfo::VECTOR;
        meshVelocityElem->definedOn = ResultInfo::ELEMENT;
        availResults_.insert( meshVelocityElem );
        DefineFieldResult( gridVelCoef_, meshVelocityElem );

        // === FLUID-MECHANIC MESH VELOCITY===
        shared_ptr<ResultInfo> meshVelocity(new ResultInfo);
        meshVelocity->resultType = FLUIDMECH_MESH_VELOCITY;
        meshVelocity->dofNames = dispDofNames;
        meshVelocity->unit =  MapSolTypeToUnit(FLUIDMECH_MESH_VELOCITY);;
        meshVelocity->entryType = ResultInfo::VECTOR;
        meshVelocity->definedOn = ResultInfo::NODE;
        availResults_.insert( meshVelocity );
        DefineFieldResult( gridVelCoef_, meshVelocity );
      }

      // === FLUID-MECHANIC TOTAL VELOCITY ELEM (for testing purposes) ===
      shared_ptr<ResultInfo> totalVelocityElem(new ResultInfo);
      totalVelocityElem->resultType = FLUIDMECH_TOTAL_VELOCITY_ELEM;
      totalVelocityElem->dofNames = dispDofNames;
      totalVelocityElem->unit =  MapSolTypeToUnit(FLUIDMECH_TOTAL_VELOCITY_ELEM);;
      totalVelocityElem->entryType = ResultInfo::VECTOR;
      totalVelocityElem->definedOn = ResultInfo::ELEMENT;
      availResults_.insert( totalVelocityElem );
      PtrCoefFct velocityFuncElem = this->GetCoefFct( FLUIDMECH_VELOCITY );
      PtrCoefFct totalVelCoefElem;

      if ( hasGridVel )  {
        //! Total velocity = flow velocity + grid velocity (for testing purposes)
        totalVelCoefElem =
                CoefFunction::Generate( mp_, part,
                CoefXprBinOp(mp_,velocityFuncElem,gridVelCoef_,CoefXpr::OP_ADD));
      }
      else {
        totalVelCoefElem = velocityFuncElem;
      }
      DefineFieldResult( totalVelCoefElem, totalVelocityElem );
      

      // === FLUID-MECHANIC VISCOUS STRESS ===
      // mu ( grad(v)+grad(v)^T ) + (lambda - 2/3 mu) div(v) I (compressible formulation)
      // mu ( grad(v)+grad(v)^T ) (incompressible formulation)
      shared_ptr<ResultInfo> stressVisc(new ResultInfo);
      stressVisc->resultType = FLUIDMECH_VISC_STRESS;
      stressVisc->dofNames = stressComponents;
      stressVisc->unit = MapSolTypeToUnit(FLUIDMECH_VISC_STRESS);
      stressVisc->entryType = ResultInfo::TENSOR;
      stressVisc->definedOn = ResultInfo::ELEMENT;
      stressVisc->SetFeFunction(feFunctions_[FLUIDMECH_VELOCITY]);
      availResults_.insert( stressVisc );
      PtrCoefFct stressViscCoef;
      // we define the stressViscCoef in the next section since we have to differ between compressible / incompressible formulation


      PtrCoefFct constZeroSolDepend = CoefFunction::Generate( mp_, part, "0.0");
      constZeroSolDepend->SetSolDependent(); //we use this coefFunction later on to fill in the gaps for coefFunctionDiagTensorFromScalar - but here we need similar depend types, hence we use this trick
      if ( isCompressible_ ) {
        // === FLUID-MECHANIC STRESS (COMPRESSIBLE) ===
        // (lambda - 2/3 mu) div(v) I
        shared_ptr<ResultInfo> stressComp(new ResultInfo);
        stressComp->resultType = FLUIDMECH_COMP_STRESS;
        stressComp->dofNames = stressComponents;
        stressComp->unit = MapSolTypeToUnit(FLUIDMECH_COMP_STRESS);
        stressComp->entryType = ResultInfo::TENSOR;
        stressComp->definedOn = ResultInfo::ELEMENT;
        stressComp->SetFeFunction(feFunctions_[FLUIDMECH_VELOCITY]);
        availResults_.insert( stressComp );
        shared_ptr<CoefFunctionFormBased> stressCompFunc;
        if( isComplex_ ) {
          stressCompFunc.reset(new CoefFunctionFlux<Complex>(velFeFct, stressComp));
        } else {
          stressCompFunc.reset(new CoefFunctionFlux<Double>(velFeFct, stressComp));
        }
        // we add a specific integrator name to ensure that only this integrator gets passed to the coefFunction during the assignment in the SinglePDE during FinalizePostProcResults
        stressCompFunc->SetIntegratorName("LinFlowStiffIntBulkViscous"); // coefFunctionFlux acts on strain part
        StdVector<PtrCoefFct> stressCompDiagValues;
        if( subType_ == "axi" ) {
          stressCompDiagValues = StdVector<PtrCoefFct>(2*dim_);
        } else {
          stressCompDiagValues = StdVector<PtrCoefFct>(dim_);
        }
        for(UInt i = 0; i < dim_; i++){
          stressCompDiagValues[i] = stressCompFunc;
        }
        if( subType_ == "axi" ) {
          stressCompDiagValues[2] = constZeroSolDepend;
          stressCompDiagValues[3] = stressCompFunc;
        }
        shared_ptr<CoefFunction> stressCompTensor (new CoefFunctionDiagTensorFromScalar(stressCompDiagValues,subType_));
        DefineFieldResult( stressCompTensor, stressComp );
        stiffFormCoefs_.insert(stressCompFunc);


        // insert the missing coefFunction for the viscous stress result
        stressViscCoef =
              CoefFunction::Generate( mp_, part,
              CoefXprBinOp(mp_,sigmaFunc,stressCompFunc,CoefXpr::OP_ADD));
      } else {
        // insert the missing coefFunction for the viscous stress result
        stressViscCoef = sigmaFunc;
      }
      DefineFieldResult( stressViscCoef, stressVisc );


      // === FLUID-MECHANIC PRESSURE TENSOR (needed for further postprocessing) ===
      // p I
      shared_ptr<ResultInfo> presTens(new ResultInfo);
      presTens->resultType = FLUIDMECH_PRES_TENS;
      presTens->dofNames = stressComponents;
      presTens->unit = MapSolTypeToUnit(FLUIDMECH_PRES_TENS);
      presTens->entryType = ResultInfo::TENSOR;
      presTens->definedOn = ResultInfo::ELEMENT;
      presTens->SetFeFunction(feFunctions_[FLUIDMECH_PRESSURE]);
      availResults_.insert( presTens );
      PtrCoefFct presFnc = this->GetCoefFct( FLUIDMECH_PRESSURE );
      StdVector<PtrCoefFct> presTensDiagValues;
      PtrCoefFct zeroPres = CoefFunction::Generate( mp_,  part,
          CoefXprBinOp(mp_, presFnc, CoefFunction::Generate( mp_, part, "0"), CoefXpr::OP_MULT));
      if( subType_ == "axi" ) {
        presTensDiagValues = StdVector<PtrCoefFct>(2*dim_);
      } else {
        presTensDiagValues = StdVector<PtrCoefFct>(dim_);
      }
      for(UInt i = 0; i < dim_; i++){
        presTensDiagValues[i] = presFnc;
      }
      
      if( subType_ == "axi" ) {
        presTensDiagValues[2] = zeroPres; //evaluates to zero but is needed to get the same depend type
        presTensDiagValues[3] = presFnc;
      }
      shared_ptr<CoefFunction> presTensCoef (new CoefFunctionDiagTensorFromScalar(presTensDiagValues,subType_));
      DefineFieldResult( presTensCoef, presTens );


      // === FLUID-MECHANIC TOTAL STRESS ===
      // -p I + mu ( grad(v)+grad(v)^T ) + (lambda - 2/3 mu) div(v) I (compressible formulation)
      // -p I + mu ( grad(v)+grad(v)^T ) (incompressible formulation)
      shared_ptr<ResultInfo> stressTotal(new ResultInfo);
      stressTotal->resultType = FLUIDMECH_TOTAL_STRESS;
      stressTotal->dofNames = stressComponents;
      stressTotal->unit = MapSolTypeToUnit(FLUIDMECH_TOTAL_STRESS);
      stressTotal->entryType = ResultInfo::TENSOR;
      stressTotal->definedOn = ResultInfo::ELEMENT;
      availResults_.insert( stressTotal );
      PtrCoefFct stressTotalCoef;
      stressTotalCoef =
              CoefFunction::Generate( mp_, part,
              CoefXprBinOp(mp_,stressViscCoef,presTensCoef,CoefXpr::OP_SUB));
      DefineFieldResult( stressTotalCoef, stressTotal );


      // === FLUID-MECHANIC SURFACE TRACTION ===
      shared_ptr<ResultInfo> surfaceTractionInfo;
      shared_ptr<CoefFunctionSurf> surfaceTractionFct;
      surfaceTractionInfo.reset(new ResultInfo);
      surfaceTractionInfo->resultType = FLUIDMECH_SURFACE_TRACTION;
      surfaceTractionInfo->dofNames = dispDofNames;
      surfaceTractionInfo->unit = MapSolTypeToUnit(FLUIDMECH_SURFACE_TRACTION);
      surfaceTractionInfo->entryType = ResultInfo::VECTOR;
      surfaceTractionInfo->definedOn = ResultInfo::SURF_ELEM;

      surfaceTractionFct.reset(new CoefFunctionSurf(true, 1.0, surfaceTractionInfo));
      DefineFieldResult(surfaceTractionFct, surfaceTractionInfo);
      surfCoefFcts_[surfaceTractionFct] = stressTotalCoef;

      // === FLUID-MECHANIC REACTION FORCE (= integral of surface traction over the surface region ) ===
      shared_ptr<ResultInfo> reactionForceInfo;
      reactionForceInfo.reset(new ResultInfo);
      reactionForceInfo->resultType = FLUIDMECH_FORCE;
      reactionForceInfo->dofNames = dispDofNames;
      reactionForceInfo->unit = MapSolTypeToUnit(FLUIDMECH_FORCE);
      reactionForceInfo->entryType = ResultInfo::VECTOR;
      reactionForceInfo->definedOn = ResultInfo::SURF_REGION;
      // Integrate surface traction
      shared_ptr<ResultFunctor> reactionForceFct;
      if(isComplex_)
          reactionForceFct.reset(new ResultFunctorIntegrate<Complex>(surfaceTractionFct, feFct, reactionForceInfo));
      else
          reactionForceFct.reset(new ResultFunctorIntegrate<Double>(surfaceTractionFct, feFct, reactionForceInfo));
      resultFunctors_[FLUIDMECH_FORCE] = reactionForceFct;
      availResults_.insert(reactionForceInfo);

      // === FLUID-MECHANIC VISCOUS DISSIPATION POWER DENSITY DIVERGENCE : 1/2(\lambda -2/3 \mu)(Div v)^2===
      shared_ptr<ResultInfo> vdpd1(new ResultInfo);
      vdpd1->resultType = FLUIDMECH_VISCOUS_DISS_POWER_DENS_DIV;
      vdpd1->dofNames = "";
      vdpd1->unit =  MapSolTypeToUnit(FLUIDMECH_VISCOUS_DISS_POWER_DENS_DIV);
      vdpd1->entryType = ResultInfo::SCALAR;
      vdpd1->definedOn = ResultInfo::ELEMENT;
      vdpd1->SetFeFunction(feFunctions_[FLUIDMECH_VELOCITY]);
      availResults_.insert(vdpd1);  
      shared_ptr<CoefFunctionFormBased> vdpd1Func;    
      if( isComplex_ ) {
        vdpd1Func.reset(new CoefFunctionBdBKernel<Complex>(velFeFct, balanceOfMomentumSign_*0.5));
      } else {
        vdpd1Func.reset(new CoefFunctionBdBKernel<Double>(velFeFct, balanceOfMomentumSign_));
      }
      // we add a specific integrator name to ensure that only this integrator gets passed to the coefFunction during the assignment in the SinglePDE during FinalizePostProcResults
      vdpd1Func->SetIntegratorName("LinFlowStiffIntBulkViscous"); // coefFunctionBOp acts on strain part
      DefineFieldResult( vdpd1Func, vdpd1 );
      stiffFormCoefs_.insert( vdpd1Func );

      // === FLUID-MECHANIC VISCOUS DISSIPATION POWER DENSITY TOTAL STRAIN : 1/2 \mu (Grad v: Grad v + Grad v: (Grad v)^T )===
      shared_ptr<ResultInfo> vdpd2(new ResultInfo);
      vdpd2->resultType = FLUIDMECH_VISCOUS_DISS_POWER_DENS_STRAIN;
      vdpd2->dofNames = "";
      vdpd2->unit =  MapSolTypeToUnit(FLUIDMECH_VISCOUS_DISS_POWER_DENS_STRAIN);
      vdpd2->entryType = ResultInfo::SCALAR;
      vdpd2->definedOn = ResultInfo::ELEMENT;
      vdpd2->SetFeFunction(feFunctions_[FLUIDMECH_VELOCITY]);
      availResults_.insert( vdpd2 );
      shared_ptr<CoefFunctionFormBased> vdpd2Func;
      if( isComplex_ ) {
        vdpd2Func.reset(new CoefFunctionBdBKernel<Complex>(velFeFct, balanceOfMomentumSign_*0.5));
      } else {
        vdpd2Func.reset(new CoefFunctionBdBKernel<Double>(velFeFct, balanceOfMomentumSign_));
      }
      // we add a specific integrator name to ensure that only this integrator gets passed to the coefFunction during the assignment in the SinglePDE during FinalizePostProcResults
      vdpd2Func->SetIntegratorName("LinFlowStiffIntViscous");
      DefineFieldResult( vdpd2Func, vdpd2 );
      stiffFormCoefs_.insert(vdpd2Func);

      // === VISCOUS DISSIPATION POWER DENSITY - TOTAL : 1/2(\lambda -2/3 \mu)(Div v) +1/2 \mu (Grad v: Grad v + Grad v: (Grad v)^T ) ===
      // The actual summation of Term1 and Term2
      shared_ptr<ResultInfo> vdpd(new ResultInfo);
      vdpd->resultType = FLUIDMECH_VISCOUS_DISS_POWER_DENS;
      vdpd->dofNames = "";
      vdpd->unit =  MapSolTypeToUnit(FLUIDMECH_VISCOUS_DISS_POWER_DENS);
      vdpd->entryType = ResultInfo::SCALAR;
      vdpd->definedOn = ResultInfo::ELEMENT;
      shared_ptr<CoefFunctionMulti> vdpdFunc(new CoefFunctionMulti(CoefFunction::SCALAR,dim_,1, isComplex_));
      vdpd->SetFeFunction(feFunctions_[FLUIDMECH_VELOCITY]);
      DefineFieldResult( vdpdFunc, vdpd );

      // Summation of FLUIDMECH_VISCOUS_DISS_POWER_DENS_DIV and FLUIDMECH_VISCOUS_DISS_POWER_DENS_STRAIN
      shared_ptr<CoefFunctionMulti> vdpdCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[FLUIDMECH_VISCOUS_DISS_POWER_DENS]);
      StdVector<RegionIdType>::iterator regIt = regions_.Begin();
      regIt = regions_.Begin();
      for( ; regIt != regions_.End(); ++regIt ){
        std::string regionName = ptGrid_->GetRegion().ToString(*regIt);
        PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
        PtrCoefFct h = CoefFunction::Generate( mp_, part, CoefXprBinOp( mp_,
                                               GetCoefFct(FLUIDMECH_VISCOUS_DISS_POWER_DENS_DIV),
                                               GetCoefFct(FLUIDMECH_VISCOUS_DISS_POWER_DENS_STRAIN), CoefXpr::OP_ADD ) );
        vdpdCoef->AddRegion(*regIt, h);
      }

      // === FLUID-MECHANIC DISSIPATION POWER (= integral of dissipation power density over the volume ) ===
      shared_ptr<ResultInfo> vdp;
      vdp.reset(new ResultInfo);
      vdp->resultType = FLUIDMECH_VISCOUS_DISS_POWER;
      vdp->dofNames = "";
      vdp->unit = MapSolTypeToUnit(FLUIDMECH_VISCOUS_DISS_POWER);
      vdp->entryType = ResultInfo::SCALAR;
      vdp->definedOn = ResultInfo::REGION;
      // Integrate surface traction
      shared_ptr<ResultFunctor> vdpFct;
      if(isComplex_)
          vdpFct.reset(new ResultFunctorIntegrate<Complex>(vdpdCoef, feFct, vdp));
      else
          vdpFct.reset(new ResultFunctorIntegrate<Double>(vdpdCoef, feFct, vdp));
      resultFunctors_[FLUIDMECH_VISCOUS_DISS_POWER] = vdpFct;
      availResults_.insert(vdp);


      // === FLUID-MECHAINC INTENSITY SIMPLE (PRESSURE ONLY) ===
      // Intensity I = p  conj(v)
      shared_ptr<ResultInfo> intensitySimple;
      intensitySimple.reset(new ResultInfo);
      intensitySimple->resultType = FLUIDMECH_INTENSITY_PRESSURE_ONLY;
      intensitySimple->dofNames = dispDofNames;
      intensitySimple->unit = MapSolTypeToUnit(FLUIDMECH_INTENSITY_PRESSURE_ONLY);
      intensitySimple->entryType = ResultInfo::VECTOR;
      intensitySimple->definedOn = ResultInfo::ELEMENT;
      PtrCoefFct intensFctSimple;
      intensFctSimple = 
          CoefFunction::Generate( mp_, part,
                                 CoefXprBinOp(mp_, feFunctions_[FLUIDMECH_PRESSURE], feFunctions_[FLUIDMECH_VELOCITY], CoefXpr::OP_MULT_CONJ ) );
      DefineFieldResult(intensFctSimple, intensitySimple);


      // === FLUID-MECHAINC INTENSITY ===
      // Intensity I = sigma \cdot conj(v)
      shared_ptr<ResultInfo> intensity;
      intensity.reset(new ResultInfo);
      intensity->resultType = FLUIDMECH_INTENSITY;
      intensity->dofNames = dispDofNames;
      intensity->unit = MapSolTypeToUnit(FLUIDMECH_INTENSITY);
      intensity->entryType = ResultInfo::VECTOR;
      intensity->definedOn = ResultInfo::ELEMENT;
      PtrCoefFct intensFct;
      PtrCoefFct velFnc = this->GetCoefFct( FLUIDMECH_VELOCITY );
      // define temporary function, without the -1 sign
      PtrCoefFct intensTmp = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, stressTotalCoef, velFnc, CoefXpr::OP_MULT_VOIGT_TENSOR_VEC_CONJ));
      intensFct = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_,  "-1.0", intensTmp , CoefXpr::OP_MULT));
      DefineFieldResult( intensFct, intensity );

      
      // === FLUID-MECHAINC SURFACE INTENSITY SIMPLE (PRESSURE ONLY) ===
      // Intensity I = p  conj(v) (evaluated at the boundary)
      shared_ptr<ResultInfo> surfIntensitySimple;
      surfIntensitySimple.reset(new ResultInfo);
      surfIntensitySimple->resultType = FLUIDMECH_SURFINTENSITY_PRESSURE_ONLY;
      surfIntensitySimple->dofNames = dispDofNames;
      surfIntensitySimple->unit = MapSolTypeToUnit(FLUIDMECH_SURFINTENSITY_PRESSURE_ONLY);
      surfIntensitySimple->entryType = ResultInfo::VECTOR;
      surfIntensitySimple->definedOn = ResultInfo::SURF_ELEM;
      shared_ptr<CoefFunctionSurf> sIntensSimple;
      sIntensSimple.reset(new CoefFunctionSurf(false, 1.0, surfIntensitySimple));
      DefineFieldResult(sIntensSimple, surfIntensitySimple);
      surfCoefFcts_[sIntensSimple] = intensFctSimple;

      
      // === FLUID-MECHAINC SURFACE INTENSITY ===
      // Intensity I = sigma \cdot conj(v) (evaluated at the boundary)
      shared_ptr<ResultInfo> surfIntensity;
      surfIntensity.reset(new ResultInfo);
      surfIntensity->resultType = FLUIDMECH_SURFINTENSITY;
      surfIntensity->dofNames = dispDofNames;
      surfIntensity->unit = MapSolTypeToUnit(FLUIDMECH_SURFINTENSITY);
      surfIntensity->entryType = ResultInfo::VECTOR;
      surfIntensity->definedOn = ResultInfo::SURF_ELEM;
      shared_ptr<CoefFunctionSurf> sIntens;
      sIntens.reset(new CoefFunctionSurf(false, 1.0, surfIntensity));
      DefineFieldResult(sIntens, surfIntensity);
      surfCoefFcts_[sIntens] = intensFct;


      // === FLUID-MECHAINC NORMAL_INTENSITY SIMPLE (PRESSURE ONLY) ===
      // Normal intensity I_n = (p  conj(v)) \cdot n
      shared_ptr<ResultInfo> intensNormalSimple;
      intensNormalSimple.reset(new ResultInfo);
      intensNormalSimple->resultType = FLUIDMECH_NORMAL_INTENSITY_PRESSURE_ONLY;
      intensNormalSimple->dofNames = "";
      intensNormalSimple->unit = MapSolTypeToUnit(FLUIDMECH_NORMAL_INTENSITY_PRESSURE_ONLY);
      intensNormalSimple->entryType = ResultInfo::SCALAR;
      intensNormalSimple->definedOn = ResultInfo::SURF_ELEM;  
      shared_ptr<CoefFunctionSurf> sNormIntensSimple;
      sNormIntensSimple.reset(new CoefFunctionSurf(true, 1.0, intensNormalSimple));
      DefineFieldResult( sNormIntensSimple, intensNormalSimple );
      surfCoefFcts_[sNormIntensSimple] = intensFctSimple;


      // === FLUID-MECHAINC NORMAL_INTENSITY ===
      // Normal intensity I_n = (\sigma \cdot conj(v)) \cdot n
      shared_ptr<ResultInfo> intensNormal;
      intensNormal.reset(new ResultInfo);
      intensNormal->resultType = FLUIDMECH_NORMAL_INTENSITY;
      intensNormal->dofNames = "";
      intensNormal->unit = MapSolTypeToUnit(FLUIDMECH_NORMAL_INTENSITY);
      intensNormal->entryType = ResultInfo::SCALAR;
      intensNormal->definedOn = ResultInfo::SURF_ELEM;  
      shared_ptr<CoefFunctionSurf> sNormIntens;
      sNormIntens.reset(new CoefFunctionSurf(true, 1.0, intensNormal));
      DefineFieldResult( sNormIntens, intensNormal );
      surfCoefFcts_[sNormIntens] = intensFct;


      // === FLUID-MECHAINC POWER SIMPLE (PRESSURE ONLY) ===
      // Power P = \int_Gamma (p  conj(v)) \cdot n dGamma
      // Integrate normal intensity
      shared_ptr<ResultInfo> powerSimple;
      powerSimple.reset(new ResultInfo);
      powerSimple->resultType = FLUIDMECH_POWER_PRESSURE_ONLY;
      powerSimple->dofNames = "";
      powerSimple->unit = MapSolTypeToUnit(FLUIDMECH_POWER_PRESSURE_ONLY);
      powerSimple->entryType = ResultInfo::SCALAR;
      powerSimple->definedOn = ResultInfo::SURF_REGION;
      shared_ptr<ResultFunctor> powerFctSimple;
      if( isComplex_ ) {
        powerFctSimple.reset(new ResultFunctorIntegrate<Complex>(sNormIntensSimple, 
                                                            feFct, powerSimple ) );
      } else {
        powerFctSimple.reset(new ResultFunctorIntegrate<Double>(sNormIntensSimple, 
                                                           feFct, powerSimple ) );
      }
      resultFunctors_[FLUIDMECH_POWER_PRESSURE_ONLY] = powerFctSimple;
      availResults_.insert(powerSimple);


      // === FLUID-MECHAINC POWER ===
      // Power P = \int_Gamma (sigma \cdot conj(v)) \cdot n dGamma
      // Integrate normal intensity
      shared_ptr<ResultInfo> power;
      power.reset(new ResultInfo);
      power->resultType = FLUIDMECH_POWER;
      power->dofNames = "";
      power->unit = MapSolTypeToUnit(FLUIDMECH_POWER);
      power->entryType = ResultInfo::SCALAR;
      power->definedOn = ResultInfo::SURF_REGION;
      shared_ptr<ResultFunctor> powerFct;
      if( isComplex_ ) {
        powerFct.reset(new ResultFunctorIntegrate<Complex>(sNormIntens, 
                                                            feFct, power ) );
      } else {
        powerFct.reset(new ResultFunctorIntegrate<Double>(sNormIntens, 
                                                           feFct, power ) );
      }
      resultFunctors_[FLUIDMECH_POWER] = powerFct;
      availResults_.insert(power);

      // === FLUID-MECHAINC VELOCITY ===
      // Normal velocity v_{\mathrm{lf,n}} = \bm{v} \cdot \bm{n}
      shared_ptr<ResultInfo> velNormal(new ResultInfo);
      velNormal->resultType = FLUIDMECH_NORMAL_VELOCITY;
      velNormal->dofNames = "";
      velNormal->unit = MapSolTypeToUnit(FLUIDMECH_NORMAL_VELOCITY);
      velNormal->entryType = ResultInfo::SCALAR;
      velNormal->definedOn = ResultInfo::SURF_ELEM;
      shared_ptr<CoefFunctionSurf> velFncNormal;
      velFncNormal.reset(new CoefFunctionSurf(true, 1.0, velNormal));
      DefineFieldResult(velFncNormal, velNormal);
      surfCoefFcts_[velFncNormal] = velFnc;
      
      
      // === FLUID-MECHAINC SURFIMPEDANCE ===
      // Surface impedance Z_{\mathrm{lf,surf}} = v_{\mathrm{lf,n}} / p
      shared_ptr<ResultInfo> surfImpedance(new ResultInfo);
      surfImpedance->resultType = FLUIDMECH_SURFIMPEDANCE;
      surfImpedance->dofNames = "";
      surfImpedance->unit = MapSolTypeToUnit(FLUIDMECH_SURFIMPEDANCE);
      surfImpedance->entryType = ResultInfo::SCALAR;
      surfImpedance->definedOn = ResultInfo::SURF_ELEM;
      PtrCoefFct surfImpFct = CoefFunction::Generate(
          mp_, part, CoefXprBinOp(mp_, presFnc, velFncNormal, CoefXpr::OP_DIV));
      DefineFieldResult(surfImpFct, surfImpedance);
      
      // === FLUID-MECHAINC IMPEDANCE ===
      // Impedance Z_{\mathrm{lf}} = \int_{\mathrm{\Gamma}} v_{\mathrm{lf,n}} / p dGamma
      shared_ptr<ResultInfo> impedance(new ResultInfo);
      impedance->resultType = FLUIDMECH_IMPEDANCE;
      impedance->dofNames = "";
      impedance->unit = MapSolTypeToUnit(FLUIDMECH_IMPEDANCE);
      impedance->entryType = ResultInfo::SCALAR;
      impedance->definedOn = ResultInfo::SURF_REGION;
      // TODO area weighted average
      shared_ptr<ResultFunctor> impedanceFct;
      if( isComplex_ ) {
        impedanceFct.reset(
            new ResultFunctorIntegrate<Complex>(surfImpFct, feFct, impedance));
        dynamic_pointer_cast<ResultFunctorIntegrate<Complex>>(impedanceFct)
            ->SetAveraged(true);
      } else {
        impedanceFct.reset(
            new ResultFunctorIntegrate<Double>(surfImpFct, feFct, impedance));
        dynamic_pointer_cast<ResultFunctorIntegrate<Double>>(impedanceFct)
            ->SetAveraged(true);
      }
      resultFunctors_[FLUIDMECH_IMPEDANCE] = impedanceFct;
      availResults_.insert(impedance);

    }
  
  void LinFlowPDE::FinalizePostProcResults() {

    SinglePDE::FinalizePostProcResults();
    meanFlowFeFct_->Finalize();
    meanFlowFeFct_->ApplyExternalData();

    shared_ptr<CoefFunctionFormBased> sigmaFunc = 
        dynamic_pointer_cast<CoefFunctionFormBased>(
          fieldCoefs_[FLUIDMECH_STRESS]
          );
    shared_ptr<CoefFunctionFormBased> strainFunc = 
        dynamic_pointer_cast<CoefFunctionFormBased>(
          fieldCoefs_[FLUIDMECH_STRAINRATE]
          );
  }

  
  std::map<SolutionType, shared_ptr<FeSpace> >
  LinFlowPDE::CreateFeSpaces(const std::string& formulation,
                                   PtrParamNode infoNode) 
  {

    // Set sign of integrators
    // The originally implemented sign was 1.0, this was globally changed to
    // -1.0 except for the case of the eigenValue analysis type
    if (analysistype_ == AnalysisType::EIGENVALUE) {
      if( isHeatPDECoupled_ && isCouplingFormulationSymmetric_ ) {
          EXCEPTION("Eigenvalue analysis type is currently not implemented for symmetric LinFlowHeatCoupling!" \
            "\nPlease use the asymetric LinFlowHeatCoupling or consider implementing it for the symmetric LinFlowHeatCoupling.");
      }
      balanceOfMomentumSign_ = 1.0;
    } else {
      balanceOfMomentumSign_ = -1.0;
    }

    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;

    if(formulation == "default" || formulation == "H1"){
      PtrParamNode spaceNode;
      spaceNode = infoNode->Get(SolutionTypeEnum.ToString(FLUIDMECH_VELOCITY));

      crSpaces[FLUIDMECH_VELOCITY] =
        FeSpace::CreateInstance(myParam_,spaceNode,FeSpace::H1, ptGrid_);
      crSpaces[FLUIDMECH_VELOCITY]->Init(solStrat_);

      // Create same FeSpace for mean velocity as for the unknown velocity.
      meanFlowFeSpace_ =
        FeSpace::CreateInstance(myParam_,spaceNode,FeSpace::H1, ptGrid_);
      meanFlowFeSpace_->Init(solStrat_);

      spaceNode = infoNode->Get(SolutionTypeEnum.ToString(FLUIDMECH_PRESSURE));

      crSpaces[FLUIDMECH_PRESSURE] =
        FeSpace::CreateInstance(myParam_,spaceNode,FeSpace::H1, ptGrid_);
      crSpaces[FLUIDMECH_PRESSURE]->Init(solStrat_);

      if ( useLagrangeMultVec_ ) {
        // Create space for the lagrange multiplier
        spaceNode = infoNode->Get(SolutionTypeEnum.ToString(LAGRANGE_MULT));
        crSpaces[LAGRANGE_MULT] = FeSpace::CreateInstance(myParam_, spaceNode, FeSpace::H1, ptGrid_);
        crSpaces[LAGRANGE_MULT]->Init(solStrat_);
      }
      if ( useLagrangeMultScal_ ) {
        // Create space for the lagrange multiplier 1
        spaceNode = infoNode->Get(SolutionTypeEnum.ToString(LAGRANGE_MULT_1));
        crSpaces[LAGRANGE_MULT_1] = FeSpace::CreateInstance(myParam_, spaceNode, FeSpace::H1, ptGrid_);
        crSpaces[LAGRANGE_MULT_1]->Init(solStrat_);
      }

    }else{
      EXCEPTION("The formulation " << formulation << 
                "of fluid perturbed PDE is not known!");
    }
    return crSpaces;
  }


  void LinFlowPDE::CreateMeanFlowFunction(StdVector<std::string> dofNames)
  {
    //// === MEAN FLUIDMECH VELOCITY ===
    shared_ptr<ResultInfo> flowvelocity( new ResultInfo );
    flowvelocity->resultType = MEAN_FLUIDMECH_VELOCITY;
    flowvelocity->dofNames = dofNames;
    flowvelocity->unit = "m/s";

    flowvelocity->definedOn = ResultInfo::NODE;
    flowvelocity->entryType = ResultInfo::VECTOR;
    results_.Push_back( flowvelocity );
    availResults_.insert( flowvelocity );
    
    meanVelocityCoef_.reset(
      new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1, isComplex_ )
      );

    // Since we also need the derivative of the mean flow velocity, we need to
    // create a FeFunction from the CoefFunction.
    if(isComplex_) {
      meanFlowFeFct_.reset( new FeFunction<Complex>(domain_->GetMathParser()) );
    }
    else {
      meanFlowFeFct_.reset( new FeFunction<Double>(domain_->GetMathParser()) );    
    }    
    flowvelocity->SetFeFunction(meanFlowFeFct_);
    meanFlowFeSpace_->AddFeFunction(meanFlowFeFct_);
    meanFlowFeFct_->SetFeSpace( meanFlowFeSpace_ );
    meanFlowFeFct_->SetPDE( this );
    meanFlowFeFct_->SetGrid(ptGrid_);
    meanFlowFeFct_->SetResultInfo(flowvelocity);
    meanFlowFeFct_->SetFctId(PSEUDO_FCT_ID);

    DefineFieldResult(meanVelocityCoef_, flowvelocity);
  }


  void LinFlowPDE::AddSurfaceIntegratorFromPartialIntegration(UInt i, StdVector<std::string> volumeRegions,
                                                              StdVector<shared_ptr<EntityList> > ent, std::string integrator)
  {
    // Get FESpace and FeFunction of fluid velocity
    shared_ptr<BaseFeFunction> velFct = feFunctions_[FLUIDMECH_VELOCITY];
    shared_ptr<FeSpace> velSpace = velFct->GetFeSpace();

    // Get FESpace and FeFunction of fluid pressure
    shared_ptr<BaseFeFunction> presFct = feFunctions_[FLUIDMECH_PRESSURE];
    shared_ptr<FeSpace> presSpace = presFct->GetFeSpace();
    
    // get the volume region for defining the correct normal direction
    RegionIdType aRegion = ptGrid_->GetRegion().Parse(volumeRegions[i]);
    std::set<RegionIdType> volRegion;
    volRegion.insert(aRegion);
    std::string regionName = ent[i]->GetName();
    shared_ptr<EntityList> actSDList = ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );

    // define coefFunctions for the material
    PtrCoefFct shearViscosity = materials_[aRegion]->GetScalCoefFnc(FLUID_DYNAMIC_VISCOSITY, Global::REAL); // mu
    PtrCoefFct bulkViscosity;
    if ( isCompressible_ ) { // we need to subtract 2*mu/3 Div(v)I and add the bulk part lambda*Div(v)I
      bulkViscosity = materials_[aRegion]->GetScalCoefFnc(FLUID_BULK_VISCOSITY, Global::REAL);
    }
    PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
    PtrCoefFct constTwo = CoefFunction::Generate( mp_, Global::REAL, "2.0");

    PtrCoefFct shearViscosityDouble = CoefFunction::Generate( mp_,  Global::REAL,
        CoefXprBinOp(mp_, shearViscosity, CoefFunction::Generate( mp_, Global::REAL, "2"), CoefXpr::OP_MULT));
    PtrCoefFct coefZero = CoefFunction::Generate( mp_, Global::REAL, "0");
    StdVector<PtrCoefFct> tensorComponents(dim_ == 2 ? (isaxi_ == true ? 16 : 9) : 36);
    tensorComponents.Init(coefZero);
    PtrCoefFct coefBB;
    if( dim_ == 2 ) {
      if( isaxi_ ) {
        tensorComponents[0] = shearViscosityDouble;
        tensorComponents[5] = shearViscosityDouble;
        tensorComponents[10] = shearViscosity;
        tensorComponents[15] = shearViscosityDouble;
        coefBB = CoefFunction::Generate(mp_,Global::REAL,4,4,tensorComponents);
      } else if( subType_ == "plane" ) {
        tensorComponents[0] = shearViscosityDouble;
        tensorComponents[4] = shearViscosityDouble;
        tensorComponents[8] = shearViscosity;
        coefBB = CoefFunction::Generate(mp_,Global::REAL,3,3,tensorComponents);
      }
    } else {
      tensorComponents[0] = shearViscosityDouble;
      tensorComponents[7] = shearViscosityDouble;
      tensorComponents[14] = shearViscosityDouble;
      tensorComponents[21] = shearViscosity;
      tensorComponents[28] = shearViscosity;
      tensorComponents[35] = shearViscosity;
      coefBB = CoefFunction::Generate(mp_,Global::REAL,6,6,tensorComponents);
    }

    PtrCoefFct coefDivDiv;
    if ( isCompressible_ ) {
      coefDivDiv = CoefFunction::Generate( mp_,  Global::REAL,
        CoefXprBinOp(mp_,
            bulkViscosity,
            CoefXprBinOp(mp_,shearViscosityDouble,CoefFunction::Generate( mp_, Global::REAL, "3"),CoefXpr::OP_DIV),
        CoefXpr::OP_SUB ));
    }

    // change subType to be able to use mech-integrators for the 2D case
    std::string subType;
    if( subType_ == "plane" ) {
      subType = "planeStrain";
    } else {
      subType = subType_;
    }

    // check which integrators we have to add
    bool setP1 = false;
    bool setP2 = false;
    bool setP3 = false;
    if( integrator == "all" ) {
      setP1 = true;
      setP2 = true;
      setP3 = true;
    } else if ( integrator == "P1" ) {
      setP1 = true;
    } else if ( integrator == "P2" ) {
      setP2 = true;
    } else if ( integrator == "P3" ) {
      setP3 = true;
    } else {
      EXCEPTION("Implementation error: Unknown integrator type specified!");
    }

    if( setP1 == true ) {
      presFct->AddEntityList( actSDList );
    }
    velFct->AddEntityList( actSDList );

    if( setP1 == true ) {
      // --------------------------------------------------------------------
      //  K_VP Surface Integrator - part 1 - pressure part
      //  (upper off-diagonal integrators - partially integrated, surface)
      // --------------------------------------------------------------------
      BiLinearForm * stiffIntVP1Surf = NULL;
      if( dim_ == 2 ) {
          stiffIntVP1Surf = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                new IdentityOperatorNormal<FeH1,2>(),
                                constOne, balanceOfMomentumSign_, volRegion, updatedGeo_);
      } else {
          stiffIntVP1Surf = new SurfaceABInt<>(new IdentityOperator<FeH1,3,3>(),
                                new IdentityOperatorNormal<FeH1,3>(),
                                constOne, balanceOfMomentumSign_, volRegion, updatedGeo_);
      }
      stiffIntVP1Surf->SetName("LinFlowStiffIntVP1Surf");
      BiLinFormContext *stiffContVP1 = NULL;
      stiffContVP1 = new BiLinFormContext(stiffIntVP1Surf, STIFFNESS );
      stiffContVP1->SetEntities( actSDList, actSDList );
      stiffContVP1->SetFeFunctions( velFct, presFct );
      // In case the LinFLowPDE is coupled to the HeatPDE in a symmetric form the counterpart needs to be set, to also
      // create symmetric counterpart
      stiffContVP1->SetCounterPart(isHeatPDECoupled_);
      assemble_->AddBiLinearForm( stiffContVP1 );
    }


    if( setP2 == true ) {
      // --------------------------------------------------------------------
      //  K_VP Surface Integrator - part 2 - strain part
      //  (upper off-diagonal integrators - partially integrated, surface)
      // --------------------------------------------------------------------
      BiLinearForm * stiffIntVP2Surf = NULL;
      BaseBOperator * sNSOp1 = NULL;
      if( dim_ == 2 ) {
        if( isComplex_ ) {
          sNSOp1 = new SurfaceNormalStressOperator<FeH1,2,2,Complex>(subType, false);
          sNSOp1->SetCoefFunction(coefBB);
          stiffIntVP2Surf = new SurfaceABInt<Complex,Complex>(new SurfaceIdentityOperator<FeH1,2,2>(),
                                sNSOp1,
                                constOne, -balanceOfMomentumSign_, volRegion, updatedGeo_);
        } else {
          sNSOp1 = new SurfaceNormalStressOperator<FeH1,2,2,Double>(subType, false);
          sNSOp1->SetCoefFunction(coefBB);
          stiffIntVP2Surf = new SurfaceABInt<Double,Double>(new SurfaceIdentityOperator<FeH1,2,2>(),
                                sNSOp1,
                                constOne, -balanceOfMomentumSign_, volRegion, updatedGeo_);
        }
      } else {
        if( isComplex_ ) {
          sNSOp1 = new SurfaceNormalStressOperator<FeH1,3,3,Complex>(subType, false);
          sNSOp1->SetCoefFunction(coefBB);
          stiffIntVP2Surf = new SurfaceABInt<Complex,Complex>(new SurfaceIdentityOperator<FeH1,2,2>(),
                                sNSOp1,
                                constOne, -balanceOfMomentumSign_, volRegion, updatedGeo_);
        } else {
          sNSOp1 = new SurfaceNormalStressOperator<FeH1,3,3,Double>(subType, false);
          sNSOp1->SetCoefFunction(coefBB);
          stiffIntVP2Surf = new SurfaceABInt<Double,Double>(new SurfaceIdentityOperator<FeH1,2,2>(),
                                sNSOp1,
                                constOne, -balanceOfMomentumSign_, volRegion, updatedGeo_);
        }
      }
      stiffIntVP2Surf->SetName("LinFlowStiffIntVP2Surf");
      // we set the volume evaluation for the BiLinForm since the context will automatically use the value of the BiLinForm
      // we have to do this before the initialization of the context!
      stiffIntVP2Surf->SetUseVolEqnA( true );
      stiffIntVP2Surf->SetUseVolEqnB( true ); 
      BiLinFormContext *stiffContVP2 = NULL;
      stiffContVP2 = new BiLinFormContext(stiffIntVP2Surf, STIFFNESS );
      stiffContVP2->SetEntities( actSDList, actSDList );
      stiffContVP2->SetFeFunctions( velFct, velFct );
      // In case the LinFLowPDE is coupled to the HeatPDE in a symmetric form the counterpart needs to be set, to also
      // create symmetric counterpart
      stiffContVP2->SetCounterPart(isHeatPDECoupled_);
      assemble_->AddBiLinearForm( stiffContVP2 );
    }


    if( setP3 == true ) {
      // --------------------------------------------------------------------
      //  K_VP Surface Integrator - part 3 - div part
      //  (upper off-diagonal integrators - partially integrated, surface)
      // --------------------------------------------------------------------
      if ( isCompressible_ ) {
        BiLinearForm * stiffIntVP3Surf = NULL;
        BaseBOperator * sNDOp1 = NULL;
        if( dim_ == 2 ) {
          if( isComplex_ ) {
            sNDOp1 = new SurfaceNormalDivOperator<FeH1,2,Complex>();
            sNDOp1->SetCoefFunction(coefDivDiv);
            stiffIntVP3Surf = new SurfaceABInt<Complex,Complex>(new SurfaceIdentityOperator<FeH1,2,2>(),
                                  sNDOp1,
                                  constOne, -balanceOfMomentumSign_, volRegion, updatedGeo_);
          } else {
            sNDOp1 = new SurfaceNormalDivOperator<FeH1,2,Double>();
            sNDOp1->SetCoefFunction(coefDivDiv);
            stiffIntVP3Surf = new SurfaceABInt<Double,Double>(new SurfaceIdentityOperator<FeH1,2,2>(),
                                  sNDOp1,
                                  constOne, -balanceOfMomentumSign_, volRegion, updatedGeo_);
          }
        } else {
          if( isComplex_ ) {
            sNDOp1 = new SurfaceNormalDivOperator<FeH1,3,Complex>();
            sNDOp1->SetCoefFunction(coefDivDiv);
            stiffIntVP3Surf = new SurfaceABInt<Complex,Complex>(new SurfaceIdentityOperator<FeH1,2,2>(),
                                  sNDOp1,
                                  constOne, -balanceOfMomentumSign_, volRegion, updatedGeo_);
          } else {
            sNDOp1 = new SurfaceNormalDivOperator<FeH1,3,Double>();
            sNDOp1->SetCoefFunction(coefDivDiv);
            stiffIntVP3Surf = new SurfaceABInt<Double,Double>(new SurfaceIdentityOperator<FeH1,2,2>(),
                                  sNDOp1,
                                  constOne, -balanceOfMomentumSign_, volRegion, updatedGeo_);
          }
        }
        stiffIntVP3Surf->SetName("LinFlowStiffIntVP3Surf");
        // we set the volume evaluation for the BiLinForm since the context will automatically use the value of the BiLinForm
        // we have to do this before the initialization of the context!
        stiffIntVP3Surf->SetUseVolEqnA( true ); 
        stiffIntVP3Surf->SetUseVolEqnB( true );
        BiLinFormContext *stiffContVP3 = NULL;
        stiffContVP3 = new BiLinFormContext(stiffIntVP3Surf, STIFFNESS );
        stiffContVP3->SetEntities( actSDList, actSDList );
        stiffContVP3->SetFeFunctions( velFct, velFct );
        // In case the LinFLowPDE is coupled to the HeatPDE in a symmetric form the counterpart needs to be set, to also
        // create symmetric counterpart
        stiffContVP3->SetCounterPart(isHeatPDECoupled_);
        assemble_->AddBiLinearForm( stiffContVP3 );
      }
    }
  }
}
