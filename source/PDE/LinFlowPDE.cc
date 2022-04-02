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
    isHeatCoupled_  = false;

    //! Always use total Lagrangian formulation 
    updatedGeo_        = true;
 
    //! check, if a compressible formulation is specified
    isCompressible_ = false;
    std::string formulation = myParam_->Get("formulation")->As<std::string>();
    if ( formulation == "compressible")
        isCompressible_ = true;

    // Check the subtype of the problem
    paramNode->GetValue("subType", subType_);

    if(ptGrid_->GetDim() == 3)
      subType_ = "3d";

    if(subType_ == "axi")
      EXCEPTION("SubType 'axi' not implemented for LinFlowPDE");
      
    //get oder of FE basis functions
    presPolyId_ = myParam_->Get("presPolyId")->As<std::string>();
    velPolyId_  = myParam_->Get("velPolyId")->As<std::string>();

    //get order of integration
    presIntegId_ = myParam_->Get("presIntegId")->As<std::string>();
    velIntegId_  = myParam_->Get("velIntegId")->As<std::string>();

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
  }
  
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

    //  Loop over all regions
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {

      // Set current region and material
      actRegion = it->first;

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

      // ====================================================================
      // stiffness integrators: conservation of mass
      //
      //  K_PV Integrator (lower off-diagonal integrator):
      //  \int_{\Omega_f} phi div(v') d\Omega
      // ===================================================================
      PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
      BiLinearForm * stiffIntPV = NULL;
      if( dim_ == 2 ) {
        stiffIntPV = new ABInt<>( new MultiIdOp<FeH1,2>(), 
                                  new DivOperator<FeH1,2>(), constOne, 1.0, updatedGeo_ );
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
          if ( isHeatCoupled_ ) {
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
                                density, 1.0, updatedGeo_ );
      } else {
        dampIntvv = new BBInt<>(new IdentityOperator<FeH1,3,3>(),
                                density, 1.0, updatedGeo_ );
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
        stiffIntVP = new ABInt<>( new DivOperator<FeH1,2>(),
                                  new MultiIdOp<FeH1,2>(), coeffKVP, -1.0, updatedGeo_ );
      } else {
        stiffIntVP = new ABInt<>( new DivOperator<FeH1,3>(),
                                  new MultiIdOp<FeH1,3>(), coeffKVP, -1.0, updatedGeo_ );
      }
      stiffIntVP->SetName("LinFlowStiffIntVP");
      BiLinFormContext *stiffContVP = NULL;
      stiffContVP = new BiLinFormContext(stiffIntVP, STIFFNESS );

      stiffContVP->SetEntities( actSDList, actSDList );
      stiffContVP->SetFeFunctions( velFct, presFct );
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
      StdVector<PtrCoefFct> tensorComponents(dim_ == 2 ? 9 : 36);
      tensorComponents.Init(coefZero);
      if( dim_ == 2 ) {
        tensorComponents[0] = shearViscosityDouble;
        tensorComponents[4] = shearViscosityDouble;
        tensorComponents[8] = shearViscosity;
        PtrCoefFct coefBB = CoefFunction::Generate(mp_,Global::REAL,3,3,tensorComponents);
        stiffIntLaplace = new BDBInt<>( new StrainOperator2D<FeH1>(), coefBB, 1.0, updatedGeo_ );
      } else {
        tensorComponents[0] = shearViscosityDouble;
        tensorComponents[7] = shearViscosityDouble;
        tensorComponents[14] = shearViscosityDouble;
        tensorComponents[21] = shearViscosity;
        tensorComponents[28] = shearViscosity;
        tensorComponents[35] = shearViscosity;
        PtrCoefFct coefBB = CoefFunction::Generate(mp_,Global::REAL,6,6,tensorComponents);
        stiffIntLaplace = new BDBInt<>( new StrainOperator3D<FeH1>(), coefBB, 1.0, updatedGeo_ );
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
          stiffIntDivDiv = new BBInt<>(new ScalarDivergenceOperator<FeH1,2,Double>(), coefDivDiv, 1.0, updatedGeo_);
        } else {
          stiffIntDivDiv = new BBInt<>(new ScalarDivergenceOperator<FeH1,3,Double>(), coefDivDiv, 1.0, updatedGeo_);
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
      if(flowId != ""){

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
            convectiveVv = new ABInt<Complex>( new IdentityOperator<FeH1,2,2>(), new ConvectiveOperator<FeH1,2,2,Complex>(), density, 1.0, coefUpdateGeo );
          }
          else
          {
            convectiveVv = new ABInt<Double>( new IdentityOperator<FeH1,2,2>(), new ConvectiveOperator<FeH1,2,2>(), density, 1.0, coefUpdateGeo );
          }
        } else {
          if(isComplex_) 
          {
            convectiveVv = new ABInt<Complex>( new IdentityOperator<FeH1,3,3>(), new ConvectiveOperator<FeH1,3,3,Complex>(), density, 1.0, coefUpdateGeo );
          }
          else
          {
            convectiveVv = new ABInt<Double>( new IdentityOperator<FeH1,3,3>(), new ConvectiveOperator<FeH1,3,3>(), density, 1.0, coefUpdateGeo );
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
            convectivevV = new BDBInt<Complex,Complex>( bOpId, rho0GradV0, 1.0 );
          } else {
            convectivevV = new BDBInt<Double,Double>( bOpId, rho0GradV0, 1.0 );
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
              convectiveVp = new ADBInt<Complex>( new MultiIdOp<FeH1,1,2,Complex>(), new MultiIdOp<FeH1,2,1,Complex>(), coefPv, 1.0, coefUpdateGeo );
            } else {
              convectiveVp = new ADBInt<Complex>( new MultiIdOp<FeH1,1,3,Complex>(), new MultiIdOp<FeH1,3,1,Complex>(), coefPv, 1.0, coefUpdateGeo );
            }
          } else {
            if(dim_==2){
              convectiveVp = new ADBInt<Double>( new MultiIdOp<FeH1,1,3,Double>(), new MultiIdOp<FeH1,3,1,Double>(), coefPv, 1.0, coefUpdateGeo );
            } else {
              convectiveVp = new ADBInt<Double>( new MultiIdOp<FeH1,1,3,Double>(), new MultiIdOp<FeH1,3,1,Double>(), coefPv, 1.0, coefUpdateGeo );
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

            if (isHeatCoupled_) {
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
            if(isComplex_)
            {
              stiffIntMovingGridVV = new ABInt<Complex>( new IdentityOperator<FeH1,2,2>(),
                                                  new ConvectiveOperator<FeH1,2,2,Complex>(),
                                                  density, -1.0, updatedGeo_ );
            }
            else
            {
              stiffIntMovingGridVV = new ABInt<Double>( new IdentityOperator<FeH1,2,2>(),
                                                new ConvectiveOperator<FeH1,2,2>(),
                                                density, -1.0, updatedGeo_ );
            }
          } else {
            if(isComplex_)
            {
              stiffIntMovingGridVV = new ABInt<Complex>( new IdentityOperator<FeH1,3,3>(),
                                                  new ConvectiveOperator<FeH1,3,3,Complex>(),
                                                  density, -1.0, updatedGeo_ );
            }
            else
            {
              stiffIntMovingGridVV = new ABInt<Double>( new IdentityOperator<FeH1,3,3>(),
                                                new ConvectiveOperator<FeH1,3,3>(),
                                                density, -1.0, updatedGeo_ );
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
                                     oneFuncs, 1.0, flowRegions, updatedGeo_);
            } else {
                stiffIntVPSurf = new SurfaceABInt<>(new IdentityOperator<FeH1,3,3>(),
                                     new IdentityOperatorNormal<FeH1,3>(),
                                     oneFuncs, 1.0, flowRegions, updatedGeo_);
            }
            stiffIntVPSurf->SetName("LinFlowStiffIntVPSurf");
            BiLinFormContext *stiffContVP = NULL;
            stiffContVP = new BiLinFormContext(stiffIntVPSurf, STIFFNESS );

            stiffContVP->SetEntities( actSDList, actSDList );
            stiffContVP->SetFeFunctions( velFct, presFct );
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
              Complex(1.0), coef[i], coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,2,2>(),
              1.0, coef[i],coefUpdateGeo);
        }
      } else  {
        if(isComplex_) {
          lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,3,3>(),
              Complex(1.0), coef[i], coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,3,3>(),
              1.0, coef[i], coefUpdateGeo);
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
              Complex(tracFac), coef[i],
              volRegions, coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double,true> ( new IdentityOperatorNormalTrans<FeH1,2>(),
              tracFac, coef[i], volRegions,
              coefUpdateGeo);
        }
      } else  {
        if(isComplex_) {
          lin = new BUIntegrator<Complex, true> ( new IdentityOperatorNormalTrans<FeH1,3>(),
              Complex(tracFac), coef[i],
              volRegions, coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double, true> ( new IdentityOperatorNormalTrans<FeH1,3>(),
              tracFac, coef[i],
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
  }

  void LinFlowPDE::DefineSurfaceIntegrators(){
    // Get FESpace and FeFunction of fluid velocity
    shared_ptr<BaseFeFunction> velFct = feFunctions_[FLUIDMECH_VELOCITY];
    shared_ptr<FeSpace> velSpace = velFct->GetFeSpace();
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
            surfMassInt = new SurfaceBBInt<Complex,Complex>(new IdentityOperatorNormalTrans<FeH1,2,2,Complex>(), kCoef[i], Complex(1.0,0), volRegion, updatedGeo_ );
          } else {
            surfMassInt = new SurfaceBBInt<Complex,Complex>(new IdentityOperatorNormalTrans<FeH1,3,3,Complex>(), kCoef[i], Complex(1.0,0), volRegion, updatedGeo_ );
          }
        } else {
          if (dim_ == 2){
            surfMassInt = new SurfaceBBInt<>(new IdentityOperatorNormalTrans<FeH1,2,2>(), kCoef[i], 1.0, volRegion, updatedGeo_ );
          } else {
            surfMassInt = new SurfaceBBInt<>(new IdentityOperatorNormalTrans<FeH1,3,3>(), kCoef[i], 1.0, volRegion, updatedGeo_ );
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
            impedanceInt = new SurfaceBBInt<Complex,Complex>(new IdentityOperatorNormalTrans<FeH1,2,2,Complex>(), kCoef[i], Complex(1.0,0), volRegion, updatedGeo_ );
          } else {
            impedanceInt = new SurfaceBBInt<Complex,Complex>(new IdentityOperatorNormalTrans<FeH1,3,3,Complex>(), kCoef[i], Complex(1.0,0), volRegion, updatedGeo_ );
          }
        } else {
          if (dim_ == 2){
            impedanceInt = new SurfaceBBInt<>(new IdentityOperatorNormalTrans<FeH1,2,2>(), kCoef[i], 1.0, volRegion, updatedGeo_ );
          } else {
            impedanceInt = new SurfaceBBInt<>(new IdentityOperatorNormalTrans<FeH1,3,3>(), kCoef[i], 1.0, volRegion, updatedGeo_ );
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
      // z=c+di, and we want to find sqrt(z)=a+bi
      // a = sqrt((c+norm(z))/2), norm(z) = sqrt(c^2+d^2)
      // b = d/abs(d)*sqrt((-c+norm(z))/2)
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
        //finding c and d
        PtrCoefFct complexWaveNumRe1 = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_, complexWaveNum1, CoefXpr::OP_RE));//real
        PtrCoefFct complexWaveNumIm1 = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_, complexWaveNum1, CoefXpr::OP_IM));//imag
        // norm(z)
        PtrCoefFct complexWaveNumAbs = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_, complexWaveNum1, CoefXpr::OP_NORM));
        //finding a = sqrt((c+norm(z))/2), norm(z) = sqrt(c^2+d^2)
        PtrCoefFct complexWaveNumRe = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_,
            CoefFunction::Generate(mp_, Global::REAL,CoefXprBinOp(mp_,
                CoefFunction::Generate(mp_,Global::REAL,CoefXprBinOp(mp_,complexWaveNumAbs,complexWaveNumRe1, CoefXpr::OP_ADD)),
                CoefFunction::Generate( mp_, Global::REAL, "2"),
                CoefXpr::OP_DIV)),
                CoefXpr::OP_SQRT));
        //d/abs(d)
        PtrCoefFct dabsd = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_,complexWaveNumIm1,
            CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_, complexWaveNumIm1, CoefXpr::OP_NORM)), CoefXpr::OP_DIV));
        //finding b = d/abs(d)*sqrt((-c+norm(z))/2)
        PtrCoefFct complexWaveNumIm = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_,dabsd,
            CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_,
                CoefFunction::Generate(mp_, Global::REAL,CoefXprBinOp(mp_,
                    CoefFunction::Generate(mp_,Global::REAL,CoefXprBinOp(mp_, complexWaveNumAbs, complexWaveNumRe1, CoefXpr::OP_SUB)),
                    CoefFunction::Generate( mp_, Global::REAL, "2"),
                    CoefXpr::OP_DIV)),
                    CoefXpr::OP_SQRT)),
                    CoefXpr::OP_MULT));
        //create a+bi
        PtrCoefFct complexWaveNum = CoefFunction::Generate(mp_, Global::COMPLEX ,CoefXprBinOp(mp_,complexWaveNumRe,
            CoefFunction::Generate( mp_, Global::COMPLEX,CoefXprBinOp(mp_,
                complexWaveNumIm, CoefFunction::Generate( mp_, Global::COMPLEX, "0", "1"),
                CoefXpr::OP_MULT)),
                CoefXpr::OP_ADD));
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
            abcInt = new SurfaceBBInt<Complex,Complex>(new IdentityOperatorNormalTrans<FeH1,2,2,Complex>(), impCoeff, Complex(1.0,0), volRegion, updatedGeo_ );
          } else {
            abcInt = new SurfaceBBInt<Complex,Complex>(new IdentityOperatorNormalTrans<FeH1,3,3,Complex>(), impCoeff, Complex(1.0,0), volRegion, updatedGeo_ );
          }
        } else {
          if (dim_ == 2){
            abcInt = new SurfaceBBInt<>(new IdentityOperatorNormalTrans<FeH1,2,2>(), impCoeff, 1.0, volRegion, updatedGeo_ );
          } else {
            abcInt = new SurfaceBBInt<>(new IdentityOperatorNormalTrans<FeH1,3,3>(), impCoeff, 1.0, volRegion, updatedGeo_ );
          }
        }
        abcInt->SetName("abcIntegrator");
        BiLinFormContext *abcContext = new BiLinFormContext(abcInt, STIFFNESS );
        abcContext->SetEntities( actSDList, actSDList);
        abcContext->SetFeFunctions( feFunctions_[FLUIDMECH_VELOCITY], feFunctions_[FLUIDMECH_VELOCITY]);
        feFunctions_[FLUIDMECH_VELOCITY]->AddEntityList( actSDList );
        assemble_->AddBiLinearForm( abcContext );
      }
    }
  }

  void LinFlowPDE::DefineSolveStep() {
    solveStep_ = new StdSolveStep(*this);
  }

  void LinFlowPDE::InitTimeStepping() {

    GLMScheme * schemeV = new Trapezoidal(0.5);
    GLMScheme * schemeP = new Trapezoidal(0.5);

    shared_ptr<BaseTimeScheme> mySchemeV(new TimeSchemeGLM(schemeV, 0) );
    shared_ptr<BaseTimeScheme> mySchemeP(new TimeSchemeGLM(schemeP, 0) );
    feFunctions_[FLUIDMECH_VELOCITY]->SetTimeScheme(mySchemeV);
    feFunctions_[FLUIDMECH_PRESSURE]->SetTimeScheme(mySchemeP);
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

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    idbcSolNameMap_[FLUIDMECH_PRESSURE] = "pressure";
    idbcSolNameMap_[FLUIDMECH_VELOCITY] = "velocity";
  
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
      availResults_.insert( stressVisc );
      PtrCoefFct stressViscCoef;
      // we define the stressViscCoef in the next section since we have to differ between compressible / incompressible formulation


      if ( isCompressible_ ) {
        // === FLUID-MECHANIC STRESS (COMPRESSIBLE) ===
        // (lambda - 2/3 mu) div(v) I
        shared_ptr<ResultInfo> stressComp(new ResultInfo);
        stressComp->resultType = FLUIDMECH_COMP_STRESS;
        stressComp->dofNames = stressComponents;
        stressComp->unit = MapSolTypeToUnit(FLUIDMECH_COMP_STRESS);
        stressComp->entryType = ResultInfo::TENSOR;
        stressComp->definedOn = ResultInfo::ELEMENT;
        availResults_.insert( stressComp );
        shared_ptr<CoefFunctionFormBased> stressCompFunc;
        if( isComplex_ ) {
          stressCompFunc.reset(new CoefFunctionFlux<Complex>(velFeFct, stressComp));
        } else {
          stressCompFunc.reset(new CoefFunctionFlux<Double>(velFeFct, stressComp));
        }
        // we add a specific integrator name to ensure that only this integrator gets passed to the coefFunction during the assignment in the SinglePDE during FinalizePostProcResults
        stressCompFunc->SetIntegratorName("LinFlowStiffIntBulkViscous"); // coefFunctionFlux acts on strain part
        StdVector<PtrCoefFct> stressCompDiagValues = StdVector<PtrCoefFct>(dim_);
        for(UInt i = 0; i < dim_; i++){
          stressCompDiagValues[i] = stressCompFunc;
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
      availResults_.insert( presTens );
      PtrCoefFct presFnc = this->GetCoefFct( FLUIDMECH_PRESSURE );
      StdVector<PtrCoefFct> presTensDiagValues = StdVector<PtrCoefFct>(dim_);
      for(UInt i = 0; i < dim_; i++){
        presTensDiagValues[i] = presFnc;
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


      // === FLUID-MECHANIC NORMAL SURFACE STRESS ===
      shared_ptr<ResultInfo> surfaceNormalStressInfo;
      shared_ptr<CoefFunctionSurf> surfaceNormalStressFct;
      surfaceNormalStressInfo.reset(new ResultInfo);
      surfaceNormalStressInfo->resultType = FLUIDMECH_NORMAL_SURFACE_STRESS;
      surfaceNormalStressInfo->dofNames = dispDofNames;
      surfaceNormalStressInfo->unit = MapSolTypeToUnit(FLUIDMECH_NORMAL_SURFACE_STRESS);
      surfaceNormalStressInfo->entryType = ResultInfo::VECTOR;
      surfaceNormalStressInfo->definedOn = ResultInfo::SURF_ELEM;

      surfaceNormalStressFct.reset(new CoefFunctionSurf(true, 1.0, surfaceNormalStressInfo));
      DefineFieldResult(surfaceNormalStressFct, surfaceNormalStressInfo);
      surfCoefFcts_[surfaceNormalStressFct] = stressTotalCoef;

      // === FLUID-MECHANIC REACTION FORCE (= integral of surface traction, i.e. normal stress from above, over the surface region ) ===
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
          reactionForceFct.reset(new ResultFunctorIntegrate<Complex>(surfaceNormalStressFct, feFct, reactionForceInfo));
      else
          reactionForceFct.reset(new ResultFunctorIntegrate<Double>(surfaceNormalStressFct, feFct, reactionForceInfo));
      resultFunctors_[FLUIDMECH_FORCE] = reactionForceFct;
      availResults_.insert(reactionForceInfo);

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

}
