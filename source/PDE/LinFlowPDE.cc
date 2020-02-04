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
                                  new DivOperator<FeH1,2>(), constOne, 1.0 );
      } else {
        stiffIntPV = new ABInt<>( new MultiIdOp<FeH1,3>(),
                                  new DivOperator<FeH1,3>(), constOne, 1.0 );
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
                                density, 1.0 );
      } else {
        dampIntvv = new BBInt<>(new IdentityOperator<FeH1,3,3>(),
                                density, 1.0 );
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
                                  new MultiIdOp<FeH1,2>(), coeffKVP, -1.0 );
      } else {
        stiffIntVP = new ABInt<>( new DivOperator<FeH1,3>(),
                                  new MultiIdOp<FeH1,3>(), coeffKVP, -1.0 );
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
      BiLinearForm * stiffIntLaplace = NULL;
      StdVector<PtrCoefFct> tensorComponents(dim_ == 2 ? 9 : 36);
      tensorComponents.Init(coefZero);
      if( dim_ == 2 ) {
        tensorComponents[0] = shearViscosityDouble;
        tensorComponents[4] = shearViscosityDouble;
        tensorComponents[8] = shearViscosity;
        PtrCoefFct coefBB = CoefFunction::Generate(mp_,Global::REAL,3,3,tensorComponents);
        stiffIntLaplace = new BDBInt<>( new StrainOperator2D<FeH1>(), coefBB, 1.0 );
      } else {
        tensorComponents[0] = shearViscosityDouble;
        tensorComponents[7] = shearViscosityDouble;
        tensorComponents[14] = shearViscosityDouble;
        tensorComponents[21] = shearViscosity;
        tensorComponents[28] = shearViscosity;
        tensorComponents[35] = shearViscosity;
        PtrCoefFct coefBB = CoefFunction::Generate(mp_,Global::REAL,6,6,tensorComponents);
        stiffIntLaplace = new BDBInt<>( new StrainOperator3D<FeH1>(), coefBB, 1.0 );
      }
      stiffIntLaplace->SetName("LinFlowStiffIntViscous");
      BiLinFormContext *stiffContLaplace;
      stiffContLaplace = new BiLinFormContext(stiffIntLaplace, STIFFNESS );

      stiffContLaplace->SetEntities( actSDList, actSDList );
      stiffContLaplace->SetFeFunctions( velFct, velFct );
      assemble_->AddBiLinearForm( stiffContLaplace );

      if ( isCompressible_ ) { // we need to subtract 2*mu/3 Div(v)I and add the bulk part lambda*Div(v)I
        PtrCoefFct bulkViscosity = materials_[actRegion]->GetScalCoefFnc(FLUID_BULK_VISCOSITY, Global::REAL);
        PtrCoefFct coefDivDiv = CoefFunction::Generate( mp_,  Global::REAL,
            CoefXprBinOp(mp_,
                bulkViscosity,
                CoefXprBinOp(mp_,shearViscosityDouble,CoefFunction::Generate( mp_, Global::REAL, "3"),CoefXpr::OP_DIV),
            CoefXpr::OP_SUB ));

        BiLinearForm * stiffIntDivDiv = NULL;
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
                                     oneFuncs, 1.0, flowRegions);
            } else {
                stiffIntVPSurf = new SurfaceABInt<>(new IdentityOperator<FeH1,3,3>(),
                                     new IdentityOperatorNormal<FeH1,3>(),
                                     oneFuncs, 1.0, flowRegions);
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
      DefineFieldResult( sigmaFunc, stress );
      
      // === FLUID-MECHANIC STRAINRATE ===
      shared_ptr<ResultInfo> strain(new ResultInfo);
      strain->resultType = FLUIDMECH_STRAINRATE;
      strain->dofNames = stressComponents;
      strain->unit =  MapSolTypeToUnit(FLUIDMECH_STRAINRATE);;
      strain->entryType = ResultInfo::TENSOR;
      strain->definedOn = ResultInfo::ELEMENT;
      availResults_.insert( strain );
      shared_ptr<CoefFunctionFormBased> strainFunc;
      if( isComplex_ ) {
        strainFunc.reset(new CoefFunctionBOp<Complex>(velFeFct, strain));
      } else {
        strainFunc.reset(new CoefFunctionBOp<Double>(velFeFct, strain));
      }
      DefineFieldResult( strainFunc, strain );

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

    // ============================
    // Initialize result functors:
    // ============================

    //  Loop over all regions
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {      
      RegionIdType region = it->first;
      BaseMaterial* actSDMat = it->second;

      // 2) pass integrators to functors
      // eFunc->AddIntegrator(stiffIntVP, region);
      BaseBDBInt *integ = GetStiffIntegrator( actSDMat, region, isComplex_ );
      integ->SetName("ViscousStrainRateInt");
      integ->SetFeSpace( feFunctions_[FLUIDMECH_VELOCITY]->GetFeSpace() );
      sigmaFunc->AddIntegrator(integ,region);
      strainFunc->AddIntegrator(integ,region);
    }
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

  BaseBDBInt* LinFlowPDE::GetStiffIntegrator( BaseMaterial* actSDMat,
                                              RegionIdType regionId,
                                              bool isComplex ) {

    // Get region name
    std::string regionName = ptGrid_->GetRegion().ToString( regionId );

    // ------------------------
    //  Obtain linear material
    // ------------------------
    shared_ptr<CoefFunction > curCoef;

    BaseMaterial* regionMat = materials_[regionId];
    SubTensorType subTensorType = NO_TENSOR;
    Global::ComplexPart complexPart = Global::REAL;

    if( isComplex ) {
      complexPart = Global::COMPLEX;
    }
    
    if( subType_ == "axi" ) {
      subTensorType = AXI;
    } else if( subType_ == "plane" ) {
      subTensorType = PLANE;
    } else if( subType_ == "3d") {
      subTensorType = FULL;
    } else {
      EXCEPTION( "Subtype '" << subType_ << "' unknown for fluid-mechanic physic" );
    }

    curCoef = regionMat->GetTensorCoefFnc( FLUID_DYNAMIC_VISCOSITY, subTensorType, complexPart );
    //curCoef = regionMat->GetScalCoefFnc(DYNAMIC_VISCOSITY, Global::REAL);;
    // ----------------------------------------
    //  Determine correct stiffness integrator 
    // ----------------------------------------
    BaseBDBInt * integ = NULL;
    if( isComplex ) {
      if( subType_ == "axi" ) {
        integ = new BDBInt<Complex>(new StrainOperatorAxi<FeH1,Complex>(), curCoef, 1.0);
      } else if( subType_ == "plane" ) {
        integ = new BDBInt<Complex>(new StrainOperator2D<FeH1,Complex>(), curCoef, 1.0);
      } else if( subType_ == "3d") {
        integ = new BDBInt<Complex>(new StrainOperator3D<FeH1,Complex>(), curCoef, 1.0);
      } else {
        EXCEPTION( "Subtype '" << subType_ << "' unknown for fluid-mechanic physic" );
      }
    }
    else {
      if( subType_ == "axi" ) {
        integ = new BDBInt<>(new StrainOperatorAxi<FeH1>(), curCoef, 1.0);
      } else if( subType_ == "plane" ) {
        integ = new BDBInt<>(new StrainOperator2D<FeH1>(), curCoef, 1.0);
      } else if( subType_ == "3d") {
        integ = new BDBInt<>(new StrainOperator3D<FeH1>(), curCoef, 1.0);
      } else {
        EXCEPTION( "Subtype '" << subType_ << "' unknown for fluid-mechanic physic" );
      }
    }

    return integ;
  }


}
