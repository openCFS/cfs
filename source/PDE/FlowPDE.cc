// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <set>

#include "FlowPDE.hh"

#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionStabParams.hh"

#include "Utils/StdVector.hh"
#include "Driver/Assemble.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/FeFunctions.hh"
#include "Domain/ElemMapping/ElemShapeMap.hh"

// new integrator concept
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/KXInt.hh"

#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/DivOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorProjected.hh"
#include "Forms/Operators/MultiIdOp.hh"
#include "Forms/Operators/LaplOp.hh"
#include "Forms/Operators/ConvectiveOperator.hh"
#include "Domain/CoefFunction/CoefFunctionMeanFlowConvection.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"
#include "Forms/Operators/StrainOperator.hh"

#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"

#include "Domain/CoefFunction/CoefXpr.hh"
#include "OLAS/algsys/SolStrategy.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"

namespace CoupledField {

  DEFINE_LOG(fluidmechpde, "pde.fluidmech")

  // ***************
  //   Constructor
  // ***************
  FlowPDE::FlowPDE( Grid* grid, PtrParamNode paramNode,
                                      PtrParamNode infoNode,
                                      shared_ptr<SimState> simState, 
                                      Domain* domain )
  :SinglePDE( grid, paramNode, infoNode,simState, domain )
  {
    // ======================================================================
    // set solution information
    // ======================================================================

    pdename_          = "fluidMech";
    pdematerialclass_ = FLOW;
 
    nonLinMaterial_ = false;

    //! Always use total Lagrangian formulation 
    updatedGeo_        = true;
 
    // Check the subtype of the problem
    paramNode->GetValue("subType", subType_);

    if(ptGrid_->GetDim() == 3)
      subType_ = "3d";
      
    // Obtain regions the pde is defined on
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

    //check for stabilization!!
    stabilizedBochev_ = false;
    stabilizedBochev_ = myParam_->Get("stabilizationBochev")->As<bool>();
    if ( stabilizedBochev_ )
          std::cerr << "\n DO STABILIZED FORMULATION OF TYPE BOCHEV!!\n" << std::endl;

    addBochevPressStab_ = false;
    addBochevPressStab_= myParam_->Get("addBochevPressStabilization")->As<bool>();

    //check for FE spcae formulation
    myParam_->GetValue("feSpaceFormulation",FEspaceFormulation_,ParamNode::EX);
  }
  
  void FlowPDE::InitNonLin() {
    nonLinMethod_ = FIXEDPOINT;
    PtrParamNode nonLinNode = solStrat_->GetNonLinNode();
    if( nonLinNode ) {
      std::string methodString;
      nonLinNode->GetValue(  "method", methodString, ParamNode::PASS );
      nonLinMethod_ = NonLinMethodTypeEnum.Parse(methodString);
    }

    nonLin_ = true;
    //nonLinTotalFormulation_ = true;
  }

  void FlowPDE::DefineIntegrators() {

    bool updateGeo = false;

    RegionIdType actRegion;

    // Get FESpace and FeFunction of fluid velocity
    shared_ptr<BaseFeFunction> velFct = feFunctions_[FLUIDMECH_VELOCITY];
    shared_ptr<FeSpace> velSpace = velFct->GetFeSpace();

    // Get FESpace and FeFunction of fluid pressure
    shared_ptr<BaseFeFunction> presFct = feFunctions_[FLUIDMECH_PRESSURE];
    shared_ptr<FeSpace> presSpace = presFct->GetFeSpace();

    // Create coefficient functions for all fluid densities
    std::map< RegionIdType, PtrCoefFct > oneFuncs;
    std::set< RegionIdType > flowRegions;

    //  Loop over all regions
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion++){
      // Set current region and material
      actRegion = regions_[iRegion];

    // Not needed at the moment. Commented out due to gcc 4.6.
#if 0
      actMat = it->second;
#endif

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

      // --------------------------
      //  Set region approximation
      // --------------------------
      // We hardcode the Taylor-Hood spaces for the moment
//      velSpace->SetRegionApproximation(actRegion, "velPolyId", "velIntegId");
//      presSpace->SetRegionApproximation(actRegion, "presPolyId", "presIntegId");
      velSpace->SetRegionApproximation(actRegion, "velPolyId", "default");
      presSpace->SetRegionApproximation(actRegion, "presPolyId", "default");

      PtrCoefFct density = materials_[actRegion]->GetScalCoefFnc(
        DENSITY,
        Global::REAL
        );
      PtrCoefFct viscosity = materials_[actRegion]->GetScalCoefFnc(
        FLUID_DYNAMIC_VISCOSITY,
        Global::REAL
        );

      // Create set of flow regions and map of density functions for surface integrators.
      flowRegions.insert(actRegion);
      oneFuncs[actRegion] = 
        CoefFunction::Generate( mp_, Global::REAL,
                                lexical_cast<std::string>(1.0) );

      // ====================================================================
      // stiffness integrators
      // ====================================================================
      // --------------------------------------------------------------------
      //  VERSION 1: K_PV Integrator (lower off-diagonal integrator)
      //  \int_{\Omega_f} phi div(v') d\Omega = 0
      // --------------------------------------------------------------------
      BiLinearForm * stiffIntPV = NULL;
      if( dim_ == 2 ) {
        stiffIntPV = new ABInt<>( new MultiIdOp<FeH1,2>(), 
                                  new DivOperator<FeH1,2>(), density, 1.0 );
      } else {
        stiffIntPV = new ABInt<>( new MultiIdOp<FeH1,3>(),
                                  new DivOperator<FeH1,3>(), density, 1.0 );
      }
      stiffIntPV->SetName("FlowStiffIntPV");
      BiLinFormContext *stiffContPV = NULL;
      stiffContPV = new BiLinFormContext(stiffIntPV, STIFFNESS );

      stiffContPV->SetEntities( actSDList, actSDList );
      stiffContPV->SetFeFunctions( presFct, velFct);
      assemble_->AddBiLinearForm( stiffContPV );

      if ( FEspaceFormulation_ == "L2" ) {
        //      // --------------------------------------------------------------------
        //      //  VERSION 2: K_VP Integrator
        //      //  (upper off-diagonal integrators - partially integrated, volume)
        //      // --------------------------------------------------------------------
        PtrCoefFct coeffKVP = CoefFunction::Generate(mp_,Global::REAL, "1.0");
        BiLinearForm * stiffIntVP = NULL;
        if( dim_ == 2 ) {
          stiffIntVP = new ABInt<>( new DivOperator<FeH1,2>(),
              new MultiIdOp<FeH1,2>(), coeffKVP, -1.0 );
        } else {
          stiffIntVP = new ABInt<>( new DivOperator<FeH1,3>(),
              new MultiIdOp<FeH1,3>(), coeffKVP, -1.0 );
        }
        stiffIntVP->SetName("FlowStiffIntVP");
        BiLinFormContext *stiffContVP = NULL;
        stiffContVP = new BiLinFormContext(stiffIntVP, STIFFNESS );

        stiffContVP->SetEntities( actSDList, actSDList );
        stiffContVP->SetFeFunctions( velFct, presFct );
        assemble_->AddBiLinearForm( stiffContVP );

        // add penalty terms
        BiLinearForm *penaltyPresOpp = NULL;
        BiLinearForm *penaltyPres = NULL;
//        BiLinearForm *penaltyPressEx = NULL;
        std::set<RegionIdType> volRegion;
        volRegion.insert(actRegion);
        if( dim_ == 2 ) {
          penaltyPresOpp = new BBInt<Double>(new IdentityOperator<FeH1,2,1,Double>(),
              density, -0.5, updatedGeo_ );
          penaltyPres = new BBInt<Double>(new IdentityOperator<FeH1,2,1,Double>(),
              density, 0.5, updatedGeo_);
//          penaltyPressEx = new BBInt<Double>(new IdentityOperator<FeH1,2,1,Double>(),
//              density, 0.5, updatedGeo_);
        }
        else {
          penaltyPresOpp = new BBInt<Double>(new IdentityOperator<FeH1,3,1,Double>(),
               density, -0.5, updatedGeo_ );
           penaltyPres = new BBInt<Double>(new IdentityOperator<FeH1,3,1,Double>(),
               density, 0.5, updatedGeo_);
//           penaltyPressEx = new BBInt<Double>(new IdentityOperator<FeH1,3,1,Double>(),
//               density, 0.5, updatedGeo_);
        }

        penaltyPresOpp->SetName("penaltyPressOpposite");
        penaltyPres->SetName("penalityPress");
//        penaltyPressEx->SetName("penaltyPressExterior");

        //now we obtain the entity lists
        shared_ptr<EntityList> list,oppositList; //,extList;
        presFct->GetFeSpace()->GetInteriorSurfaceElems(actRegion,list,oppositList);
//        presFct->GetFeSpace()->GetExteriorSurfaceElems(actRegion,extList);

        BiLinFormContext * penaltyContextPressOpp =  new BiLinFormContext(penaltyPresOpp, STIFFNESS );
        BiLinFormContext * penaltyContextPress =  new BiLinFormContext(penaltyPres, STIFFNESS );
//        BiLinFormContext * penaltyContextPressExt =  new BiLinFormContext(penaltyPressEx, STIFFNESS );

        penaltyContextPressOpp->SetEntities(  oppositList,list );
        penaltyContextPress->SetEntities( list, list );
//        penaltyContextPressExt->SetEntities(extList,extList);

        penaltyContextPressOpp->SetFeFunctions( presFct, presFct);
        penaltyContextPress->SetFeFunctions( presFct,presFct);
//        penaltyContextPressExt->SetFeFunctions( presFct,presFct);
        assemble_->AddBiLinearForm( penaltyContextPressOpp );
        assemble_->AddBiLinearForm( penaltyContextPress );
//        assemble_->AddBiLinearForm( penaltyContextPressExt );
      }
      else {
        // --------------------------------------------------------------------
        //  VERSION 2: K_VP Integrator (upper off-diagonal integrator, not integration by parts)
        // --------------------------------------------------------------------
        PtrCoefFct coeffKVP
        = CoefFunction::Generate(mp_,Global::REAL, "1.0");
        BiLinearForm * stiffIntVP = NULL;
        if( dim_ == 2 ) {
          stiffIntVP = new ABInt<> (new IdentityOperator<FeH1,2,2>,
              new GradientOperator<FeH1,2>, coeffKVP, 1.0 );
        } else {
          stiffIntVP = new ABInt<> (new IdentityOperator<FeH1,3,3>,
              new GradientOperator<FeH1,3>, coeffKVP, 1.0);
        }
        stiffIntVP->SetName("FlowStiffIntVP");
        BiLinFormContext *stiffContVP = NULL;
        stiffContVP = new BiLinFormContext(stiffIntVP, STIFFNESS );

        stiffContVP->SetEntities( actSDList, actSDList );
        stiffContVP->SetFeFunctions( velFct, presFct );
        assemble_->AddBiLinearForm( stiffContVP );
      }

      // --------------------------------------------------------------------
      //  VERSION 2: K_Laplace Integrator
      // --------------------------------------------------------------------
      BiLinearForm * stiffIntLaplace = NULL;
      if( dim_ == 2 ) {
        stiffIntLaplace = new BBInt<>( new LaplOperator<FeH1,2>(),
                                       viscosity, 1.0 );
      } else {
        stiffIntLaplace = new BBInt<>( new LaplOperator<FeH1,3>(),
                                       viscosity, 1.0 );
      }
      stiffIntLaplace->SetName("FlowStiffIntViscous");
      BiLinFormContext *stiffContLaplace;
      stiffContLaplace = new BiLinFormContext(stiffIntLaplace, STIFFNESS );

      stiffContLaplace->SetEntities( actSDList, actSDList );
      stiffContLaplace->SetFeFunctions( velFct, velFct );
      assemble_->AddBiLinearForm( stiffContLaplace );

      // --------------------------------------------------------------------
      //  Convective Operator
      // --------------------------------------------------------------------
      BiLinearForm *convectiveVv = NULL;

      if( dim_ == 2 ) {
        convectiveVv = new ABInt<Double>( new IdentityOperator<FeH1,2,2>(),
            new ConvectiveOperator<FeH1,2,2>(),
            density, 1.0, updateGeo );
      } else {
        convectiveVv = new ABInt<Double>( new IdentityOperator<FeH1,3,3>(),
            new ConvectiveOperator<FeH1,3,3>(),
            density, 1.0, updateGeo );
      }

      PtrCoefFct velCoef = this->GetCoefFct(FLUIDMECH_VELOCITY);

      convectiveVv->SetBCoefFunctionOpB(velCoef);
      convectiveVv->SetSolDependent(true);

      convectiveVv->SetName("FlowStiffIntConvectiveVv");
      BiLinFormContext *convectiveContextVv = NULL;
      convectiveContextVv = new BiLinFormContext(convectiveVv, STIFFNESS );
      convectiveContextVv->SetEntities( actSDList, actSDList );
      convectiveContextVv->SetFeFunctions( feFunctions_[FLUIDMECH_VELOCITY],
                                           feFunctions_[FLUIDMECH_VELOCITY] );
      assemble_->AddBiLinearForm( convectiveContextVv );

      if( nonLinMethod_ == NEWTON ) {
        BaseBOperator* bOpGrad;
        BaseBOperator* bOpId;
        if( dim_ == 2 )
        bOpId = new IdentityOperator<FeH1,2,2>();
        else
          bOpId = new IdentityOperator<FeH1,3,3>();

        //now create the integrators
        BiLinearForm *convectivevV = NULL;
        PtrCoefFct coeffConvec;
        if( dim_ == 2 ) {
          bOpGrad = new GradientOperator<FeH1,2, 1, Double>();
          coeffConvec.reset(
              new CoefFunctionMeanFlowConvection<Double,2>( density, bOpGrad, feFunctions_[FLUIDMECH_VELOCITY]) );
        } else {
          bOpGrad = new GradientOperator<FeH1,3, 1, Double>();
          coeffConvec.reset(
              new CoefFunctionMeanFlowConvection<Double,3>( density, bOpGrad, feFunctions_[FLUIDMECH_VELOCITY]) );
        }
        convectivevV = new BDBInt<Double,Double>( bOpId, coeffConvec, 1.0 );

        convectivevV->SetName("FlowStiffIntConvectiveNewtonVv");
        //! mark the bi-linear form to be a Newton part
        convectivevV->SetNewtonBiLinearForm();

        BiLinFormContext *convectiveContextvV = NULL;
        convectiveContextvV = new BiLinFormContext(convectivevV, STIFFNESS );

        convectiveContextvV->SetEntities( actSDList, actSDList );
        convectiveContextvV->SetFeFunctions( feFunctions_[FLUIDMECH_VELOCITY],
                                             feFunctions_[FLUIDMECH_VELOCITY] );
        assemble_->AddBiLinearForm( convectiveContextvV );
      }

      if ( addBochevPressStab_ ) {
        //we add the pressure stabilization of Bochev
        PtrCoefFct coeffKPPstab = CoefFunction::Generate( mp_, Global::REAL, "1.0");

        BiLinearForm * stiffIntPPstab = NULL;

        if( dim_ == 2 ) {
          stiffIntPPstab = new BBInt<>( new IdentityOperatorProjected<FeH1,2>(),
              coeffKPPstab,1.0);
        } else {
          stiffIntPPstab = new BBInt<>( new IdentityOperatorProjected<FeH1,3>(),
              coeffKPPstab,1.0);
        }
        stiffIntPPstab->SetName("FlowStiffIntPPstab");
        //stiffIntPPstab->SetNewtonBiLinearForm();

        BiLinFormContext *stiffContPPstab = NULL;
        stiffContPPstab = new BiLinFormContext(stiffIntPPstab, STIFFNESS );

        stiffContPPstab->SetEntities( actSDList, actSDList );
        stiffContPPstab->SetFeFunctions( presFct, presFct );
        assemble_->AddBiLinearForm( stiffContPPstab );
      }

      // ====================================================================
      // damping integrators
      // ====================================================================
      BiLinearForm *dampIntvv = NULL;
      if( dim_ == 2 ) {
        dampIntvv = new BBInt<>(new IdentityOperator<FeH1,2,2>(),
                                density, 1.0 );
      } else {
        dampIntvv = new BBInt<>(new IdentityOperator<FeH1,3,3>(),
                                density, 1.0 );
      }
      dampIntvv->SetName("FlowDampInt");

      BiLinFormContext *dampContextvv = NULL;
      dampContextvv = new BiLinFormContext(dampIntvv, DAMPING );

      dampContextvv->SetEntities( actSDList, actSDList );
      dampContextvv->SetFeFunctions( feFunctions_[FLUIDMECH_VELOCITY],
                                     feFunctions_[FLUIDMECH_VELOCITY] );
      assemble_->AddBiLinearForm( dampContextvv );

      //======================================================================
      // DO STABILIZATION due to pressure
      //======================================================================
      if ( stabilizedBochev_) {
        //
        //   stabilization of pressure
        //
        PtrCoefFct coeffKPPstab = CoefFunction::Generate( mp_, Global::REAL, "1.0");

//        PtrCoefFct coeffKPPstab =
//            CoefFunction::Generate( mp_, Global::REAL,
//                CoefXprBinOp( mp_, density, viscosity, CoefXpr::OP_DIV ) );

        BiLinearForm * stiffIntPPstab = NULL;

        if( dim_ == 2 ) {
          stiffIntPPstab = new BBInt<>( new IdentityOperatorProjected<FeH1,2>(),
              coeffKPPstab,1.0);
        } else {
          stiffIntPPstab = new BBInt<>( new IdentityOperatorProjected<FeH1,3>(),
              coeffKPPstab,1.0);
        }
        stiffIntPPstab->SetName("FlowStiffIntPPstab");
        BiLinFormContext *stiffContPPstab = NULL;
        stiffContPPstab = new BiLinFormContext(stiffIntPPstab, STIFFNESS );
        
        stiffContPPstab->SetEntities( actSDList, actSDList );
        stiffContPPstab->SetFeFunctions( presFct, presFct );
        assemble_->AddBiLinearForm( stiffContPPstab );

        //
        //   stabilization of convective term
        //
        Double densityVal;
        LocPointMapped map;
        density->GetScalar(densityVal, map);
        BaseBOperator* bOpGrad;
        if( dim_ == 2 ) {
          bOpGrad = new GradientOperator<FeH1,2, 1, Double>();
        }
        else {
          bOpGrad = new GradientOperator<FeH1,3, 1, Double>();
        }

        PtrCoefFct velCoef = this->GetCoefFct(FLUIDMECH_VELOCITY);

        PtrCoefFct coeffConvecStab;
        coeffConvecStab.reset(
          new CoefFunctionStabParams( density, viscosity, velCoef,
                                      bOpGrad, presFct,
                                      CoefFunctionStabParams::SUPG, isComplex_ ) );
        BiLinearForm *convecStiffVVstab = NULL;
        if( dim_ == 2 ) {
          convecStiffVVstab = new BBInt<>(
              new ConvectiveOperator<FeH1,2,2>(),
                coeffConvecStab, densityVal, updateGeo);
        } else {
          convecStiffVVstab = new BBInt<>(
              new ConvectiveOperator<FeH1,3,3>(),
              coeffConvecStab, densityVal, updateGeo);
        }

        convecStiffVVstab->SetBCoefFunctionOpB(velCoef);
        convecStiffVVstab->SetSolDependent(true);
        convecStiffVVstab->SetName("FlowStiffIntVVConvectiveStab");

        BiLinFormContext *stiffContextVVStab = NULL;
        stiffContextVVStab = new BiLinFormContext( convecStiffVVstab,
                                                   STIFFNESS );
        stiffContextVVStab->SetEntities( actSDList, actSDList );
        stiffContextVVStab->SetFeFunctions( velFct, velFct);
        assemble_->AddBiLinearForm( stiffContextVVStab );
      }
    }
    
    StdVector<RegionIdType>::iterator surfIt, surfEnd;
    surfIt = presSurfaces_.Begin();
    surfEnd = presSurfaces_.End();

    for( ; surfIt != surfEnd; surfIt++ ) {
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( *surfIt );

      velFct->AddEntityList( actSDList );
      velSpace->SetRegionApproximation(*surfIt, "velPolyId", "default");
      presFct->AddEntityList( actSDList );
      presSpace->SetRegionApproximation(*surfIt, "presPolyId", "default");

      // --------------------------------------------------------------------
      //  VERSION 2: K_VP Integrator
      //  (upper off-diagonal integrators - partially integrated, surface)
      // --------------------------------------------------------------------
      BiLinearForm * stiffIntVPSurf = NULL;
      if( dim_ == 2 ) {
        stiffIntVPSurf = new SurfaceABInt<>(
          new IdentityOperator<FeH1,2,2>(),
          new IdentityOperatorNormal<FeH1,2>(),
          oneFuncs, 1.0, flowRegions
          );
      } else {
        stiffIntVPSurf = new SurfaceABInt<>(
          new IdentityOperator<FeH1,3,3>(),
          new IdentityOperatorNormal<FeH1,3>(),
          oneFuncs, 1.0, flowRegions
          );
      }
      stiffIntVPSurf->SetName("FlowStiffIntVPSurf");
      BiLinFormContext *stiffContVP = NULL; 
      stiffContVP = new BiLinFormContext(stiffIntVPSurf, STIFFNESS );

      stiffContVP->SetEntities( actSDList, actSDList );
      stiffContVP->SetFeFunctions( velFct, presFct );
      assemble_->AddBiLinearForm( stiffContVP );
    }

  }
  
  void FlowPDE::DefineSolveStep() {
    solveStep_ = new StdSolveStep(*this);
  }


  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================

  void FlowPDE::InitTimeStepping() {

    GLMScheme * schemeV = new Trapezoidal(0.5);
    GLMScheme * schemeP = new Trapezoidal(0.5);

//    Double gamma = 0.5;
//    GLMScheme * schemeV = new Trapezoidal(gamma);
//    GLMScheme * schemeP = new Trapezoidal(gamma);

    if ( nonLinTotalFormulation_ ) {
      shared_ptr<BaseTimeScheme> mySchemeV(new TimeSchemeGLM(schemeV, 0 ) );
      shared_ptr<BaseTimeScheme> mySchemeP(new TimeSchemeGLM(schemeP, 0 ) );
      feFunctions_[FLUIDMECH_VELOCITY]->SetTimeScheme(mySchemeV);
      feFunctions_[FLUIDMECH_PRESSURE]->SetTimeScheme(mySchemeP);
    }
    else {
      TimeSchemeGLM::NonLinType nlType = (nonLin_)? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
      shared_ptr<BaseTimeScheme> mySchemeV(new TimeSchemeGLM(schemeV, 0, nlType ) );
      shared_ptr<BaseTimeScheme> mySchemeP(new TimeSchemeGLM(schemeP, 0, nlType ) );
      feFunctions_[FLUIDMECH_VELOCITY]->SetTimeScheme(mySchemeV);
      feFunctions_[FLUIDMECH_PRESSURE]->SetTimeScheme(mySchemeP);
    }
  }

  void FlowPDE::ReadSpecialBCs( ) 
  {
  }
  
  void FlowPDE::DefinePrimaryResults() {
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
    pressure->unit = "Pa";

    pressure->definedOn = ResultInfo::MapSolTypeToDefinedOn(FLUIDMECH_PRESSURE);
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
    velocity->unit = "m/s";

    velocity->definedOn = ResultInfo::MapSolTypeToDefinedOn(FLUIDMECH_VELOCITY);
    velocity->entryType = ResultInfo::VECTOR;
    feFunctions_[FLUIDMECH_VELOCITY]->SetResultInfo(velocity);
    results_.Push_back( velocity );
    availResults_.insert( velocity );

    velocity->SetFeFunction(feFunctions_[FLUIDMECH_VELOCITY]);
    DefineFieldResult( feFunctions_[FLUIDMECH_VELOCITY], velocity );

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    idbcSolNameMap_[FLUIDMECH_PRESSURE] = "pressure";
    idbcSolNameMap_[FLUIDMECH_VELOCITY] = "velocity";
  
    hdbcSolNameMap_[FLUIDMECH_PRESSURE] = "noPressure";
    hdbcSolNameMap_[FLUIDMECH_VELOCITY] = "noSlip";
  }
  
  
  void FlowPDE::DefinePostProcResults() {

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
      stress->definedOn = ResultInfo::MapSolTypeToDefinedOn(FLUIDMECH_STRESS);
      stress->SetFeFunction(feFunctions_[FLUIDMECH_VELOCITY]);
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
      strain->definedOn = ResultInfo::MapSolTypeToDefinedOn(FLUIDMECH_STRAINRATE);
      strain->SetFeFunction(feFunctions_[FLUIDMECH_VELOCITY]);
      availResults_.insert( strain );
      shared_ptr<CoefFunctionFormBased> strainFunc;
      if( isComplex_ ) {
        strainFunc.reset(new CoefFunctionBOp<Complex>(velFeFct, strain));
      } else {
        strainFunc.reset(new CoefFunctionBOp<Double>(velFeFct, strain));
      }
      DefineFieldResult( strainFunc, strain );

  }
  
  void FlowPDE::FinalizePostProcResults() {

    SinglePDE::FinalizePostProcResults();

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
      sigmaFunc->AddIntegrator(
        GetStiffIntegrator( actSDMat, region, isComplex_ ), 
        region
        );
      strainFunc->AddIntegrator(
        GetStiffIntegrator( actSDMat, region, isComplex_ ),
        region
        );
    }

  }

  
  std::map<SolutionType, shared_ptr<FeSpace> >
  FlowPDE::CreateFeSpaces(const std::string& formulation,
                                   PtrParamNode infoNode) 
  {
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;

    std::cout << "Creating FESpaces, formulation " << formulation << std::endl;

    if(formulation == "default" || formulation == "H1"){
      PtrParamNode spaceNode;
      spaceNode = infoNode->Get(SolutionTypeEnum.ToString(FLUIDMECH_VELOCITY));

      crSpaces[FLUIDMECH_VELOCITY] =
        FeSpace::CreateInstance(myParam_,spaceNode,FeSpace::H1, ptGrid_);
      crSpaces[FLUIDMECH_VELOCITY]->Init(solStrat_);

      spaceNode = infoNode->Get(SolutionTypeEnum.ToString(FLUIDMECH_PRESSURE));

      crSpaces[FLUIDMECH_PRESSURE] =
        FeSpace::CreateInstance(myParam_,spaceNode,FeSpace::H1, ptGrid_);
      crSpaces[FLUIDMECH_PRESSURE]->Init(solStrat_);

    } else if ( formulation == "L2" ) {
      std::string form = SolutionTypeEnum.ToString(FLUIDMECH_VELOCITY);
      PtrParamNode potSpaceNode = infoNode->Get(form);
      crSpaces[FLUIDMECH_VELOCITY] =
          FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);

      form = SolutionTypeEnum.ToString(FLUIDMECH_PRESSURE);
      potSpaceNode = infoNode->Get(form);
      crSpaces[FLUIDMECH_PRESSURE] =
              FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::L2, ptGrid_);

      crSpaces[FLUIDMECH_VELOCITY]->Init(solStrat_);
      crSpaces[FLUIDMECH_PRESSURE]->Init(solStrat_);
    }
    else{
      EXCEPTION("The formulation " << formulation << 
                "  of fluidmech PDE is not known!");
    }
    return crSpaces;
  }


  BaseBDBInt *
  FlowPDE::GetStiffIntegrator( BaseMaterial* actSDMat,
                                        RegionIdType regionId,
                                        bool isComplex )
  {

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

    integ->SetName("LinElastInt");
    integ->SetFeSpace( feFunctions_[FLUIDMECH_VELOCITY]->GetFeSpace() );

    return integ;
  }


}
