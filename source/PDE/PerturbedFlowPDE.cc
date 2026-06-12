// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <set>

#include "PerturbedFlowPDE.hh"

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
#include "Forms/Operators/StrainOperator.hh"

#include "Domain/CoefFunction/CoefXpr.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
//#include "Domain/Results/ExternalFieldFunctors.hh"

namespace CoupledField {

  DEFINE_LOG(fluidmechpertpde, "pde.fluidmechpert")

  // ***************
  //   Constructor
  // ***************
  PerturbedFlowPDE::PerturbedFlowPDE( Grid* grid, PtrParamNode paramNode,
                                      PtrParamNode infoNode,
                                      shared_ptr<SimState> simState, 
                                      Domain* domain )
  :SinglePDE( grid, paramNode, infoNode,simState, domain )
  {
    // ======================================================================
    // set solution information
    // ======================================================================
    pdename_          = "fluidMechPerturbed";
    pdematerialclass_ = FLOW;
 
    nonLin_    = false;
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
    stabilized_ = false;
    stabilized_ = myParam_->Get("stabilization")->As<bool>();
    if ( stabilized_ )
    	std::cerr << "\n DO STABILIZED FORMULATION!!\n" << std::endl;

    stabilizedBochev_ = false;
    stabilizedBochev_ = myParam_->Get("stabilizationBochev")->As<bool>();
    if ( stabilizedBochev_ )
          std::cerr << "\n DO STABILIZED FORMULATION OF TYPE BOCHEV!!\n" << std::endl;

    if ( stabilizedBochev_ && stabilized_)
      EXCEPTION("Standad stabilization AND stabilization of Bochev not possible");
    
    enableC2_ = false;
    enableC2_ = myParam_->Get("enableC2")->As<bool>();
    if ( enableC2_ )
          std::cerr << "\n ENABLING SECOND CONVECTIVE TERM 2 (C2)\n" << std::endl;

    factorC1_ = 2.0;
    factorC1_ = myParam_->Get("factorC1")->As<Double>();
  }
  
  void PerturbedFlowPDE::InitNonLin() {
    SinglePDE::InitNonLin();
  }

  void PerturbedFlowPDE::DefineIntegrators() {
    // Not needed at the moment. Commented out due to gcc 4.6.
#if 0
    //transform the type
    SubTensorType tensorType;

    if ( dim_ == 3 ) {
      tensorType = FULL;
    }
    else {
      if ( isaxi_ == true ) {
        tensorType = AXI;
      }
      else {
        // 2d: plane case
        tensorType = PLANE_STRAIN;
      }
    }

    BaseMaterial * actMat = NULL;
#endif

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
      meanVelFct->AddEntityList( actSDList );

      // --------------------------
      //  Set region approximation
      // --------------------------
      // We hardcode the Taylor-Hood spaces for the moment
      velSpace->SetRegionApproximation(actRegion, "velPolyId", "velIntegId");
      //velSpace->SetRegionApproximation(actRegion, "default", "default");
      meanVelSpace->SetRegionApproximation(actRegion, "velPolyId", "velIntegId");
      // meanVelSpace->SetRegionApproximation(actRegion, "default", "default");
      presSpace->SetRegionApproximation(actRegion, "presPolyId", "presIntegId");

      PtrCoefFct density = materials_[actRegion]->GetScalCoefFnc(
        DENSITY,
        Global::REAL
        );
      PtrCoefFct viscosity = materials_[actRegion]->GetScalCoefFnc(
        FLUID_DYNAMIC_VISCOSITY,
        Global::REAL
        );

      WARN("fluid density: " << density->ToString() << "\n " <<
           "dynamic viscosity: " << viscosity->ToString());

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
      stiffIntPV->SetName("PerturbedStiffIntPV");
      BiLinFormContext *stiffContPV = NULL;
      stiffContPV = new BiLinFormContext(stiffIntPV, STIFFNESS );

      stiffContPV->SetEntities( actSDList, actSDList );
      stiffContPV->SetFeFunctions( presFct, velFct);
      assemble_->AddBiLinearForm( stiffContPV );

#ifdef NOT_PARTIALLY_INTEGRATED
      // --------------------------------------------------------------------
      //  VERSION 2: K_VP Integrator (upper off-diagonal integrator)
      // --------------------------------------------------------------------
      PtrCoefFct coeffKVP
                = CoefFunction::Generate(mp_,Global::REAL, "1.0");
      BiLinearForm * stiffIntVP = NULL;
      if( dim_ == 2 ) {
        stiffIntVP = new ABInt< IdentityOperator<FeH1,2,2>,
                                GradientOperator<FeH1,2> > (coeffKVP,
                                                            1.0 / density );
      } else {
        stiffIntVP = new ABInt< IdentityOperator<FeH1,3,3>,
                                GradientOperator<FeH1,3> > (coeffKVP,
                                                            1.0 / density );
      }
      stiffIntVP->SetName("PerturbedStiffIntVP");
      BiLinFormContext *stiffContVP = NULL;
      stiffContVP = new BiLinFormContext(stiffIntVP, STIFFNESS );

      stiffContVP->SetEntities( actSDList, actSDList );
      stiffContVP->SetFeFunctions( velFct, presFct );
      assemble_->AddBiLinearForm( stiffContVP );
#endif


      // --------------------------------------------------------------------
      //  VERSION 2: K_VP Integrator
      //  (upper off-diagonal integrators - partially integrated, volume)
      // --------------------------------------------------------------------
      PtrCoefFct coeffKVP = CoefFunction::Generate(mp_,Global::REAL, "1.0");
      BiLinearForm * stiffIntVP = NULL;
      if( dim_ == 2 ) {
        stiffIntVP = new ABInt<>( new DivOperator<FeH1,2>(),
                                  new MultiIdOp<FeH1,2>(), coeffKVP, -1.0 );
      } else {
        stiffIntVP = new ABInt<>( new DivOperator<FeH1,3>(),
                                  new MultiIdOp<FeH1,3>(), coeffKVP, -1.0 );
      }
      stiffIntVP->SetName("PerturbedStiffIntVP");
      BiLinFormContext *stiffContVP = NULL;
      stiffContVP = new BiLinFormContext(stiffIntVP, STIFFNESS );

      stiffContVP->SetEntities( actSDList, actSDList );
      stiffContVP->SetFeFunctions( velFct, presFct );
      assemble_->AddBiLinearForm( stiffContVP );

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
      stiffIntLaplace->SetName("PerturbedStiffIntViscous");
      BiLinFormContext *stiffContLaplace;
      stiffContLaplace = new BiLinFormContext(stiffIntLaplace, STIFFNESS );

      stiffContLaplace->SetEntities( actSDList, actSDList );
      stiffContLaplace->SetFeFunctions( velFct, velFct );
      assemble_->AddBiLinearForm( stiffContLaplace );


      //======================================================================
      // CHECK FOR FLOW
      //======================================================================
      std::string flowId = curRegNode->Get("flowId")->As<std::string>();
      if(flowId != ""){

        // Get result info object for flow
        shared_ptr<ResultInfo> flowInfo;
        flowInfo = GetResultInfo(MEAN_FLUIDMECH_VELOCITY);
        // Read coefficient flow coefficient function for this region
        bool coefUpdateGeo = false;
        PtrCoefFct regionFlow;
        std::set<UInt> definedDofs;
        
        //Add the region information
        PtrParamNode flowNode = 
          myParam_->Get("flowList")->GetByVal("flow",
                                              "name",
                                              flowId.c_str());
        
        ReadUserFieldValues( actSDList, flowNode, flowInfo->dofNames,
                             flowInfo->entryType, isComplex_, regionFlow,
                             definedDofs, coefUpdateGeo );
        meanFlowCoef_->AddRegion( actRegion, regionFlow );
          
        
        meanVelFct->AddEntityList( actSDList );
        meanVelFct->AddExternalDataSource( regionFlow,
                                           actSDList );

        //now create the integrators
        BiLinearForm *convectiveVv = NULL;

        if( dim_ == 2 ) {
          if(isComplex_) 
          {
            convectiveVv = new ABInt<Complex>( new IdentityOperator<FeH1,2,2>(),
                                               new ConvectiveOperator<FeH1,2,2,Complex>(),
                                               density, factorC1_, coefUpdateGeo );
          }
          else
          {
            convectiveVv = new ABInt<Double>( new IdentityOperator<FeH1,2,2>(),
                                              new ConvectiveOperator<FeH1,2,2>(),
                                              density, factorC1_, coefUpdateGeo );
          }
        } else {
          if(isComplex_) 
          {
            convectiveVv = new ABInt<Complex>( new IdentityOperator<FeH1,3,3>(),
                                               new ConvectiveOperator<FeH1,3,3,Complex>(),
                                               density, factorC1_, coefUpdateGeo );
          }
          else
          {
            convectiveVv = new ABInt<Double>( new IdentityOperator<FeH1,3,3>(),
                                              new ConvectiveOperator<FeH1,3,3>(),
                                              density, factorC1_, coefUpdateGeo );
          }
        }

        convectiveVv->SetBCoefFunctionOpB(meanVelFct);
        
        convectiveVv->SetName("PerturbedStiffIntConvectiveVv");

        BiLinFormContext *convectiveContextVv = NULL;
        convectiveContextVv = new BiLinFormContext(convectiveVv, STIFFNESS );

        convectiveContextVv->SetEntities( actSDList, actSDList );
        convectiveContextVv->SetFeFunctions( feFunctions_[FLUIDMECH_VELOCITY],
                                             feFunctions_[FLUIDMECH_VELOCITY] );
        assemble_->AddBiLinearForm( convectiveContextVv );

        if(enableC2_) 
        {
          // Second convective term. Derivative tensor of mean flow field is a
          // factor computed by a CoefFunction.
          BaseBOperator* bOpGrad;
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
            if( dim_ == 2 ) {
              bOpGrad = new GradientOperator<FeH1,2, 1, Complex>();
              coeffConvec.reset(
                new CoefFunctionMeanFlowConvection<Complex,2>( density, bOpGrad, meanVelFct )
                );
            } else {
              bOpGrad = new GradientOperator<FeH1,3, 1, Complex>();
              coeffConvec.reset(
                new CoefFunctionMeanFlowConvection<Complex,3>( density, bOpGrad, meanVelFct )
                );
            }          
            convectivevV = new BDBInt<Complex,Complex>( bOpId, coeffConvec, 1.0 );
            
          } else {
            if( dim_ == 2 ) {
              bOpGrad = new GradientOperator<FeH1,2, 1, Double>();
              coeffConvec.reset(
                new CoefFunctionMeanFlowConvection<Double,2>( density, bOpGrad, meanVelFct )
                );
            } else {
              bOpGrad = new GradientOperator<FeH1,3, 1, Double>();
              coeffConvec.reset(
                new CoefFunctionMeanFlowConvection<Double,3>( density, bOpGrad, meanVelFct )
                );
            }
            convectivevV = new BDBInt<Double,Double>( bOpId, coeffConvec, 1.0 );
          }
          
          convectivevV->SetName("PerturbedStiffIntConvectivevV");
          
          BiLinFormContext *convectiveContextvV = NULL;
          convectiveContextvV = new BiLinFormContext(convectivevV, STIFFNESS );
          
          convectiveContextvV->SetEntities( actSDList, actSDList );
          convectiveContextvV->SetFeFunctions( feFunctions_[FLUIDMECH_VELOCITY],
                                               feFunctions_[FLUIDMECH_VELOCITY] );
          assemble_->AddBiLinearForm( convectiveContextvV );
        }
        
        
        //======================================================================
        // DO STABILIZATION due to convection
        //======================================================================
        if ( stabilized_) {
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
          PtrCoefFct coeffConvecStab;
          coeffConvecStab.reset(
            new CoefFunctionStabParams( density, viscosity, meanFlowCoef_,
                                        bOpGrad, presFct,
                                        CoefFunctionStabParams::SUPG, isComplex_ )
            );

          //stabilization: velocity mass matrix
          BiLinearForm *convecMassVVstab = NULL;
          if( dim_ == 2 ) {
            if(isComplex_) 
            {
              convecMassVVstab = new ABInt<Complex>(
                new ConvectiveOperator<FeH1,2,2,Complex>(),
                new IdentityOperator<FeH1,2,2>(),
                coeffConvecStab, densityVal, coefUpdateGeo
                );
            }
            else
            {
              convecMassVVstab = new ABInt<Double>(
                new ConvectiveOperator<FeH1,2,2>(),
                new IdentityOperator<FeH1,2,2>(),
                coeffConvecStab, densityVal, coefUpdateGeo
                );
            }
          } else {
            if(isComplex_) 
            {
              convecMassVVstab = new ABInt<Complex>(
                new ConvectiveOperator<FeH1,3,3,Complex>(),
                new IdentityOperator<FeH1,3,3>(),
                coeffConvecStab, densityVal, coefUpdateGeo
                );
            }
            else
            {
              convecMassVVstab = new ABInt<>(
                new ConvectiveOperator<FeH1,3,3>(),
                new IdentityOperator<FeH1,3,3>(),
                coeffConvecStab, densityVal, coefUpdateGeo
                );
            }
          }
          convecMassVVstab->SetBCoefFunctionOpA(meanFlowCoef_);
          convecMassVVstab->SetName("PerturbedDampIntVVConvectiveStab");

          BiLinFormContext *dampContextVVStab = NULL;
          dampContextVVStab = new BiLinFormContext( convecMassVVstab,
                                                    DAMPING );
          dampContextVVStab->SetEntities( actSDList, actSDList );
          dampContextVVStab->SetFeFunctions( velFct, velFct);
          assemble_->AddBiLinearForm( dampContextVVStab );

          //stabilization: velocity stiffness matrix
          BiLinearForm *convecStiffVVstab = NULL;
          if( dim_ == 2 ) {
            if(isComplex_)
            {
              convecStiffVVstab = new BBInt<Complex>(
                new ConvectiveOperator<FeH1,2,2,Complex>(),
                coeffConvecStab, densityVal, coefUpdateGeo
                );
            }
            else
            {
              convecStiffVVstab = new BBInt<>(
                new ConvectiveOperator<FeH1,2,2>(),
                coeffConvecStab, densityVal, coefUpdateGeo
                );
            }
          } else {
            if(isComplex_) 
            {
              convecStiffVVstab = new BBInt<Complex>(
                new ConvectiveOperator<FeH1,3,3,Complex>(),
                coeffConvecStab, densityVal, coefUpdateGeo
                );
            }
            else
            {
              convecStiffVVstab = new BBInt<>(
                new ConvectiveOperator<FeH1,3,3>(),
                coeffConvecStab, densityVal, coefUpdateGeo
                );
            }
          }
          convecStiffVVstab->SetBCoefFunctionOpB(meanFlowCoef_);
          convecStiffVVstab->SetName("PerturbedStiffIntVVConvectiveStab");

          BiLinFormContext *stiffContextVVStab = NULL;
          stiffContextVVStab = new BiLinFormContext( convecStiffVVstab,
                                                     STIFFNESS );
          stiffContextVVStab->SetEntities( actSDList, actSDList );
          stiffContextVVStab->SetFeFunctions( velFct, velFct);
          assemble_->AddBiLinearForm( stiffContextVVStab );

          //stabilization: pressure stiffness matrix, in velocity equation
          BiLinearForm *convecStiffVPstab = NULL;
          if( dim_ == 2 ) {
            if(isComplex_) 
            {
              convecStiffVPstab = new ABInt<Complex>(
                new ConvectiveOperator<FeH1,2,2,Complex>(),
                new GradientOperator<FeH1,2>(),
                coeffConvecStab, 1.0, coefUpdateGeo
                );
            }
            else
            {
              convecStiffVPstab = new ABInt<>(
                new ConvectiveOperator<FeH1,2,2>(),
                new GradientOperator<FeH1,2>(),
                coeffConvecStab, 1.0, coefUpdateGeo
                );
            }
          } else {
            if(isComplex_) 
            {
              convecStiffVPstab = new ABInt<Complex>(
                new ConvectiveOperator<FeH1,3,3,Complex>(),
                new GradientOperator<FeH1,3>(),
                coeffConvecStab, 1.0, coefUpdateGeo
                );
            }
            else
            {
              convecStiffVPstab = new ABInt<>(
                new ConvectiveOperator<FeH1,3,3>(),
                new GradientOperator<FeH1,3>(),
                coeffConvecStab, 1.0, coefUpdateGeo
                );
            }
          }
          
          convecStiffVPstab->SetBCoefFunctionOpA(meanFlowCoef_);
          convecStiffVPstab->SetName("PerturbedStiffIntVPConvectiveStab");

          BiLinFormContext *stiffContextVPStab = NULL;
          stiffContextVPStab = new BiLinFormContext( convecStiffVPstab,
                                                     STIFFNESS );
          stiffContextVPStab->SetEntities( actSDList, actSDList );
          stiffContextVPStab->SetFeFunctions( velFct, presFct);
          assemble_->AddBiLinearForm( stiffContextVPStab );

          //pressure stabilization: velocity stiffness matrix, in pressure equation
          PtrCoefFct coeffPressStab;
          coeffPressStab.reset(
            new CoefFunctionStabParams(density, viscosity, meanFlowCoef_,
                                       bOpGrad, presFct,
                                       CoefFunctionStabParams::PSPG, isComplex_ )
            );

          BiLinearForm *convecStiffPVstab = NULL;
          if( dim_ == 2 ) {
            if(isComplex_) 
            {
              convecStiffPVstab = new ABInt<Complex>( new GradientOperator<FeH1,2>(),
                                                      new ConvectiveOperator<FeH1,2,2,Complex>(),
                                                      coeffPressStab, densityVal,
                                                      coefUpdateGeo);
            }
            else
            {
              convecStiffPVstab = new ABInt<>( new GradientOperator<FeH1,2>(),
                                               new ConvectiveOperator<FeH1,2,2>(),
                                               coeffPressStab, densityVal,
                                               coefUpdateGeo);
            }
          } else {
            if(isComplex_) 
            {
              convecStiffPVstab = new ABInt<Complex>( new GradientOperator<FeH1,3>(),
                                                      new ConvectiveOperator<FeH1,3,3,Complex>(),
                                                      coeffPressStab, densityVal,
                                                      coefUpdateGeo);
            }
            else
            {
              convecStiffPVstab = new ABInt<>( new GradientOperator<FeH1,3>(),
                                               new ConvectiveOperator<FeH1,3,3>(),
                                               coeffPressStab, densityVal,
                                               coefUpdateGeo);
            }
          }
          convecStiffPVstab->SetBCoefFunctionOpB(meanFlowCoef_);
          convecStiffPVstab->SetName("PerturbedStiffIntPVConvectiveStab");

          BiLinFormContext *stiffContextPVStab = NULL;
          stiffContextPVStab = new BiLinFormContext( convecStiffPVstab,
                                                     STIFFNESS );
          stiffContextPVStab->SetEntities( actSDList, actSDList );
          stiffContextPVStab->SetFeFunctions( presFct, velFct );
          assemble_->AddBiLinearForm( stiffContextPVStab );

          //stabilization: Least Squares InCompressibility
          PtrCoefFct coeffLsicStab;
          coeffLsicStab.reset(
            new CoefFunctionStabParams( density, viscosity, meanFlowCoef_,
                                        bOpGrad, presFct,
                                        CoefFunctionStabParams::LSIC, isComplex_ )
            );
          BiLinearForm *LsicStiffVVstab = NULL;
          if( dim_ == 2 ) {
            LsicStiffVVstab = new BBInt<>( new DivOperator<FeH1,2>(),
                                           coeffConvecStab,
                                           densityVal,
                                           coefUpdateGeo);
          } else {
            LsicStiffVVstab = new BBInt<>( new DivOperator<FeH1,3>(),
                                           coeffConvecStab,
                                           densityVal,
                                           coefUpdateGeo);
          }
          LsicStiffVVstab->SetName("PerturbedLSICStiffIntVVStab");

          BiLinFormContext *stiffContextLsicStab = NULL;
          stiffContextLsicStab =new BiLinFormContext( LsicStiffVVstab,
                                                      STIFFNESS );
          stiffContextLsicStab->SetEntities( actSDList, actSDList );
          stiffContextLsicStab->SetFeFunctions( velFct, velFct);
          assemble_->AddBiLinearForm( stiffContextLsicStab );
        } //stabilized
      } // is flow
      
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
      dampIntvv->SetName("PerturbedDampInt");

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
        //stabilization of pressure
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
        stiffIntPPstab->SetName("PerturbedStiffIntPPstab");
        BiLinFormContext *stiffContPPstab = NULL;
        stiffContPPstab = new BiLinFormContext(stiffIntPPstab, STIFFNESS );
        
        stiffContPPstab->SetEntities( actSDList, actSDList );
        stiffContPPstab->SetFeFunctions( presFct, presFct );
        assemble_->AddBiLinearForm( stiffContPPstab );
      }

      if ( stabilized_) {
        PtrCoefFct coeffKPPstab;
        coeffKPPstab.reset(
          new CoefFunctionStabParams( density,
                                      viscosity,
                                      CoefFunctionStabParams::PSPG, isComplex_ )
          );
        BiLinearForm * stiffIntPPstab = NULL;
        if( dim_ == 2 ) {
          stiffIntPPstab = new BBInt<>( new GradientOperator<FeH1,2>(),
                                        coeffKPPstab,1.0);
        } else {
          stiffIntPPstab = new BBInt<>( new GradientOperator<FeH1,3>(),
                                        coeffKPPstab,1.0);
        }
        stiffIntPPstab->SetName("PerturbedStiffIntPPstab");
        BiLinFormContext *stiffContPPstab = NULL;
        stiffContPPstab = new BiLinFormContext(stiffIntPPstab, STIFFNESS );

        stiffContPPstab->SetEntities( actSDList, actSDList );
        stiffContPPstab->SetFeFunctions( presFct, presFct );
        assemble_->AddBiLinearForm( stiffContPPstab );

        //stabilization of pressure
        Double densityVal;
        LocPointMapped map;
        density->GetScalar(densityVal, map);
        BiLinearForm *dampIntpvStab = NULL;
        if( dim_ == 2 ) {
          dampIntpvStab = new ABInt<>( new GradientOperator<FeH1,2>(),
                                       new IdentityOperator<FeH1,2,2>(),
                                       coeffKPPstab, densityVal );
        } else {
          dampIntpvStab = new ABInt<>( new GradientOperator<FeH1,3>(),
                                       new IdentityOperator<FeH1,3,3>(),
                                       coeffKPPstab, densityVal );
        }
        dampIntpvStab->SetName("PerturbedDampIntStab");
        BiLinFormContext *dampContextpvStab = NULL;
        dampContextpvStab = new BiLinFormContext( dampIntpvStab, DAMPING );

        dampContextpvStab->SetEntities( actSDList, actSDList );
        dampContextpvStab->SetFeFunctions( presFct, velFct);
        assemble_->AddBiLinearForm( dampContextpvStab );
      }      
    }
    
    StdVector<RegionIdType>::iterator surfIt, surfEnd;
    surfIt = presSurfaces_.Begin();
    surfEnd = presSurfaces_.End();

    for( ; surfIt != surfEnd; surfIt++ ) {
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( *surfIt );

      velFct->AddEntityList( actSDList );
      velSpace->SetRegionApproximation(*surfIt, "velPolyId", "velIntegId");
      presFct->AddEntityList( actSDList );
      presSpace->SetRegionApproximation(*surfIt, "presPolyId", "presIntegId");
      meanVelFct->AddEntityList( actSDList );
      meanVelSpace->SetRegionApproximation(*surfIt, "velPolyId", "velIntegId");

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
      stiffIntVPSurf->SetName("PerturbedStiffIntVPSurf");
      BiLinFormContext *stiffContVP = NULL; 
      stiffContVP = new BiLinFormContext(stiffIntVPSurf, STIFFNESS );

      stiffContVP->SetEntities( actSDList, actSDList );
      stiffContVP->SetFeFunctions( velFct, presFct );
      assemble_->AddBiLinearForm( stiffContVP );
    }

    // Since we manage the mean flow FeFunction, we also have to finalize here.
    // c.f. SinglePDE::InitStage3
    meanVelSpace->Finalize();
    meanVelSpace->PreCalcShapeFncs();
    meanVelFct->Finalize();

    //    meanVelSpace->PrintEqnMap();

  }
  
  void PerturbedFlowPDE::DefineSolveStep() {
    solveStep_ = new SolveStepElec(*this);
  }

  void PerturbedFlowPDE::ReadSpecialBCs( ) 
  {
  }
  
  void PerturbedFlowPDE::DefinePrimaryResults() {
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
    velocity->unit = MapSolTypeToUnit(velocity->resultType);

    velocity->definedOn = ResultInfo::MapSolTypeToDefinedOn(FLUIDMECH_VELOCITY);
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
  
  
  void PerturbedFlowPDE::DefinePostProcResults() {

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

  #if 0
      // === PRESSURE GRADIENT (just for debugging purposes) ===
      shared_ptr<ResultInfo> ef ( new ResultInfo );
      ef->resultType = FLUIDMECH_PRES_GRADIENT;
      ef->SetVectorDOFs(dim_, isaxi_);
      ef->unit = MapSolTypeToUnit(FLUIDMECH_PRES_GRADIENT);
      ef->definedOn = ResultInfo::MapSolTypeToDefinedOn(FLUIDMECH_PRES_GRADIENT);
      ef->entryType = ResultInfo::VECTOR;
      availResults_.insert( ef );
      shared_ptr<CoefFunctionFormBased> eFunc;
      if( isComplex_ ) {
        eFunc.reset(new CoefFunctionBOp<Complex>(feFct, ef, -1.0));
      } else {
        eFunc.reset(new CoefFunctionBOp<Double>(feFct, ef, -1.0));
      }

        PtrCoefFct coeffKVP
                  = CoefFunction::Generate(mp_, Global::REAL, "1.0");
        BaseBDBInt * stiffIntVP = NULL;
        if( dim_ == 2 ) {
          stiffIntVP = new ABInt<>(new IdentityOperator<FeH1,2,2>(), 
                                   new GradientOperator<FeH1,2>(), coeffKVP, 1.0 );
        } else {
          stiffIntVP = new ABInt<>(new IdentityOperator<FeH1,3,3>(), 
                                   new GradientOperator<FeH1,3>(), coeffKVP, 1.0 );
        }
        
      DefineFieldResult( eFunc, ef );
  #endif

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
  
  void PerturbedFlowPDE::FinalizePostProcResults() {

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
  PerturbedFlowPDE::CreateFeSpaces(const std::string& formulation,
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


  void PerturbedFlowPDE::CreateMeanFlowFunction(StdVector<std::string> dofNames)
  {
    //// === MEAN FLUIDMECH VELOCITY ===
    shared_ptr<ResultInfo> flowvelocity( new ResultInfo );
    flowvelocity->resultType = MEAN_FLUIDMECH_VELOCITY;
    flowvelocity->dofNames = dofNames;
    flowvelocity->unit = "m/s";

    flowvelocity->definedOn = ResultInfo::MapSolTypeToDefinedOn(MEAN_FLUIDMECH_VELOCITY);
    flowvelocity->entryType = ResultInfo::VECTOR;
    results_.Push_back( flowvelocity );
    availResults_.insert( flowvelocity );
    
    meanFlowCoef_.reset(
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

    DefineFieldResult(meanFlowFeFct_, flowvelocity);
  }

  BaseBDBInt *
  PerturbedFlowPDE::GetStiffIntegrator( BaseMaterial* actSDMat,
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
