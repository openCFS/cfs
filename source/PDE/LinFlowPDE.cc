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


#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
//#include "Domain/Results/ExternalFieldFunctors.hh"

namespace CoupledField {

  DECLARE_LOG(linfluidmechpde)
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
    if ( enableC2_ )
          std::cerr << "\n ENABLING SECOND CONVECTIVE TERM 2 (C2)\n" << std::endl;

    factorC1_ = 2.0;
    factorC1_ = myParam_->Get("factorC1")->As<Double>();
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
      PtrCoefFct viscosity =
    		  materials_[actRegion]->GetScalCoefFnc(DYNAMIC_VISCOSITY, Global::REAL);

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

      if ( isCompressible_ ) {
        // ====================================================================
        // K_PP: stiffness integrator, conservation of mass
        // add time derivative of density expressed by pressure according to
    	  // thermodynamic relation
        // ====================================================================
    	  PtrCoefFct refPres = materials_[actRegion]->GetScalCoefFnc(
    	      			  REF_PRESSURE, Global::REAL);
    	  PtrCoefFct fnc;
    	  if ( isHeatCoupled_ ) {
    		  fnc = CoefFunction::Generate( mp_, Global::REAL,
    		  				  CoefXprBinOp(mp_, constOne, refPres, CoefXpr::OP_DIV ) );
    	  }
    	  else {
    		  PtrCoefFct adiabaticExp = materials_[actRegion]->GetScalCoefFnc(
    			  ADIABATIC_EXPONENT, Global::REAL);
    		  PtrCoefFct help1 = CoefFunction::Generate( mp_,  Global::REAL,
    			  CoefXprBinOp(mp_, adiabaticExp, refPres, CoefXpr::OP_MULT ) );
    		  fnc = CoefFunction::Generate( mp_, Global::REAL,
				  CoefXprBinOp(mp_, constOne, help1, CoefXpr::OP_DIV ) );
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
      //  (upper off-diagonal integrators - partially integrated, volume)
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
      //  Vector Laplace term
      // ====================================================================
      BiLinearForm * stiffIntLaplace = NULL;
      if( dim_ == 2 ) {
        stiffIntLaplace = new BBInt<>( new LaplOperator<FeH1,2>(),
                                       viscosity, 1.0 );
      } else {
        stiffIntLaplace = new BBInt<>( new LaplOperator<FeH1,3>(),
                                       viscosity, 1.0 );
      }
      stiffIntLaplace->SetName("LinFlowStiffIntViscous");
      BiLinFormContext *stiffContLaplace;
      stiffContLaplace = new BiLinFormContext(stiffIntLaplace, STIFFNESS );

      stiffContLaplace->SetEntities( actSDList, actSDList );
      stiffContLaplace->SetFeFunctions( velFct, velFct );
      assemble_->AddBiLinearForm( stiffContLaplace );

      if ( isCompressible_ ) {
        // ====================================================================
        // K_VV Integrator: conservation of momentum
        //  diagonal integrator - partially integrated
        //  Grad Div term
        // ====================================================================
    	  PtrCoefFct bulkViscosity1 = materials_[actRegion]->GetScalCoefFnc(BULK_VISCOSITY, Global::REAL);
          // we will add -2/3 of dynamic viscosity to bulk viscosity in order to only have viscous effect
          // on shear layer
          // 1/3 is a coef for total bulk viscosity, therefore we ignore that here and only add -2

          PtrCoefFct coeftwo = CoefFunction::Generate( mp_, Global::REAL, "-2.0");
          PtrCoefFct coeftwodv = CoefFunction::Generate( mp_,  Global::REAL, CoefXprBinOp(mp_, coeftwo, viscosity, CoefXpr::OP_MULT ) );
          PtrCoefFct bulkViscosity = CoefFunction::Generate( mp_,  Global::REAL, CoefXprBinOp(mp_, coeftwodv, bulkViscosity1, CoefXpr::OP_ADD ) );
    	  PtrCoefFct coef = CoefFunction::Generate( mp_,  Global::REAL,
    	      			  CoefXprBinOp(mp_, viscosity, bulkViscosity, CoefXpr::OP_ADD ) );

    	  BiLinearForm * stiffIntDivDiv = NULL;
    	  if( dim_ == 2 ) {
    		  stiffIntDivDiv = new BBInt<>(new ScalarDivergenceOperator<FeH1,2,Double>(),
    				                       coef, 1.0, updatedGeo_);
    	  } else {
    		  stiffIntDivDiv = new BBInt<>(new ScalarDivergenceOperator<FeH1,3,Double>(),
    				                       coef, 1.0, updatedGeo_);
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
        
        convectiveVv->SetName("LinFlowStiffIntConvectiveVv");

        BiLinFormContext *convectiveContextVv = NULL;
        convectiveContextVv = new BiLinFormContext(convectiveVv, STIFFNESS );

        convectiveContextVv->SetEntities( actSDList, actSDList );
        convectiveContextVv->SetFeFunctions( velFct, velFct );
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
                new CoefFunctionMeanFlowConvection<Complex,2>( density, viscosity,
                                                               bOpGrad, meanVelFct )
                );
            } else {
              bOpGrad = new GradientOperator<FeH1,3, 1, Complex>();
              coeffConvec.reset(
                new CoefFunctionMeanFlowConvection<Complex,3>( density, viscosity,
                                                               bOpGrad, meanVelFct )
                );
            }          
            convectivevV = new BDBInt<Complex,Complex>( bOpId, coeffConvec, 1.0 );
            
          } else {
            if( dim_ == 2 ) {
              bOpGrad = new GradientOperator<FeH1,2, 1, Double>();
              coeffConvec.reset(
                new CoefFunctionMeanFlowConvection<Double,2>( density, viscosity,
                                                              bOpGrad, meanVelFct )
                );
            } else {
              bOpGrad = new GradientOperator<FeH1,3, 1, Double>();
              coeffConvec.reset(
                new CoefFunctionMeanFlowConvection<Double,3>( density, viscosity,
                                                              bOpGrad, meanVelFct )
                );
            }
            convectivevV = new BDBInt<Double,Double>( bOpId, coeffConvec, 1.0 );
          }
          
          convectivevV->SetName("LinFlowStiffIntConvectivevV");
          
          BiLinFormContext *convectiveContextvV = NULL;
          convectiveContextvV = new BiLinFormContext(convectivevV, STIFFNESS );
          
          convectiveContextvV->SetEntities( actSDList, actSDList );
          convectiveContextvV->SetFeFunctions( velFct, velFct );
          assemble_->AddBiLinearForm( convectiveContextvV );
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
    shared_ptr<FeSpace> mySpace = velFct->GetFeSpace();

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

    curCoef = regionMat->GetTensorCoefFnc( DYNAMIC_VISCOSITY, subTensorType, complexPart );
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
