// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
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
#include "Utils/StdVector.hh"
#include "Driver/SolveSteps/SolveStepElec.hh"
#include "CoupledPDE/PDECoupling.hh"
#include "Driver/Assemble.hh"
#include "Utils/ApproxData.hh"
#include "Utils/SmoothSpline.hh"
#include "Materials/Models/Hysteresis.hh"
#include "Materials/Models/Preisach.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/FeFunctions.hh"

// new integrator concept
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/BUInt.hh"

#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/DivOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/MultiIdOp.hh"
#include "Forms/Operators/LaplOp.hh"
#include "Forms/Operators/ConvectiveOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"
#include "Forms/Operators/StrainOperator.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
//#include "Domain/Results/ExternalFieldFunctors.hh"

namespace CoupledField {

  DECLARE_LOG(fluidmechpertpde)
  DEFINE_LOG(fluidmechpertpde, "pde.fluidmechpert")

  // ***************
  //   Constructor
  // ***************
  PerturbedFlowPDE::PerturbedFlowPDE( Grid* grid, PtrParamNode paramNode )
    :SinglePDE( grid, paramNode ) {


    // ======================================================================
    // set solution information
    // ======================================================================
    pdename_          = "fluidMechPerturbed";
    pdematerialclass_ = FLOW;
    maxTimeDerivOrder_ = 1;
 
    nonLin_    = false;
    nonLinMaterial_ = false;
    isAlwaysStatic_ = false;

    needSolPrev_ = true;
 
    // Check the subtype of the problem
    paramNode->GetValue("subType", subType_);

    if(ptGrid_->GetDim() == 3)
      subType_ = "3d";
      
    // Obtain regions the pde is defined on
    PtrParamNode presSurfaceNode = myParam_->Get("presSurfaceList", ParamNode::PASS);
    ParamNodeList regionNodes;
    if(presSurfaceNode) {
      regionNodes = presSurfaceNode->GetList("presSurf");
    }

    // output and set regions_
    for( UInt i = 0; i < regionNodes.GetSize(); i++ )
    {
      RegionIdType actRegionId = ptGrid_->GetRegion().Parse(regionNodes[i]->Get("name")->As<std::string>());

      // Check, if region was already defined an issue a warning otherwise
      if( std::find(presSurfaces_.Begin(), presSurfaces_.End(), actRegionId )
          != presSurfaces_.End() )  {
        WARN( "The pressure surface '" << regionNodes[i]->Get("name")->As<std::string>()
              << "' was already defined for PDE '" << pdename_
              << "'. Please remove duplicate entries." );
      }

      presSurfaces_.Push_back( actRegionId );
    }

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

      // --------------------------
      //  Set region approximation
      // --------------------------
#if 0
      // --- Set the approximation for the current region ---
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
#endif
      // We hardcode the Taylor-Hood spaces for the moment
      velSpace->SetRegionApproximation(actRegion, "velPolyId", "velIntegId");
      presSpace->SetRegionApproximation(actRegion, "presPolyId", "presIntegId");

      PtrCoefFct density = materials_[actRegion]->GetScalCoefFnc( DENSITY, Global::REAL );
      PtrCoefFct viscosity = materials_[actRegion]->GetScalCoefFnc( DYNAMIC_VISCOSITY, Global::REAL );


      WARN("fluid density: " << density->ToString() << "\n" << " dynamic viscosity: " << viscosity->ToString());
      

      // Create set of flow regions and map of density functions for surface integrators.
      flowRegions.insert(actRegion);
      oneFuncs[actRegion] = CoefFunction::Generate(Global::REAL,
                                                       lexical_cast<std::string>(1.0));

      // ====================================================================
      // stiffness integrators
      // ====================================================================
      // --------------------------------------------------------------------
      //  VERSION 1: K_PV Integrator (lower off-diagonal integrator)
      //  \int_{\Omega_f} phi div(v') d\Omega = 0
      // --------------------------------------------------------------------

//      PtrCoefFct coeffKPV
//                = CoefFunction::Generate(Global::REAL, "1.0");
                
      BiLinearForm * stiffIntPV = NULL;
      if( dim_ == 2 ) {
        stiffIntPV = new ABInt<>( new MultiIdOp<FeH1,2>(), 
                                  new DivOperator<FeH1,2>(), density, 1.0 );
      } else {
        stiffIntPV = new ABInt<>(new MultiIdOp<FeH1,3>() , new DivOperator<FeH1,3>(), density, 1.0 );
      }
      stiffIntPV->SetName("PerturbedStiffIntPV");
      //stiffIntPV->SetFeSpace( feFunctions_[ACOU_PRESSURE]->GetFeSpace(), feFunctions_[ACOU_VELOCITY]->GetFeSpace() );
      BiLinFormContext *stiffContPV = new BiLinFormContext(stiffIntPV, STIFFNESS );

      stiffContPV->SetEntities( actSDList, actSDList );
      stiffContPV->SetFeFunctions( presFct, velFct);
//      stiffContPV->SetCounterPart(true);
      assemble_->AddBiLinearForm( stiffContPV );

#ifdef NOT_PARTIALLY_INTEGRATED
      // --------------------------------------------------------------------
      //  VERSION 2: K_VP Integrator (upper off-diagonal integrator)
      // --------------------------------------------------------------------
      PtrCoefFct coeffKVP
                = CoefFunction::Generate(Global::REAL, "1.0");
      BiLinearForm * stiffIntVP = NULL;
      if( dim_ == 2 ) {
        stiffIntVP = new ABInt< IdentityOperator<FeH1,2,2> , GradientOperator<FeH1,2> >(coeffKVP, 1.0 / density );
      } else {
        stiffIntVP = new ABInt< IdentityOperator<FeH1,3,3> , GradientOperator<FeH1,3> >(coeffKVP, 1.0 / density );
      }
      stiffIntVP->SetName("PerturbedStiffIntVP");
      //stiffIntVP->SetFeSpace( feFunctions_[ACOU_PRESSURE]->GetFeSpace(), feFunctions_[ACOU_VELOCITY]->GetFeSpace() );
      BiLinFormContext *stiffContVP = new BiLinFormContext(stiffIntVP, STIFFNESS );

      stiffContVP->SetEntities( actSDList, actSDList );
      stiffContVP->SetFeFunctions( velFct, presFct );
      //stiffContVP->SetCounterPart(true);
      assemble_->AddBiLinearForm( stiffContVP );
#endif

      // --------------------------------------------------------------------
      //  VERSION 2: K_VP Integrator
      //  (upper off-diagonal integrators - partially integrated, volume)
      // --------------------------------------------------------------------
      PtrCoefFct coeffKVP
                = CoefFunction::Generate(Global::REAL, "1.0");
      BiLinearForm * stiffIntVP = NULL;
      if( dim_ == 2 ) {
        stiffIntVP = new ABInt<>(new  DivOperator<FeH1,2>(), new MultiIdOp<FeH1,2>(), coeffKVP, -1.0 );
      } else {
        stiffIntVP = new ABInt<>(new DivOperator<FeH1,3>() , new MultiIdOp<FeH1,3>(), coeffKVP, -1.0 );
      }
      stiffIntVP->SetName("PerturbedStiffIntVP");
      //stiffIntVP->SetFeSpace( feFunctions_[ACOU_PRESSURE]->GetFeSpace(), feFunctions_[ACOU_VELOCITY]->GetFeSpace() );
      BiLinFormContext *stiffContVP = new BiLinFormContext(stiffIntVP, STIFFNESS );

      stiffContVP->SetEntities( actSDList, actSDList );
      stiffContVP->SetFeFunctions( velFct, presFct );
      //stiffContVP->SetCounterPart(true);
      assemble_->AddBiLinearForm( stiffContVP );

      // --------------------------------------------------------------------
      //  VERSION 2: K_Laplace Integrator
      // --------------------------------------------------------------------
//      PtrCoefFct coeffKvv
//                = CoefFunction::Generate(Global::REAL,
//                                         lexical_cast<std::string>(viscosity));
      BiLinearForm * stiffIntLaplace = NULL;
      if( dim_ == 2 ) {
        stiffIntLaplace = new BBInt<>(new LaplOperator<FeH1,2>(), viscosity, 1.0);
      } else {
        stiffIntLaplace = new BBInt<>(new LaplOperator<FeH1,3>(), viscosity, 1.0);
      }
      stiffIntLaplace->SetName("PerturbedStiffIntViscous");
      BiLinFormContext *stiffContLaplace;
      stiffContLaplace = new BiLinFormContext(stiffIntLaplace, STIFFNESS );

      stiffContLaplace->SetEntities( actSDList, actSDList );
      stiffContLaplace->SetFeFunctions( velFct, velFct );
      assemble_->AddBiLinearForm( stiffContLaplace );


//      PtrCoefFct coeffMVV
//                = CoefFunction::Generate(Global::REAL,
//                                         lexical_cast<std::string>(density));

      //======================================================================
      // CHECK FOR FLOW
      //======================================================================
      std::string flowId = curRegNode->Get("flowId")->As<std::string>();
      if(flowId != ""){
        // Get result info object for flow
        shared_ptr<ResultInfo> flowInfo = GetResultInfo(MEAN_FLUIDMECH_VELOCITY); 
        
        //Add the region information
        PtrParamNode flowNode = myParam_->Get("flowList")->GetByVal("flow","name",flowId.c_str());
        
        // Read coefficient flow coefficient function for this region
        PtrCoefFct regionFlow;
        std::set<UInt> definedDofs;
        ReadUserFieldValues( regionName, flowNode, flowInfo->dofNames, flowInfo->entryType, 
                             false, regionFlow, definedDofs );
        meanFlowCoef_->AddRegion( actRegion, regionFlow );
        
//        if(isComplex_){
//          shared_ptr<FieldFunctor<Complex> > fct = dynamic_pointer_cast<FieldFunctor<Complex> >(meanFlowFunctor_);
//          fct->AddRegion(actRegion,flowNode);
//        }else{
//          shared_ptr<FieldFunctor<Double> > fct = dynamic_pointer_cast<FieldFunctor<Double> >(meanFlowFunctor_);
//          fct->AddRegion(actRegion,flowNode);
//        }

        //now create the integrators
        BiLinearForm *convectiveVv = NULL;
        if( dim_ == 2 ) {
          convectiveVv = new ABInt<>(new IdentityOperator<FeH1,2,2>(),
                                     new ConvectiveOperator<FeH1,2,2>(),
                                     density, 1.0);
        } else {
          convectiveVv = new ABInt<>(new IdentityOperator<FeH1,3,3>(),
                                     new ConvectiveOperator<FeH1,3,3>(), density, 1.0);
        }

        convectiveVv->SetBCoefFunctionOpB(meanFlowCoef_);

        convectiveVv->SetName("PerturbedStiffIntConvective");

        BiLinFormContext *convectiveContextVv =  new BiLinFormContext(convectiveVv, STIFFNESS );

        convectiveContextVv->SetEntities( actSDList, actSDList );
        convectiveContextVv->SetFeFunctions( feFunctions_[FLUIDMECH_VELOCITY],feFunctions_[FLUIDMECH_VELOCITY]);
        assemble_->AddBiLinearForm( convectiveContextVv );
      }

      // ====================================================================
      // damping integrators
      // ====================================================================
//      PtrCoefFct coeffDvv
//                = CoefFunction::Generate(Global::REAL, lexical_cast<std::string>(density));

      BiLinearForm *dampIntvv = NULL;
      if( dim_ == 2 ) {
        dampIntvv = new BBInt<>(new IdentityOperator<FeH1,2,2>(), density, 1.0 );
      } else {
        dampIntvv = new BBInt<>(new IdentityOperator<FeH1,3,3>(), density, 1.0 );
      }
      dampIntvv->SetName("PerturbedDampInt");
      //massIntPP->SetFeSpace( feFunctions_[ACOU_PRESSURE]->GetFeSpace() );

      BiLinFormContext *dampContextvv =  new BiLinFormContext(dampIntvv, DAMPING );

      dampContextvv->SetEntities( actSDList, actSDList );
      dampContextvv->SetFeFunctions( feFunctions_[FLUIDMECH_VELOCITY],feFunctions_[FLUIDMECH_VELOCITY]);
      assemble_->AddBiLinearForm( dampContextvv );
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
      stiffIntVPSurf->SetName("PerturbedStiffIntVPSurf");
      //stiffIntVP->SetFeSpace( feFunctions_[ACOU_PRESSURE]->GetFeSpace(), feFunctions_[ACOU_VELOCITY]->GetFeSpace() );
      BiLinFormContext *stiffContVP = new BiLinFormContext(stiffIntVPSurf, STIFFNESS );

      stiffContVP->SetEntities( actSDList, actSDList );
      stiffContVP->SetFeFunctions( velFct, presFct );
      //stiffContVP->SetCounterPart(true);
      assemble_->AddBiLinearForm( stiffContVP );
    }
  }
  
  void PerturbedFlowPDE::DefineSolveStep() {
    solveStep_ = new SolveStepElec(*this);
  }

  void PerturbedFlowPDE::ReadSpecialBCs( ) 
  {
  }
  
   
  // ========================================================================
  // COUPLING SECTION
  // ========================================================================

  void PerturbedFlowPDE::InitCoupling(PDECoupling * Coupling)
  {
    REFACTOR;
  }
  


  void PerturbedFlowPDE::CalcOutputCoupling()
  {
    REFACTOR;
  }


  bool PerturbedFlowPDE::HasOutput(SolutionType output)
  {
    switch (output)
      {
      case FLUIDMECH_VELOCITY:
        return true;
        break;
      case FLUIDMECH_PRESSURE:
        return true;
        break;
      case MEAN_FLUIDMECH_VELOCITY:
        return true;
        break;
      default:
        return false;
        break;
      }
    return false;
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
    pressure->unit = "Pa";

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
    velocity->unit = "m/s";

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
    ef->definedOn = ResultInfo::ELEMENT;
    ef->entryType = ResultInfo::VECTOR;
    availResults_.insert( ef );
    shared_ptr<CoefFunctionFormBased> eFunc;
    if( isComplex_ ) {
      eFunc.reset(new CoefFunctionBOp<Complex>(feFct, ef, -1.0));
    } else {
      eFunc.reset(new CoefFunctionBOp<Double>(feFct, ef, -1.0));
    }

      PtrCoefFct coeffKVP
                = CoefFunction::Generate(Global::REAL, "1.0");
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
      sigmaFunc->AddIntegrator(GetStiffIntegrator( actSDMat, region, isComplex_ ), region);
      strainFunc->AddIntegrator(GetStiffIntegrator( actSDMat, region, isComplex_ ), region);

    }

  }

  
  std::map<SolutionType, shared_ptr<FeSpace> >
  PerturbedFlowPDE::CreateFeSpaces(const std::string& formulation,
		  PtrParamNode infoNode) {
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

    }else{
      EXCEPTION("The formulation " << formulation << "of fluid perturbed PDE is not known!");
    }
    return crSpaces;
  }


  void PerturbedFlowPDE::CreateMeanFlowFunction(StdVector<std::string> dofNames){
    //// === MEAN FLUIDMECH VELOCITY ===
    shared_ptr<ResultInfo> flowvelocity( new ResultInfo);
    flowvelocity->resultType = MEAN_FLUIDMECH_VELOCITY;
    flowvelocity->dofNames = dofNames;
    flowvelocity->unit = "m/s";

    flowvelocity->definedOn = ResultInfo::NODE;
    flowvelocity->entryType = ResultInfo::VECTOR;
    results_.Push_back( flowvelocity );
    availResults_.insert( flowvelocity );
    
    meanFlowCoef_.reset(new CoefFunctionMulti());
    DefineFieldResult( meanFlowCoef_, flowvelocity );
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

    curCoef = regionMat->GetTensorCoefFnc( DYNAMIC_VISCOSITY, subTensorType, complexPart );

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
