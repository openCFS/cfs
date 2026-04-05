// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "FluidMechCoupling.hh"

#include "PDE/SinglePDE.hh"
#include "PDE/PerturbedFlowPDE.hh"
#include "PDE/MechPDE.hh"
#include "CoupledPDE/BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "Materials/BaseMaterial.hh"
#include "Driver/FormsContexts.hh"
#include "Driver/Assemble.hh"

// include fespaces
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"


namespace CoupledField {



  // ***************
  //   Constructor
  // ***************
  FluidMechCoupling::FluidMechCoupling( SinglePDE *pde1, SinglePDE *pde2,
                                PtrParamNode paramNode, 
                                PtrParamNode infoNode,
                                shared_ptr<SimState> simState,
                                Domain* domain)
    : BasePairCoupling( pde1, pde2, paramNode, infoNode, simState, domain ),
      lmOrderSameAsVel_(true)
  {
    couplingName_ = "FluidMechDirect";
    materialClass_ = FLOW;

    formulation_ = NO_SOLUTION_TYPE;

    // determine subtype from mechanic pde
    pde1_->GetParamNode()->GetValue( "subType", subType_ );

    nonLin_ = false;
    
    // Initialize nonlinearities
    InitNonLin();

    if(paramNode->Has("lmOrderSameAsVel")) 
    {
      lmOrderSameAsVel_ =  paramNode->Get("lmOrderSameAsVel")->As<bool>();
    }    
  }


  // **************
  //   Destructor
  // **************
  FluidMechCoupling::~FluidMechCoupling() {
  }


  // *********************
  //   DefineIntegrators
  // *********************
  void FluidMechCoupling::DefineIntegrators() {

    // get hold of both feFunctions
    MechPDE* mechPDE = dynamic_cast<MechPDE*>(pde1_);
    shared_ptr<BaseFeFunction> dispFct = mechPDE->GetFeFunction(MECH_DISPLACEMENT);
//    std::map<RegionIdType, BaseMaterial*> mechMaterials;
//    mechMaterials = mechPDE->getPDEMaterialData();
    MathParser * mp = domain_->GetMathParser();

    PerturbedFlowPDE* flowPDE = dynamic_cast<PerturbedFlowPDE*>(pde2_);
    shared_ptr<BaseFeFunction> velFct = flowPDE->GetFeFunction(FLUIDMECH_VELOCITY);
    shared_ptr<BaseFeFunction> presFct = flowPDE->GetFeFunction(FLUIDMECH_PRESSURE);

    shared_ptr<BaseFeFunction> lagrangeMultFct = feFunctions_[LAGRANGE_MULT];

    std::map<RegionIdType, BaseMaterial*> flowMaterials;
    flowMaterials = flowPDE->GetMaterialData();

    // Create coefficient functions for all fluid densities
    std::map< RegionIdType, PtrCoefFct > oneFuncs;
    std::set< RegionIdType > flowRegions;
    std::map<RegionIdType, BaseMaterial*>::iterator it, end;
    it = flowMaterials.begin();
    end = flowMaterials.end();

    shared_ptr<FeSpace> dispSpace = dispFct->GetFeSpace();
    shared_ptr<FeSpace> velSpace = velFct->GetFeSpace();
    shared_ptr<FeSpace> presSpace = presFct->GetFeSpace();
    shared_ptr<FeSpace> lagrangeMultSpace = lagrangeMultFct->GetFeSpace();

    for( ; it != end; it++ ) {
      RegionIdType volRegId = it->first;

      flowRegions.insert(volRegId);

      // Get bulk density for acoustics
      BaseMaterial * flowMat = flowMaterials[volRegId];
//      Double density = 1.0;
//      Double viscosity = 1.0;
//      flowMat->GetScalar(density,DENSITY,Global::REAL);
//      flowMat->GetScalar(viscosity,DYNAMIC_VISCOSITY,Global::REAL);

      PtrCoefFct density = flowMat->GetScalCoefFnc( DENSITY, Global::REAL );
      PtrCoefFct viscosity = flowMat->GetScalCoefFnc( FLUID_DYNAMIC_VISCOSITY, Global::REAL );
      
      oneFuncs[volRegId] = CoefFunction::Generate(mp, Global::REAL, "1.0");

      WARN("fluid density: " << density->ToString() << " dynamic viscosity " << viscosity->ToString());

      shared_ptr<ElemList> actSDList( new ElemList( ptGrid_ ) );
      actSDList->SetRegion( volRegId );

      //      if(actSDList->GetType() != EntityList::SURF_ELEM_LIST) {
      //        lagrangeMultFct->b( actSDList );
      //        lagrangeMultSpace->SetRegionApproximation(volRegId, "velPolyId", "velIntegId");
      //      }
    }

    for ( UInt actSD = 0, n = entityLists_.GetSize(); actSD < n; actSD++ ) {

      shared_ptr<SurfElemList> actSDList = 
          dynamic_pointer_cast<SurfElemList>(entityLists_[actSD]);
      RegionIdType region = actSDList->GetRegion();

      velFct->AddEntityList(actSDList);
      presFct->AddEntityList(actSDList);
      dispFct->AddEntityList(actSDList);
      lagrangeMultFct->AddEntityList(actSDList);

      velSpace->SetRegionApproximation(region, "velPolyId", "velIntegId");
      presSpace->SetRegionApproximation(region, "presPolyId", "presIntegId");
      dispSpace->SetRegionApproximation(region, "default", "default");
      if(lmOrderSameAsVel_)
      {
        lagrangeMultSpace->SetRegionApproximation(region, "velPolyId", "velIntegId");
      }
      else
      {
        lagrangeMultSpace->SetRegionApproximation(region, "presPolyId", "velIntegId");
      }

      // This integrator gets assembled into the damping (first time deriv.) matrix in the row of the LM
      DefineDampingIntegrators("FluidMechDampingLMVelCouplingInt",
                            dispFct,
                            lagrangeMultFct,
                            actSDList,
                            oneFuncs,
                            flowRegions);

      // This integrator gets assembled into the stiffness matrix of the mechanic PDE
      DefineStiffnessIntegrators("FluidMechStiff",
                                 dispFct,
                                 velFct,
                                 lagrangeMultFct,
                                 actSDList,
                                 oneFuncs,
                                 oneFuncs,
                                 oneFuncs,
                                 flowRegions);
    }
  }

  void FluidMechCoupling::DefineDampingIntegrators(const std::string& name,
                                                shared_ptr<BaseFeFunction>& dispFct,
                                                shared_ptr<BaseFeFunction>& lmFct,
                                                shared_ptr<SurfElemList>& actSDList,
                                                const std::map< RegionIdType, PtrCoefFct >& oneFuncs,
                                                const std::set< RegionIdType >& flowRegions){
    BiLinearForm * dampInt = NULL;

    if( subType_ == "axi" ) {
      dampInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                   new IdentityOperator<FeH1,2,2>(),
                                   oneFuncs, 1.0, flowRegions);
    } else if( subType_ == "planeStrain" ) {
      dampInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                   new IdentityOperator<FeH1,2,2>(),
                                   oneFuncs, 1.0, flowRegions);
    } else if( subType_ == "planeStress" ) {
      dampInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                   new IdentityOperator<FeH1,2,2>(),
                                   oneFuncs, 1.0, flowRegions);
    } else if( subType_ == "3d") {
      dampInt = new SurfaceABInt<>(new IdentityOperator<FeH1,3,3>(),
                                   new IdentityOperator<FeH1,3,3>(),
                                   oneFuncs, 1.0, flowRegions);
    } else {
      EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
    }

    dampInt->SetName(name);
    BiLinFormContext * context =
        new BiLinFormContext(dampInt, DAMPING );
    
    context->SetEntities( actSDList, actSDList );
    context->SetFeFunctions( lmFct, dispFct );
    context->SetCounterPart(false);

    assemble_->AddBiLinearForm( context );
  }

  void FluidMechCoupling::DefineStiffnessIntegrators(const std::string& name,
                                                shared_ptr<BaseFeFunction>& dispFct,
                                                shared_ptr<BaseFeFunction>& velFct,
                                                shared_ptr<BaseFeFunction>& lmFct,
                                                shared_ptr<SurfElemList>& actSDList,
                                                const std::map< RegionIdType, PtrCoefFct >& densityFuncs,
                                                const std::map< RegionIdType, PtrCoefFct >& muFuncs,
                                                const std::map< RegionIdType, PtrCoefFct >& oneFuncs,
                                                const std::set< RegionIdType >& flowRegions){
    BiLinearForm * stiffInt = NULL;

    // LM-velocity integrator in row of LM and column of velocity
    std::string intName  = name + "LMVelCouplingInt";
    if( subType_ == "axi" ) {
        stiffInt = new SurfaceABInt<>( new IdentityOperator<FeH1,2,2>(),
                                       new IdentityOperator<FeH1,2,2>(),
                                       oneFuncs, -1.0, flowRegions);
    } else if( subType_ == "planeStrain" ) {
      stiffInt = new SurfaceABInt<>(new  IdentityOperator<FeH1,2,2>(),
                                    new IdentityOperator<FeH1,2,2>(),
                                    oneFuncs, -1.0, flowRegions);
    } else if( subType_ == "planeStress" ) {
      stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                    new IdentityOperator<FeH1,2,2>(),
                                    oneFuncs, -1.0, flowRegions);
    } else if( subType_ == "3d") {
      stiffInt = new SurfaceABInt<>( new IdentityOperator<FeH1,3,3>(),
                                     new IdentityOperator<FeH1,3,3>(),
                                     oneFuncs, -1.0, flowRegions);
    } else {
      EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
    }

    stiffInt->SetName(intName);
    BiLinFormContext * context =
        new BiLinFormContext( stiffInt, STIFFNESS );

    context->SetEntities( actSDList, actSDList );
    context->SetFeFunctions( lmFct, velFct );
    context->SetCounterPart(false);

    assemble_->AddBiLinearForm( context );

    // Displacement-LM integrator in row of displacement and column of LM
    intName  = name + "DispLMCouplingInt";
    if( subType_ == "axi" ) {
        stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                      new IdentityOperator<FeH1,2,2>(),
                                      oneFuncs, 1.0, flowRegions);
    } else if( subType_ == "planeStrain" ) {
      stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                    new IdentityOperator<FeH1,2,2>(),
                                    oneFuncs, 1.0, flowRegions);
    } else if( subType_ == "planeStress" ) {
      stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                    new IdentityOperator<FeH1,2,2>(),
                                    oneFuncs, 1.0, flowRegions);
    } else if( subType_ == "3d") {
      stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,3,3>(),
                                    new IdentityOperator<FeH1,3,3>(),
                                    oneFuncs, 1.0, flowRegions);
    } else {
      EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
    }

    stiffInt->SetName(intName);
    context = new BiLinFormContext( stiffInt, STIFFNESS );

    context->SetEntities( actSDList, actSDList );
    context->SetFeFunctions( dispFct, lmFct );
    context->SetCounterPart(false);

    assemble_->AddBiLinearForm( context );

    // Velocity-LM integrator in row of velocity and column of LM
    intName  = name + "VelLMCouplingInt";
    if( subType_ == "axi" ) {
        stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                      new IdentityOperator<FeH1,2,2>(),
                                      oneFuncs, -1.0, flowRegions);
        //          (densityFuncs, 1.0, flowRegions);
    } else if( subType_ == "planeStrain" ) {
      stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                    new IdentityOperator<FeH1,2,2>(),
                                    oneFuncs, -1.0, flowRegions);
                 //        (densityFuncs, 1.0, flowRegions);
    } else if( subType_ == "planeStress" ) {
      stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
                                    new IdentityOperator<FeH1,2,2>(),
                                    oneFuncs, -1.0, flowRegions );
                 //        (densityFuncs, 1.0, flowRegions);
    } else if( subType_ == "3d") {
      stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,3,3>(),
                                    new IdentityOperator<FeH1,3,3>(),
                                    oneFuncs, -1.0, flowRegions );
                 //        (densityFuncs, 1.0, flowRegions);
    } else {
      EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
    }

    stiffInt->SetName(intName);
    context = new BiLinFormContext( stiffInt, STIFFNESS );

    context->SetEntities( actSDList, actSDList );
    context->SetFeFunctions( velFct, lmFct );
    context->SetCounterPart(false);

    assemble_->AddBiLinearForm( context );

  }

  void FluidMechCoupling::DefineAvailResults() {
    REFACTOR  
  }

  void FluidMechCoupling::DefinePrimaryResults() {
    // Check for subType
    StdVector<std::string> velDofNames;

    std::string geometryType;
    domain_->GetParamRoot()->Get("domain")->GetValue("geometryType", geometryType );

    if( geometryType == "3d" ) {
      velDofNames = "x", "y", "z";
    } else if( geometryType == "plane" ) {
      velDofNames = "x", "y";
    } else if( geometryType == "axi" ) {
      velDofNames = "r", "z";
    }

    // === LAGRANGE MULTIPLIER ===
    shared_ptr<ResultInfo> res1( new ResultInfo);
    res1->resultType = LAGRANGE_MULT;
    
    res1->dofNames = velDofNames;
    res1->unit = "Pa";
    res1->definedOn = ResultInfo::NODE;
    res1->entryType = ResultInfo::VECTOR;
    feFunctions_[LAGRANGE_MULT]->SetResultInfo(res1);
    results_.Push_back( res1 );
    availResults_.insert( res1 );
    
    res1->SetFeFunction(feFunctions_[LAGRANGE_MULT]);
    
    DefineFieldResult( feFunctions_[LAGRANGE_MULT], res1 );
  }


 void FluidMechCoupling::CreateFeSpaces( const std::string&  type,
                                         PtrParamNode infoNode,
                                         std::map<SolutionType, shared_ptr<FeSpace> >& crSpaces) {

   //we need a lagrange multiplier
   formulation_ = LAGRANGE_MULT;

   PtrParamNode spaceNode;
   if(lmOrderSameAsVel_)
   {
     spaceNode = infoNode->Get(SolutionTypeEnum.ToString(FLUIDMECH_VELOCITY));
   }
   else
   {
     spaceNode = infoNode->Get(SolutionTypeEnum.ToString(FLUIDMECH_PRESSURE));
   }
   
   
   crSpaces[formulation_] =
       FeSpace::CreateInstance(myParam_, spaceNode, FeSpace::H1, ptGrid_);

   crSpaces[formulation_]->SetLagrSurfSpace();

   crSpaces[formulation_]->Init(solStrat_);
  }
  
}

