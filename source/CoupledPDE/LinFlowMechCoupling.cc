// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "LinFlowMechCoupling.hh"

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
  LinFlowMechCoupling::LinFlowMechCoupling( SinglePDE *pde1, SinglePDE *pde2,
                                PtrParamNode paramNode, 
                                PtrParamNode infoNode,
                                shared_ptr<SimState> simState,
                                Domain* domain)
    : BasePairCoupling( pde1, pde2, paramNode, infoNode, simState, domain ),
      lmOrderSameAsVel_(true) {

    couplingName_ = "linFlowMechDirect";
    materialClass_ = FLOW;

    formulation_ = NO_SOLUTION_TYPE;

    // determine subtype from mechanic pde
    pde2_->GetParamNode()->GetValue( "subType", subType_ );

    nonLin_ = false;
    
    // Initialize nonlinearities
    InitNonLin();

    if( paramNode->Has("lmOrderSameAsVel") ) {
      lmOrderSameAsVel_ =  paramNode->Get("lmOrderSameAsVel")->As<bool>();
    }
  }


  // **************
  //   Destructor
  // **************
  LinFlowMechCoupling::~LinFlowMechCoupling() {
  }


  // *********************
  //   DefineIntegrators
  // *********************
  void LinFlowMechCoupling::DefineIntegrators() {

    MathParser * mp = domain_->GetMathParser();

    // fixed in Domain::CreateDirectCoupledPDEs: pde1 is fluidMechLin and
    // pde2 is mechanic
    //
    shared_ptr<BaseFeFunction> velFct  = pde1_->GetFeFunction(FLUIDMECH_VELOCITY);
    shared_ptr<BaseFeFunction> presFct = pde1_->GetFeFunction(FLUIDMECH_PRESSURE);
    shared_ptr<BaseFeFunction> dispFct = pde2_->GetFeFunction(MECH_DISPLACEMENT);

    shared_ptr<BaseFeFunction> lagrangeMultFct = feFunctions_[LAGRANGE_MULT];

    // Create coefficient functions
    std::map< RegionIdType, PtrCoefFct > oneFuncs;

    shared_ptr<FeSpace> dispSpace = dispFct->GetFeSpace();
    shared_ptr<FeSpace> velSpace = velFct->GetFeSpace();
    shared_ptr<FeSpace> presSpace = presFct->GetFeSpace();
    shared_ptr<FeSpace> lagrangeMultSpace = lagrangeMultFct->GetFeSpace();

    std::set< RegionIdType > flowRegions;
    std::map<RegionIdType, BaseMaterial*> flowMaterials;
    flowMaterials = pde1_->GetMaterialData();
    std::map<RegionIdType, BaseMaterial*>::iterator it, end;
    it = flowMaterials.begin();
    end = flowMaterials.end();

    for ( ; it != end; it++ ) {
      RegionIdType volRegId = it->first;
      flowRegions.insert(volRegId);

      oneFuncs[volRegId] = CoefFunction::Generate(mp, Global::REAL,
                                                   lexical_cast<std::string>(1.0));
      shared_ptr<ElemList> actSDList( new ElemList( ptGrid_ ) );
      actSDList->SetRegion( volRegId );
    }

    //go over all coupled surfaces
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
      if(lmOrderSameAsVel_) {
        lagrangeMultSpace->SetRegionApproximation(region, "velPolyId", "velIntegId");
      }
      else {
        lagrangeMultSpace->SetRegionApproximation(region, "presPolyId", "velIntegId");
      }

      // Integrator being assembled into damping (first time deriv.) matrix; first part
      // of additional equation guaranteeing continuity of velocities
      DefineDampingIntegrators("LinFlowMechDampingLMVelCouplingInt", dispFct, lagrangeMultFct,
                                actSDList, oneFuncs, flowRegions);

      // These integrators gets assembled into the stiffness matrix of the mechanic PDE,
      // LinFlowPDE and the additional equation guaranteeing continuity of velocities
      // equation for continuity of velocities)
      DefineStiffnessIntegrators("LinFlowMechStiff", dispFct, velFct, lagrangeMultFct,
                                  actSDList, oneFuncs, oneFuncs, oneFuncs, flowRegions);
    }
  }

  void LinFlowMechCoupling::DefineDampingIntegrators(const std::string& name,
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

  void LinFlowMechCoupling::DefineStiffnessIntegrators(const std::string& name,
                                                shared_ptr<BaseFeFunction>& dispFct,
                                                shared_ptr<BaseFeFunction>& velFct,
                                                shared_ptr<BaseFeFunction>& lmFct,
                                                shared_ptr<SurfElemList>& actSDList,
                                                const std::map< RegionIdType, PtrCoefFct >& densityFuncs,
                                                const std::map< RegionIdType, PtrCoefFct >& muFuncs,
                                                const std::map< RegionIdType, PtrCoefFct >& oneFuncs,
                                                const std::set< RegionIdType >& flowRegions){
    BiLinearForm * stiffInt = NULL;

    // LM-velocity integrator in row of LM and column of lin. flow velocity
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

    // Velocity-LM integrator in row of flow velocity and column of LM
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

  void LinFlowMechCoupling::DefineAvailResults() {
    REFACTOR  
  }

  void LinFlowMechCoupling::DefinePrimaryResults() {
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


 void LinFlowMechCoupling::CreateFeSpaces( const std::string&  type,
                                         PtrParamNode infoNode,
                                         std::map<SolutionType, shared_ptr<FeSpace> >& crSpaces) {

   //we need a lagrange multiplier
   formulation_ = LAGRANGE_MULT;

   PtrParamNode spaceNode;
   if(lmOrderSameAsVel_) {
     spaceNode = infoNode->Get(SolutionTypeEnum.ToString(FLUIDMECH_VELOCITY));
   }
   else {
     spaceNode = infoNode->Get(SolutionTypeEnum.ToString(FLUIDMECH_PRESSURE));
   }

   crSpaces[formulation_] =
       FeSpace::CreateInstance(myParam_, spaceNode, FeSpace::H1, ptGrid_);

   crSpaces[formulation_]->SetLagrSurfSpace();

   crSpaces[formulation_]->Init(solStrat_);
  }
  
}

