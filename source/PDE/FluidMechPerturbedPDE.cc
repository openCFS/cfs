#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <string>
#include <set>

#include "FluidMechPerturbedPDE.hh"

#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Utils/StdVector.hh"
#include "Driver/SolveSteps/SolveStepElec.hh"
#include "CoupledPDE/PDECoupling.hh"
#include "Driver/Assemble.hh"
#include "Utils/ApproxData.hh"
#include "Utils/SmoothSpline.hh"
#include "Materials/Models/Hysteresis.hh"
#include "Materials/Models/Preisach.hh"
#include "FeBasis/H1/FeSpaceH1.hh"
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
#include "Forms/Operators/ConvectiveOperator.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
#include "Domain/Results/ExternalFieldFunctors.hh"

namespace CoupledField {

  DECLARE_LOG(fluidmechpertpde)
  DEFINE_LOG(fluidmechpertpde, "pde.fluidmechpert")

  // ***************
  //   Constructor
  // ***************
  FluidMechPerturbedPDE::FluidMechPerturbedPDE( Grid* grid, PtrParamNode paramNode )
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
  }
  
  void FluidMechPerturbedPDE::InitNonLin() {
    SinglePDE::InitNonLin();
  }

  void FluidMechPerturbedPDE::DefineIntegrators() {
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

    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

    // Get FESpace and FeFunction of fluid velocity
    shared_ptr<BaseFeFunction> velFct = feFunctions_[FLUIDMECH_VELOCITY];
    shared_ptr<FeSpace> velSpace = velFct->GetFeSpace();

    // Get FESpace and FeFunction of fluid pressure
    shared_ptr<BaseFeFunction> presFct = feFunctions_[FLUIDMECH_PRESSURE];
    shared_ptr<FeSpace> presSpace = presFct->GetFeSpace();

    //  Loop over all regions
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {

      // Set current region and material
      actRegion = it->first;
      actMat = it->second;

      // Get current region name
      std::string regionName = ptgrid_->GetRegion().ToString(actRegion);

      // Get ParamNode for current region
      PtrParamNode curRegNode = myParam_->Get("regionList")
          ->GetByVal("region","name",regionName.c_str());

      // create new entity list and add it fefunction
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
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

#if 0
//      BaseBDBInt * stiffInt = NULL;
      shared_ptr<CoefFunction> density =
          actMat->GetCoefFunction( DENSITY, tensorType,
                                   Global::REAL, false);
      shared_ptr<CoefFunction> viscosity =
          actMat->GetCoefFunction( DYNAMIC_VISCOSITY, tensorType,
                                   Global::REAL, false);
#endif

      Double density;
      Double viscosity;
      materials_[actRegion]->GetScalar(density,DENSITY,Global::REAL);
      materials_[actRegion]->GetScalar(viscosity,DYNAMIC_VISCOSITY,Global::REAL);

      // ====================================================================
      // stiffness integrators
      // ====================================================================
      // --------------------------------------------------------------------
      //  VERSION 1: K_PV Integrator (upper off-diagonal integrator)
      // --------------------------------------------------------------------

      shared_ptr<CoefFunction> coeffKPV
                = CoefFunction::Generate(Global::REAL, "1.0");
      BiLinearForm * stiffIntPV = NULL;
      if( dim_ == 2 ) {
        stiffIntPV = new ABInt< GradientOperator<FeH1,2> , IdentityOperator<FeH1,2,2> >(coeffKPV,-density );
      } else {
        stiffIntPV = new ABInt< GradientOperator<FeH1,3> , IdentityOperator<FeH1,3,3> >(coeffKPV,-density );
      }
      stiffIntPV->SetName("PerturbedStiffIntPV");
      //stiffIntPV->SetFeSpace( feFunctions_[ACOU_PRESSURE]->GetFeSpace(), feFunctions_[ACOU_VELOCITY]->GetFeSpace() );
      BiLinFormContext *stiffContPV = new BiLinFormContext(stiffIntPV, STIFFNESS );

      stiffContPV->SetEntities( actSDList, actSDList );
      stiffContPV->SetFeFunctions( presFct, velFct);
//      stiffContPV->SetCounterPart(true);
      assemble_->AddBiLinearForm( stiffContPV );

      // --------------------------------------------------------------------
      //  VERSION 2: K_VP = K_PV^T Integrator (lower off-diagonal integrator)
      // --------------------------------------------------------------------
      shared_ptr<CoefFunction> coeffKVP
                = CoefFunction::Generate(Global::REAL, "1.0");
      BiLinearForm * stiffIntVP = NULL;
      if( dim_ == 2 ) {
        stiffIntVP = new ABInt< IdentityOperator<FeH1,2,2> , GradientOperator<FeH1,2> >(coeffKVP,1.0 );
      } else {
        stiffIntVP = new ABInt< IdentityOperator<FeH1,3,3> , GradientOperator<FeH1,3> >(coeffKVP,1.0 );
      }
      stiffIntVP->SetName("PerturbedStiffIntVP");
      //stiffIntVP->SetFeSpace( feFunctions_[ACOU_PRESSURE]->GetFeSpace(), feFunctions_[ACOU_VELOCITY]->GetFeSpace() );
      BiLinFormContext *stiffContVP = new BiLinFormContext(stiffIntVP, STIFFNESS );

      stiffContVP->SetEntities( actSDList, actSDList );
      stiffContVP->SetFeFunctions( velFct, presFct);
      //stiffContVP->SetCounterPart(true);
      assemble_->AddBiLinearForm( stiffContVP );

      // --------------------------------------------------------------------
      //  VERSION 2: K_Laplace Integrator
      // --------------------------------------------------------------------
      shared_ptr<CoefFunction> coeffKvv
                = CoefFunction::Generate(Global::REAL, "1.0");
      BiLinearForm * stiffIntLaplace = NULL;
      if( dim_ == 2 ) {
        stiffIntLaplace = new BBInt< DivOperator<FeH1,2> >(coeffKvv, viscosity);
      } else {
        stiffIntLaplace = new BBInt< DivOperator<FeH1,3> >(coeffKvv, viscosity);
      }
      stiffIntLaplace->SetName("PerturbedStiffIntViscous");
      BiLinFormContext *stiffContLaplace = new BiLinFormContext(stiffIntLaplace, STIFFNESS );

      stiffContLaplace->SetEntities( actSDList, actSDList );
      stiffContLaplace->SetFeFunctions( feFunctions_[FLUIDMECH_VELOCITY],feFunctions_[FLUIDMECH_VELOCITY]);
      assemble_->AddBiLinearForm( stiffContLaplace );


      shared_ptr<CoefFunction> coeffMVV
                = CoefFunction::Generate(Global::REAL, lexical_cast<std::string>(1.0));

      //======================================================================
      // CHECK FOR FLOW
      //=====================================================================
      std::string flowId = curRegNode->Get("flowId")->As<std::string>();
      if(flowId != ""){
        //Add the region information
        PtrParamNode flowNode = myParam_->Get("flowList")->GetByVal("flow","name",flowId.c_str());
        if(isComplex_){
          shared_ptr<FieldFunctor<Complex> > fct = dynamic_pointer_cast<FieldFunctor<Complex> >(meanFlowFunctor_);
          fct->AddRegion(actRegion,flowNode);
        }else{
          shared_ptr<FieldFunctor<Double> > fct = dynamic_pointer_cast<FieldFunctor<Double> >(meanFlowFunctor_);
          fct->AddRegion(actRegion,flowNode);
        }

        //now create the integrators
        BiLinearForm *convectiveVv = NULL;
        if( dim_ == 2 ) {
          convectiveVv = new ABInt<IdentityOperator<FeH1,2,2>,ConvectiveOperator<FeH1,2,2> >(coeffMVV, 1.0);
        } else {
          convectiveVv = new ABInt<IdentityOperator<FeH1,3,3>,ConvectiveOperator<FeH1,3,3>  >(coeffMVV, 1.0);
        }

        convectiveVv->SetBCoefFunctionOpB(meanFlowFunctor_);

        convectiveVv->SetName("PerturbedStiffIntConvective");

        BiLinFormContext *convectiveContextVv =  new BiLinFormContext(convectiveVv, STIFFNESS );

        convectiveContextVv->SetEntities( actSDList, actSDList );
        convectiveContextVv->SetFeFunctions( feFunctions_[FLUIDMECH_VELOCITY],feFunctions_[FLUIDMECH_VELOCITY]);
        assemble_->AddBiLinearForm( convectiveContextVv );
      }

      // ====================================================================
      // damping integrators
      // ====================================================================
      shared_ptr<CoefFunction> coeffDvv
                = CoefFunction::Generate(Global::REAL, lexical_cast<std::string>(1.0));

      BiLinearForm *dampIntvv = NULL;
      if( dim_ == 2 ) {
        dampIntvv = new BBInt<IdentityOperator<FeH1,2,2> >(coeffDvv, 1.0 / density );
      } else {
        dampIntvv = new BBInt<IdentityOperator<FeH1,3,3> >(coeffDvv, 1.0 / density );
      }
      dampIntvv->SetName("PerturbedDampInt");
      //massIntPP->SetFeSpace( feFunctions_[ACOU_PRESSURE]->GetFeSpace() );

      BiLinFormContext *dampContextvv =  new BiLinFormContext(dampIntvv, DAMPING );

      dampContextvv->SetEntities( actSDList, actSDList );
      dampContextvv->SetFeFunctions( feFunctions_[FLUIDMECH_VELOCITY],feFunctions_[FLUIDMECH_VELOCITY]);
      assemble_->AddBiLinearForm( dampContextvv );

      
    }
  }
  
  void FluidMechPerturbedPDE::DefineSolveStep() {
    solveStep_ = new SolveStepElec(*this);
  }

  void FluidMechPerturbedPDE::ReadSpecialBCs( ) 
  {
     //ReadImpedances();
  }
  
  // ========================================================================
  // POSTPROCESSING SECTION
  // ========================================================================

  // ************************************************************************
  //   CalcResults
  // ************************************************************************
  void FluidMechPerturbedPDE::CalcResults( shared_ptr<BaseResult> res ) {

    switch (res->GetResultInfo()->resultType ) {
      case FLUIDMECH_VELOCITY:
        feFunctions_[FLUIDMECH_VELOCITY]->ExtractResult( res );
        break;

      case FLUIDMECH_PRESSURE:
        feFunctions_[FLUIDMECH_PRESSURE]->ExtractResult( res );
        break;

      case MEAN_FLUIDMECH_VELOCITY:
        meanFlowFunctor_->EvalResult( res );
        break;

      default:
        WARN( "Result '" << 
              SolutionTypeEnum.ToString(res->GetResultInfo()->resultType)
              << "' type not computable by electric PDE" );
        break;
    }

  }
   
  // ========================================================================
  // COUPLING SECTION
  // ========================================================================

  void FluidMechPerturbedPDE::InitCoupling(PDECoupling * Coupling)
  {
  REFACTOR;
  }
  


  void FluidMechPerturbedPDE::CalcOutputCoupling()
  {
    REFACTOR;
  }


  bool FluidMechPerturbedPDE::HasOutput(SolutionType output)
  {
    switch (output)
      {
      case FLUIDMECH_VELOCITY:
        return true;
        break;
      case FLUIDMECH_PRESSURE:
        return true;
        break;
      default:
        return false;
        break;
      }
    return false;
  }

  void FluidMechPerturbedPDE::DefinePrimaryResults() {
    shared_ptr<BaseFeFunction> feFct = feFunctions_[FLUIDMECH_VELOCITY];

    // Check for subType
    StdVector<std::string> velDofNames;

    std::string geometryType;
    param->Get("domain")->GetValue("geometryType", geometryType );

    if( geometryType == "3d" ) {
      velDofNames = "x", "y", "z";
    } else if( geometryType == "plane" ) {
      velDofNames = "x", "y";
    } else if( geometryType == "axi" ) {
      velDofNames = "r", "z";
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

    // MEAN VELOCITY
    CreateMeanFlowFunction(velDofNames);
  }
  
  
  void FluidMechPerturbedPDE::DefinePostProcResults() {

  }

  
  std::map<SolutionType, shared_ptr<FeSpace> >
  FluidMechPerturbedPDE::CreateFeSpaces(const std::string& formulation,
		  PtrParamNode infoNode) {
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;

    std::cout << "Creating FESpaces, formulation " << formulation << std::endl;

    if(formulation == "default" || formulation == "H1"){
      PtrParamNode spaceNode;
      spaceNode = infoNode->Get(SolutionTypeEnum.ToString(FLUIDMECH_VELOCITY));

      crSpaces[FLUIDMECH_VELOCITY] =
        FeSpace::CreateInstance(myParam_,spaceNode,FeSpace::H1);
      crSpaces[FLUIDMECH_VELOCITY]->Init();

      spaceNode = infoNode->Get(SolutionTypeEnum.ToString(FLUIDMECH_PRESSURE));

      crSpaces[FLUIDMECH_PRESSURE] =
        FeSpace::CreateInstance(myParam_,spaceNode,FeSpace::H1);
      crSpaces[FLUIDMECH_PRESSURE]->Init();

    }else{
      EXCEPTION("The formulation " << formulation << "of fluid perturbed PDE is not known!");
    }
    return crSpaces;
  }

  LinearFormContext* FluidMechPerturbedPDE::CreateRhsLinearForm(SolutionType rhsType,
	  shared_ptr<CoefFunction > rhsCoef){
#if 0
    LinearFormContext * mContext = NULL;
    switch(rhsType){
    case ELEC_CHARGE_DENSITY:
      BUIntegrator<IdentityOperator,FeH1,Double>* curInt;
      curInt = new BUIntegrator<IdentityOperator,FeH1,Double>(1.0,rhsCoef);

      mContext = new LinearFormContext(curInt);
      mContext->SetFeFunction( feFunctions_[FLUIDMECH_VELOCITY]);
      curInt->SetFeSpace( feFunctions_[FLUIDMECH_VELOCITY]->GetFeSpace());
      break;
    default:
      Exception("Right hand side quantity not known for elecPDE");
      break;
    }
    return mContext;
#endif
    return NULL;
  }

  void FluidMechPerturbedPDE::CreateMeanFlowFunction(StdVector<std::string> dofNames){
    //// === MEAN FLUIDMECH VELOCITY ===
    shared_ptr<ResultInfo> flowvelocity( new ResultInfo);
    flowvelocity->resultType = MEAN_FLUIDMECH_VELOCITY;
    flowvelocity->dofNames = dofNames;
    flowvelocity->unit = "m/s";

    flowvelocity->definedOn = ResultInfo::NODE;
    flowvelocity->entryType = ResultInfo::VECTOR;



    shared_ptr<BaseFeFunction> meanFunction;
    std::string form = SolutionTypeEnum.ToString(MEAN_FLUIDMECH_VELOCITY);
    PtrParamNode feSpaceNode = infoNode_->Get("feSpaces");
    PtrParamNode potSpaceNode = feSpaceNode->Get(form);
    shared_ptr<FeSpace> tmpSpace = FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1);

    if(isComplex_){
      meanFunction.reset(new FeFunction<Complex>());
      meanFunction->SetFeSpace(tmpSpace);
      meanFunction->SetResultInfo(flowvelocity);
      meanFunction->SetGrid(ptgrid_);
      meanFunction->SetPDE(this);
      flowvelocity->SetFeFunction(meanFunction);
      if(dim_==2)
        meanFlowFunctor_.reset(new ExternalFieldFunctor<IdentityOperator<FeH1,2,2,Complex>,Complex >(meanFunction,flowvelocity));
      else
        meanFlowFunctor_.reset(new ExternalFieldFunctor<IdentityOperator<FeH1,3,3,Complex>,Complex >(meanFunction,flowvelocity));
    }else{
      meanFunction.reset(new FeFunction<Double>());
      meanFunction->SetFeSpace(tmpSpace);
      meanFunction->SetResultInfo(flowvelocity);
      meanFunction->SetGrid(ptgrid_);
      meanFunction->SetPDE(this);
      flowvelocity->SetFeFunction(meanFunction);
      if(dim_==2)
        meanFlowFunctor_.reset(new ExternalFieldFunctor<IdentityOperator<FeH1,2,2,Double>,Double >(meanFunction,flowvelocity));
      else
        meanFlowFunctor_.reset(new ExternalFieldFunctor<IdentityOperator<FeH1,3,3,Double>,Double >(meanFunction,flowvelocity));
    }

    results_.Push_back( flowvelocity );
    availResults_.insert( flowvelocity );

  }

}
