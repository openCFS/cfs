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
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"


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

    //now do PDE specifics
    std::map<std::string, NonLinType>::iterator it;
    for ( it=nonLinTypes_.begin() ; it != nonLinTypes_.end(); it++ ) {
      if ( (*it).second == HYSTERESIS ) {
        isHysteresis_ = true;
      }
    }

  }

  void FluidMechPerturbedPDE::DefineIntegrators() {
#if 0
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;
    
    // flag for updatedLagrange formulation
    bool upLagrangeForm = true;
    upLagrangeForm = true;
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

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    shared_ptr<FeSpace> mySpace = feFunctions_[FLUIDMECH_VELOCITY]->GetFeSpace();
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      
      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;

      // Get current region name
      std::string regionName = ptgrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );
      
      // ====================================================================
      // New implementation
      // ====================================================================
      // --- Set the approximation for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId, integId);

      // --- standard real-valued stiffness integrator ---
      shared_ptr<CoefFunction> curCoef = actSDMat->GetCoefFunction(DENSITY,tensorType,Global::REAL);


      BDBInt< FeH1, Double, Double >* stiffInt;
      stiffInt = new BDBInt< FeH1, Double, Double >(curCoef,1.0 );
      //stiffInt->SetName("StiffnessInt");
      //linElecInt *  linElecForm = new linElecInt( actSDMat, tensorType,
       //                                           upLagrangeForm );
      //linElecForm->SetFactor( factor );
      BiLinFormContext * stiffIntDescr =
          new BiLinFormContext(stiffInt, STIFFNESS );

      feFunctions_[FLUIDMECH_VELOCITY]->AddEntityList( actSDList );

      //stiffIntDescr->SetPtPdes(this, this);
      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions( feFunctions_[FLUIDMECH_VELOCITY],feFunctions_[FLUIDMECH_VELOCITY]);
      stiffInt->SetFeSpace( feFunctions_[FLUIDMECH_VELOCITY]->GetFeSpace());

      assemble_->AddBiLinearForm( stiffIntDescr );
      // Important: Add bdb-integrator to global list, as we need them later
      // for calculation of postprocessing results
      bdbInts_[actRegion] = stiffInt;
      
    }
#endif
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

      case ELEC_RHS_LOAD:
        rhsFeFunctions_[FLUIDMECH_VELOCITY]->ExtractResult( res );
        break;

      case ELEC_FIELD_INTENSITY:
        resultFunctors_[ELEC_FIELD_INTENSITY]->EvalResult(res);
        break;
        
      case ELEC_FLUX_DENSITY:
        resultFunctors_[ELEC_FLUX_DENSITY]->EvalResult(res);
        break;

      case ELEC_ENERGY_DENSITY:
        resultFunctors_[ELEC_ENERGY_DENSITY]->EvalResult(res);
        break;
        
      case ELEC_ENERGY:
        resultFunctors_[ELEC_ENERGY]->EvalResult(res);
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
      case ELEC_FORCE_VWP:
        return true;
        break;
      case FLUIDMECH_VELOCITY:
        return true;
        break;
      case ELEC_FIELD_INTENSITY:
        return true;
        break;
      case ELEC_INTERFACE_FORCE:
        return true;
        break;
      default:
        return false;
        break;
      }
    return false;
  }

  void FluidMechPerturbedPDE::DefinePrimaryResults() {
#if 0
    shared_ptr<BaseFeFunction> feFct = feFunctions_[FLUIDMECH_VELOCITY];
    
    // Velocity field
    shared_ptr<ResultInfo> res1( new ResultInfo);
    res1->resultType = FLUIDMECH_VELOCITY;
    res1->dofNames = "x", "y";
    res1->unit = MapSolTypeToUnit(FLUIDMECH_VELOCITY);
    res1->definedOn = ResultInfo::NODE;
    res1->entryType = ResultInfo::VECTOR;
    feFunctions_[FLUIDMECH_VELOCITY]->SetResultInfo(res1);

    res1->SetFeFunction(feFunctions_[FLUIDMECH_VELOCITY]);
    results_.Push_back( res1 );
    availResults_.insert( res1 );
    shared_ptr<BaseFieldFunctor> vFunc;
    if( isComplex_ ) {
      vFunc.reset(
          new FieldInterpolFunctor<IdentityOperator,
          FeH1,
          Complex>(feFct, res1));
    } else {
      vFunc.reset(
          new FieldInterpolFunctor<IdentityOperator,
          FeH1,
          Double>(feFct, res1));
    }
    resultFunctors_[FLUIDMECH_VELOCITY] = vFunc;
    fieldFunctors_[FLUIDMECH_VELOCITY] = vFunc;
    
    

    // Pressure field
    shared_ptr<ResultInfo> res2 ( new ResultInfo );
    res2->resultType = FLUIDMECH_PRESSURE;
    res2->dofNames = "";
    res2->unit = MapSolTypeToUnit(FLUIDMECH_PRESSURE);
    res2->definedOn = results_[0]->definedOn;
    res2->entryType = ResultInfo::SCALAR;
    feFunctions_[FLUIDMECH_PRESSURE]->SetResultInfo(res1);
    feFct = feFunctions_[FLUIDMECH_PRESSURE];

    res2->SetFeFunction(feFunctions_[FLUIDMECH_PRESSURE]);
    results_.Push_back( res2 );
    availResults_.insert( res2 );
    if( isComplex_ ) {
      vFunc.reset(
          new FieldInterpolFunctor<IdentityOperator,
          FeH1,
          Complex>(feFct, res2));
    } else {
      vFunc.reset(
          new FieldInterpolFunctor<IdentityOperator,
          FeH1,
          Double>(feFct, res2));
    }
    resultFunctors_[FLUIDMECH_PRESSURE] = vFunc;
    fieldFunctors_[FLUIDMECH_PRESSURE] = vFunc;
#endif
  }
  
  
  void FluidMechPerturbedPDE::DefinePostProcResults() {
#if 0
    shared_ptr<BaseFeFunction> feFct = feFunctions_[FLUIDMECH_VELOCITY];
    
    // Electric Field Intensity
    // create new resultDof object
    shared_ptr<ResultInfo> ef ( new ResultInfo );
    ef->resultType = ELEC_FIELD_INTENSITY;
    ef->SetVectorDOFs(dim_, isaxi_);
    ef->unit = "V/m";
    ef->definedOn = ResultInfo::ELEMENT;
    ef->entryType = ResultInfo::VECTOR;
    availResults_.insert( ef );
    postProcResults_[ELEC_FIELD_INTENSITY] = FLUIDMECH_VELOCITY;
    shared_ptr<BaseFieldFunctor> eFunc;
    if( isComplex_ ) {
      eFunc.reset(new DiffFieldFunctor<Complex>(feFct, ef));
    } else {
      eFunc.reset(new DiffFieldFunctor<Double>(feFct, ef));
    }
    resultFunctors_[ELEC_FIELD_INTENSITY] = eFunc;
    fieldFunctors_[ELEC_FIELD_INTENSITY] = eFunc;
    
    // Electric Flux Density
    shared_ptr<ResultInfo> flux ( new ResultInfo );
    flux->resultType = ELEC_FLUX_DENSITY;
    flux->SetVectorDOFs(dim_, isaxi_);
    flux->unit = "C/m^2";
    flux->definedOn = ResultInfo::ELEMENT;
    flux->entryType = ResultInfo::VECTOR;
    availResults_.insert( flux );
    postProcResults_[ELEC_FLUX_DENSITY] = FLUIDMECH_VELOCITY;
    shared_ptr<BaseFieldFunctor> fluxFunc;
    if( isComplex_ ) {
      fluxFunc.reset(new FluxFieldFunctor<Complex>(feFct, ef));
    } else {
      fluxFunc.reset(new FluxFieldFunctor<Double>(feFct, ef));
    }
    resultFunctors_[ELEC_FLUX_DENSITY] = fluxFunc;
    fieldFunctors_[ELEC_FLUX_DENSITY] = fluxFunc;

    // Electric charge
    shared_ptr<ResultInfo> charge( new ResultInfo );
    charge->resultType = ELEC_CHARGE;
    charge->definedOn = ResultInfo::SURF_ELEM;
    charge->entryType = ResultInfo::SCALAR;
    charge->dofNames = "";
    charge->unit = "C";
    availResults_.insert( charge );
    postProcResults_[ELEC_CHARGE] = FLUIDMECH_VELOCITY;

    // Electric Energy Density
    shared_ptr<ResultInfo> ed ( new ResultInfo );
    ed->resultType = ELEC_ENERGY_DENSITY;
    ed->dofNames = "";
    ed->unit = "Ws/m^3";
    ed->definedOn = ResultInfo::ELEMENT;
    ed->entryType = ResultInfo::SCALAR;
    availResults_.insert( ed );
    postProcResults_[ELEC_ENERGY_DENSITY] = FLUIDMECH_VELOCITY;
    shared_ptr<BaseFieldFunctor> edFunc;
    if( isComplex_ ) {
      edFunc.reset(new EnergyDensFieldFunctor<Complex>(feFct, ed));
    } else {
      edFunc.reset(new EnergyDensFieldFunctor<Double>(feFct, ed));
    }
    resultFunctors_[ELEC_ENERGY_DENSITY] = edFunc;
    fieldFunctors_[ELEC_ENERGY_DENSITY] = edFunc;
    
    // Electric energy
    shared_ptr<ResultInfo> energy( new ResultInfo );
    energy->resultType = ELEC_ENERGY;
    energy->definedOn = ResultInfo::REGION;
    energy->entryType = ResultInfo::SCALAR;
    energy->dofNames = "";
    energy->unit = "Ws";
    availResults_.insert ( energy );
    postProcResults_[ELEC_ENERGY] = FLUIDMECH_VELOCITY;
    shared_ptr<ResultFunctor> energyFunc;
    if( isComplex_ ) {
      energyFunc.reset(new EnergyResultFunctor<Complex>(feFct, energy));
    } else {
      energyFunc.reset(new EnergyResultFunctor<Double>(feFct, energy));
    }
    resultFunctors_[ELEC_ENERGY] = energyFunc;
    

    // Electric polarization
    shared_ptr<ResultInfo> pol( new ResultInfo );
    pol->resultType = ELEC_POLARIZATION;
    pol->definedOn = ResultInfo::ELEMENT;
    pol->entryType = ResultInfo::VECTOR;
    pol->SetVectorDOFs(dim_, isaxi_);
    pol->unit = "C/m^2";
    availResults_.insert( pol );
    postProcResults_[ELEC_POLARIZATION] = FLUIDMECH_VELOCITY;

    // pesudo electric polarization for piezo simp
    shared_ptr<ResultInfo> pseudoPol( new ResultInfo );
    pseudoPol->resultType = ELEC_PSEUDO_POLARIZATION;
    pseudoPol->definedOn = ResultInfo::ELEMENT;
    pseudoPol->entryType = ResultInfo::SCALAR;
    pseudoPol->dofNames = "";
    pseudoPol->unit = "";
    availResults_.insert( pseudoPol );
    postProcResults_[ELEC_PSEUDO_POLARIZATION] = FLUIDMECH_VELOCITY;

    // ============================
    // Initialize result functors:
    // ============================
    // 1) Loop over all BDB-integrators
    std::map<RegionIdType, BaseBDBInt*>::iterator it = bdbInts_.begin();
    for( ; it != bdbInts_.end(); ++it ) {
      RegionIdType region = it->first;
      BaseBDBInt* bdb = it->second;

      // 2) pass integrators to functors
      eFunc->AddIntegrator(bdb, region);
      fluxFunc->AddIntegrator(bdb, region);
      energyFunc->AddIntegrator(bdb, region);
      edFunc->AddIntegrator(bdb, region);
    }
#endif
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
}
