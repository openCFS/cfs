#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <string>
#include <set>

#include "ElecPDE.hh"

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

  DECLARE_LOG(elecpde)
  DEFINE_LOG(elecpde, "pde.electrostatic")

  // ***************
  //   Constructor
  // ***************
  ElecPDE::ElecPDE( Grid* aptgrid, PtrParamNode paramNode )
    :SinglePDE( aptgrid, paramNode ) {


    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "electrostatic";
    pdematerialclass_ = ELECTROSTATIC;
    maxTimeDerivOrder_ = 0;
 
    nonLin_    = false;
    nonLinMaterial_ = false;
    isAlwaysStatic_ = true;
    isPiezoCoupled_ = false;
    isThermoCoupled_= false;

    needSolPrev_ = true;
 
    // Check the subtype of the problem
    paramNode->GetValue("subType", subType_);
  }
  
  void ElecPDE::InitNonLin() {

    SinglePDE::InitNonLin();

    //now do PDE specifics
    std::map<std::string, NonLinType>::iterator it;
    for ( it=nonLinTypes_.begin() ; it != nonLinTypes_.end(); it++ ) {
      if ( (*it).second == HYSTERESIS ) {
        isHysteresis_ = true;
      }
    }

  }

  void ElecPDE::DefineIntegrators() {

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

    // if the pde is piezo-coupled, the electrostatic entries
    // have to be multiplied with -1
    std::string factor = "1.0";
    if ( isPiezoCoupled_ == true || isThermoCoupled_== true )
      factor = "-1.0";  

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    shared_ptr<FeSpace> mySpace = feFunctions_[ELEC_POTENTIAL]->GetFeSpace();
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      
      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;

      // Get current region name
      std::string regionName = ptgrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );
      
      // ==========================================================
      // New implementation
      // ==========================================================
      // --- Set the approximation for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId,integId);

      // --- standard real-valued stiffness integrator ---
      shared_ptr<CoefFunction > curCoef = actSDMat->GetCoefFunction(ELEC_PERMITTIVITY,tensorType,Global::REAL);


      BDBInt< GradientOperator ,FeH1,Double,Double >* stiffInt;
      stiffInt = new BDBInt<GradientOperator,FeH1,Double,Double >(curCoef,1.0 );
      //linElecInt *  linElecForm = new linElecInt( actSDMat, tensorType,
       //                                           upLagrangeForm );
      //linElecForm->SetFactor( factor );
      BiLinFormContext * stiffIntDescr =
          new BiLinFormContext(stiffInt, STIFFNESS );

      feFunctions_[ELEC_POTENTIAL]->AddEntityList( actSDList );

      //stiffIntDescr->SetPtPdes(this, this);
      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions( feFunctions_[ELEC_POTENTIAL],feFunctions_[ELEC_POTENTIAL]);
      stiffInt->SetFeSpace( feFunctions_[ELEC_POTENTIAL]->GetFeSpace());

      assemble_->AddBiLinearForm( stiffIntDescr );
      // Important: Add bdb-integrator to global list, as we need them later
      // for calculation of postprocessing results
      bdbInts_[actRegion] = stiffInt;
      
    }
    
    // Define integrators for composite materials
    // (only for subType "flatShell")
    std::map<RegionIdType, Composite>::iterator compIt;
    for( compIt=compositeMaterials_.begin(); compIt!=compositeMaterials_.end();
         compIt++ ) {
      
      REFACTOR;
    }

    // **********************************************************************
    //   inhom. Neumann boundary condition
    // **********************************************************************
    for( UInt iBc = 0; iBc < inBcs_.GetSize(); iBc++ ) 
    {
      REFACTOR;
    }
    
    // =======================================================================
    // Integrators for NonConforming Interfaces
    // =======================================================================
    ///OUT FOR REFACTOR
    // define integrators for electric impedances
    DefineImpedanceIntegrators();
  }
  
  void ElecPDE::DefineImpedanceIntegrators() {
    
    // The definition of the impedances is taken from:

    // Jian S. Wang and Dale F. Ostergaard:
    // A Finite Element-Electric Circuit Coupled Simulation Method for
    // Piezoelectric Transducer
    // IEEE Ultrasoncics Symposium, 1999

    // =======================================================================
    // Integrators for electric impedances
    // =======================================================================
    
    for( UInt i = 0; i < impedances_.GetSize(); i++ ) {
      
      REFACTOR;
    }
  }

  void ElecPDE::DefineSolveStep() {
    solveStep_ = new SolveStepElec(*this);
  }

  void ElecPDE::ReadSpecialBCs( ) 
  {
     //ReadImpedances();
  }
  
  void ElecPDE::ReadImpedances( ) 
  {
    REFACTOR;
  }
  
  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================

  // ***************
  //   CalcResults
  // ***************
  void ElecPDE::CalcResults( shared_ptr<BaseResult> res ) {

    switch (res->GetResultInfo()->resultType ) {
      case ELEC_POTENTIAL:
        feFunctions_[ELEC_POTENTIAL]->ExtractResult( res );
        break;

      case ELEC_RHS_LOAD:
        rhsFeFunctions_[ELEC_POTENTIAL]->ExtractResult( res );
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
    }

  }
   
  void ElecPDE::CalcPolarizationField( shared_ptr<BaseResult> res )
  {
    REFACTOR;
  }



  template <class TYPE>
  void ElecPDE::CalcCharges( shared_ptr<BaseResult> res ) {
    REFACTOR;
    // This will be also moved to the resultFunctor class,
    // as soon as the surface B-operators are defined
  }


  // ======================================================
  // COUPLING SECTION
  // ======================================================

  void ElecPDE::InitCoupling(PDECoupling * Coupling)
  {
  REFACTOR;
  }
  


  void ElecPDE::CalcOutputCoupling()
  {
    REFACTOR;
  }


  bool ElecPDE::HasOutput(SolutionType output)
  {
  
    switch (output)
      {
      case ELEC_FORCE_VWP:
        return true;
        break;
      case ELEC_POTENTIAL:
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


  void ElecPDE::SetPiezoCoupling()
  {
  
    isPiezoCoupled_ = true;

  }
  
  void ElecPDE::SetThermoCoupling()
  {
  
    isThermoCoupled_ = true;

  }
  
  void ElecPDE::DefinePolarizationMatrixIntegrators(const Vector<Double> &vals,
      StdVector<LinearFormContext*> *linForms, const int num)
  {
    REFACTOR;
  }

  void ElecPDE::DefinePrimaryResults() {
    
    shared_ptr<BaseFeFunction> feFct = feFunctions_[ELEC_POTENTIAL];
    
    // Electric Potential
    shared_ptr<ResultInfo> res1( new ResultInfo);
    res1->resultType = ELEC_POTENTIAL;
    
    // check for special subtype 
    if( subType_  != "flatShell" ) {
      res1->dofNames = "";
      res1->unit = "V";
      res1->definedOn = ResultInfo::NODE;
      res1->entryType = ResultInfo::SCALAR;
    } else {

      // check number of composite materials:
      // Currently, we can handle just one electrostatic composite
      // material, as we would have a different number of electric
      // unknowns for different regions
      if( compositeMaterials_.size() > 1 ) {
        WARN("Currently the electrostatic flatShell PDE can only "
            "handle ONE type of composite material!");
      }

      Composite & actComp = compositeMaterials_.begin()->second;
      UInt numLaminas = actComp.thickness.GetSize();
      for( UInt i=0; i<numLaminas; i++ ) {
        std::string dofName = "ep";
        dofName += lexical_cast<std::string>(i+1);
        res1->dofNames.Push_back( dofName );
      }
      res1->unit = "V";
      res1->definedOn = ResultInfo::ELEMENT;
      res1->entryType = ResultInfo::SCALAR;
    }
    feFunctions_[ELEC_POTENTIAL]->SetResultInfo(res1);

    res1->SetFeFunction(feFunctions_[ELEC_POTENTIAL]);
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
    resultFunctors_[ELEC_POTENTIAL] = vFunc;
    fieldFunctors_[ELEC_POTENTIAL] = vFunc;
    
    

    // Electric RHS Load
    // create new resultDof object
    shared_ptr<ResultInfo> rhs ( new ResultInfo );
    rhs->resultType = ELEC_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = "C";
    rhs->definedOn = results_[0]->definedOn;
    rhs->entryType = ResultInfo::SCALAR;
    availResults_.insert( rhs );
    postProcResults_[ELEC_RHS_LOAD] = ELEC_POTENTIAL;

    // ===================================
    // Check for non-conforming interfaces
    // ===================================
    StdVector<std::string> ncIfaceNames, ncIfaceNamesForPDE;
    StdVector<RegionIdType> ncIfaceIds;
    
    LOG_DBG2(elecpde) << "NonMatching: Checking if nonconforming "
                      << "interfaces of PDE exist in domain.";

    PtrParamNode elecPDENCIfaceListNode;
    elecPDENCIfaceListNode = param->GetByVal("sequenceStep", std::string("index"), sequenceStep_)
    ->Get("pdeList/electrostatic/ncInterfaceList", ParamNode::PASS);
    
    if(!elecPDENCIfaceListNode)
      return;

    PtrParamNode domainNCIfaceListNode;
    domainNCIfaceListNode = param->Get("domain")->Get("ncInterfaceList", ParamNode::PASS);

    if(!domainNCIfaceListNode)
    {
      EXCEPTION("No nonmatching interfaces have been specified in domain!");
    }

    ParamNodeList pdeNCIfaceNodes;
    pdeNCIfaceNodes = elecPDENCIfaceListNode->GetList("ncInterface");

    for (UInt i = 0; i < pdeNCIfaceNodes.GetSize(); i++) {
      std::string pdeIfaceName = pdeNCIfaceNodes[i]->Get("name")->As<std::string>();

      PtrParamNode domainIfaceNode = domainNCIfaceListNode
          ->GetByVal("ncInterface", "name", pdeIfaceName, ParamNode::PASS);
      if(!domainIfaceNode)
      {
        LOG_DBG2(elecpde) << "NonMatching: Nonconforming "
        << "interface '" << ncIfaceNames[i]
                                         << "' does not exist in domain.";

        EXCEPTION( "ncInterface referenced from PDE not defined in domain!");
      }

      ncIfaceNamesForPDE.Push_back(pdeIfaceName);
    }
    ptgrid_->GetRegion().Parse(ncIfaceNamesForPDE, ncIfaceIds);

    for (UInt i = 0; i < ncIfaceIds.GetSize(); i++) {
      ncIFaces_.Push_back(ncIfaceIds[i]);
    }

    // In the case of the presence of non-conforming interfaces,
    // a second resultdof object has to be created, which describes the 
    // Lagrange multiplier
    if( ncIFaces_.GetSize() > 0 ) {
      LOG_DBG2(elecpde) << "NonMatching: Defining new ResultDof Lagrange.";
      shared_ptr<ResultInfo> lagr ( new ResultInfo );
      lagr->resultType = LAGRANGE_MULT;
      lagr->dofNames = "l";
      lagr->definedOn = results_[0]->definedOn;
      results_.Push_back( lagr );
    } 
  }
  
  
  void ElecPDE::DefinePostProcResults() {

    shared_ptr<BaseFeFunction> feFct = feFunctions_[ELEC_POTENTIAL];
    
    // Electric Field Intensity
    // create new resultDof object
    shared_ptr<ResultInfo> ef ( new ResultInfo );
    ef->resultType = ELEC_FIELD_INTENSITY;
    ef->SetVectorDOFs(dim_, isaxi_);
    ef->unit = "V/m";
    ef->definedOn = ResultInfo::ELEMENT;
    ef->entryType = ResultInfo::VECTOR;
    availResults_.insert( ef );
    postProcResults_[ELEC_FIELD_INTENSITY] = ELEC_POTENTIAL;
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
    postProcResults_[ELEC_FLUX_DENSITY] = ELEC_POTENTIAL;
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
    postProcResults_[ELEC_CHARGE] = ELEC_POTENTIAL;

    // Electric Energy Density
    shared_ptr<ResultInfo> ed ( new ResultInfo );
    ed->resultType = ELEC_ENERGY_DENSITY;
    ed->dofNames = "";
    ed->unit = "Ws/m^3";
    ed->definedOn = ResultInfo::ELEMENT;
    ed->entryType = ResultInfo::SCALAR;
    availResults_.insert( ed );
    postProcResults_[ELEC_ENERGY_DENSITY] = ELEC_POTENTIAL;
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
    postProcResults_[ELEC_ENERGY] = ELEC_POTENTIAL;
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
    postProcResults_[ELEC_POLARIZATION] = ELEC_POTENTIAL;

    // pesudo electric polarization for piezo simp
    shared_ptr<ResultInfo> pseudoPol( new ResultInfo );
    pseudoPol->resultType = ELEC_PSEUDO_POLARIZATION;
    pseudoPol->definedOn = ResultInfo::ELEMENT;
    pseudoPol->entryType = ResultInfo::SCALAR;
    pseudoPol->dofNames = "";
    pseudoPol->unit = "";
    availResults_.insert( pseudoPol );
    postProcResults_[ELEC_PSEUDO_POLARIZATION] = ELEC_POTENTIAL;

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
  }

  
  std::map<SolutionType, shared_ptr<FeSpace> >
  ElecPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode) {
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H1"){
      PtrParamNode potSpaceNode = infoNode->Get("elecPotential");
      crSpaces[ELEC_POTENTIAL] =
        FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1);
      crSpaces[ELEC_POTENTIAL]->Init();
    }else{
      EXCEPTION("The formulation " << formulation << "of electric PDE is not known!");
    }
    return crSpaces;
  }

  LinearFormContext* ElecPDE::CreateRhsLinearForm(SolutionType rhsType,shared_ptr<CoefFunction > rhsCoef){
    LinearFormContext * mContext = NULL;
    switch(rhsType){
    case ELEC_CHARGE_DENSITY:
      BUIntegrator<IdentityOperator,FeH1,Double>* curInt;
      curInt = new BUIntegrator<IdentityOperator,FeH1,Double>(1.0,rhsCoef);

      mContext = new LinearFormContext(curInt);
      mContext->SetFeFunction( feFunctions_[ELEC_POTENTIAL]);
      curInt->SetFeSpace( feFunctions_[ELEC_POTENTIAL]->GetFeSpace());
      break;
    default:
      Exception("Right hand side quantity not known for elecPDE");
    }
    return mContext;
  }
}
