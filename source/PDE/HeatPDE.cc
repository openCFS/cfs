#include <boost/lexical_cast.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

#include <typeinfo>

#include "HeatPDE.hh"

#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionApprox.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Domain/CoefFunction/CoefFunctionMapping.hh"
#include "Domain/CoefFunction/CoefFunctionComplexToReal.hh"
#include "Domain/CoefFunction/CoefFunctionSUPG.hh"
#include "Utils/StdVector.hh"

#include "Driver/Assemble.hh"
#include "Utils/ApproxData.hh"
#include "Utils/SmoothSpline.hh"
#include "Materials/Models/Hysteresis.hh"
#include "Materials/Models/Preisach.hh"
//#include "Materials/Models/VectorPreisach.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/FeFunctions.hh"

#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

#include "Driver/SolveSteps/StdSolveStep.hh"


#include "Driver/Assemble.hh"

//new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/BiLinForms/BiLinWrappedLinForm.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/ConvectiveOperator.hh"
#include "Forms/LinForms/SingleEntryInt.hh"

//new postprocessing concept
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/Results/ResultFunctor.hh"

namespace CoupledField {

  DEFINE_LOG(heatcondpde, "pde.heatcond")

// ======================================================
// SET SOLUTION INFORMATION
// ======================================================
HeatPDE::HeatPDE(Grid * aptgrid, PtrParamNode paramNode,
                 PtrParamNode infoNode,
                 shared_ptr<SimState> simState, Domain* domain)
:SinglePDE( aptgrid, paramNode, infoNode, simState, domain ) {

  pdename_           = "heatConduction";
  pdematerialclass_  = THERMIC;
  nonLin_            = false;
  
  //! Always use updated Lagrangian formulation 
  updatedGeo_        = true;
  isLinFlowPDECoupled_ = false;
  isCouplingFormulationSymmetric_ = false;

  interfaceDrivenHeatSource_ = false;

  // convert to tensor type
  tensorType_ = FULL;
  subType_ = "3d";
  if ( ptGrid_->GetDim() == 2 ) {
    if ( ptGrid_->IsAxi() ) {
      tensorType_ = AXI;
      subType_ = "axi";
      isaxi_ = true;
    } else {
      tensorType_ = PLANE;
      subType_ = "plane";
    }
  }
}




  std::map<SolutionType, shared_ptr<FeSpace> >
  HeatPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode) {
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H1"){
      PtrParamNode potSpaceNode = infoNode->Get("heatTemperature");
      crSpaces[HEAT_TEMPERATURE] =
        FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[HEAT_TEMPERATURE]->Init(solStrat_);
    }else{
      EXCEPTION("The formulation " << formulation << "of heat PDE is not known!");
    }
    return crSpaces;
  }


  void HeatPDE::SetLinFlowPDECouplingFlags(bool useSymmtericForm) {
  	// Set flag for coupling
    isLinFlowPDECoupled_ = true;
    // Set flag whether to use symmetric form or not
    isCouplingFormulationSymmetric_= useSymmtericForm;
  }


  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void HeatPDE::InitNonLin() {

    SinglePDE::InitNonLin();
 //  nonLinTotalFormulation_ = true;

  }


  void HeatPDE::ReadDampingInformation() {
    std::map<std::string, DampingType> idDampType;

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // In this case the "damping" layer is an infinite mapping of the solution,
    // since it is a coordinate transformation it can be set up using the PML-like
    // framework.
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // try to get dampingList
    PtrParamNode dampListNode = myParam_->Get( "dampingList", ParamNode::PASS );
    if( dampListNode ) {

      // get specific damping nodes
      ParamNodeList dampNodes = dampListNode->GetChildren();

      for( UInt i = 0; i < dampNodes.GetSize(); i++ ) {

        std::string dampString = dampNodes[i]->GetName();
        std::string actId = dampNodes[i]->Get("id")->As<std::string>();

        // determine type of damping
        DampingType actType;
        String2Enum( dampString, actType );

        // store damping type string
        idDampType[actId] = actType;

      }
    }

    // Run over all region
    ParamNodeList regionNodes =
        myParam_->Get("regionList")->GetChildren();

    RegionIdType actRegionId;
    std::string actRegionName, actDampingId;

    for (UInt k = 0; k < regionNodes.GetSize(); k++) {
      regionNodes[k]->GetValue( "name", actRegionName );
      regionNodes[k]->GetValue( "dampingId", actDampingId );
      if( actDampingId == "" )
        continue;

      actRegionId = ptGrid_->GetRegion().Parse( actRegionName );

      // Check actDampingId was already registerd
      if( idDampType.count( actDampingId ) == 0 ) {
        EXCEPTION( "Damping with id '" << actDampingId
                   << "' was not defined in 'dampingList'" );
      }
      dampingList_[actRegionId] = idDampType[actDampingId];
    }
  }


void HeatPDE::DefineIntegrators() {
  RegionIdType actRegion;
  BaseMaterial * actSDMat = NULL;  

  // Define integrators for "standard" materials
  std::map<RegionIdType, BaseMaterial*>::iterator it;
  shared_ptr<FeSpace> mySpace = feFunctions_[HEAT_TEMPERATURE]->GetFeSpace();
  for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion++){
      actRegion = regions_[iRegion];
      actSDMat    = materials_[actRegion];

    // Get current region name
    std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

    // create new entity list
    shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
    actSDList->SetRegion( actRegion );

    // --- Set the FE ansatz for the current region ---
    PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
    std::string polyId = curRegNode->Get("polyId")->As<std::string>();
    std::string integId = curRegNode->Get("integId")->As<std::string>();
    mySpace->SetRegionApproximation(actRegion, polyId,integId);

    // pass entitylist of fespace / fefunction
    shared_ptr<BaseFeFunction> feFunc = feFunctions_[HEAT_TEMPERATURE];
    feFunc->AddEntityList( actSDList );

    // Create coefficient function T_ref to symmetries the system, by dividing the all  biliniear and linear from integrators
    PtrCoefFct refTemp = nullptr;
    if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
      refTemp  = actSDMat->GetScalCoefFnc(HEAT_REF_TEMPERATURE, Global::REAL);
    }

    // ====================================================================
    // Take account for mapping of an infinite domain
    // ====================================================================
    // Generate scalar valued coefficient function
    PtrCoefFct val1 = CoefFunction::Generate( mp_, Global::REAL, "1.0");
    shared_ptr<CoefFunction> coeffMAPScal, coeffMAPVec;
    bool isMapping = false;
    if( dampingList_[actRegion] == MAPPING ) {
      std::string dampId;
      curRegNode->GetValue("dampingId",dampId);
      if(analysistype_ == HARMONIC){
        EXCEPTION("Harmonic analysis not allowed!");
      }else{
        PtrParamNode mapNode = myParam_->Get("dampingList")->GetByVal("mapping","id",dampId.c_str());
        coeffMAPVec.reset(new CoefFunctionMapping<Double>(mapNode,val1,actSDList,regions_,true));
        coeffMAPScal.reset(new CoefFunctionMapping<Double>(mapNode,val1,actSDList,regions_,false));
        isMapping = true;
      }
    }


    // ====================================================================
    // stiffness integrator
    // ====================================================================

    //get possible nonlinearities defined in this region
    StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];
    if ( nonLinTypes.Find(NLHEAT_CONDUCTIVITY) != -1 ) {
      // ====================================================================
      // nonlinear stiffness integrator
      // ====================================================================
      // no infinite mapping is implemented for the nonlinear case but it can be
      // analogue to the linear case
      if(isMapping){
        EXCEPTION("Infinite mapping is not implemented for the nonlinear HeatPDE!!")
      }

      PtrCoefFct heatCoef = this->GetCoefFct(HEAT_TEMPERATURE);
      PtrCoefFct condNL = actSDMat->GetScalCoefFncNonLin( HEAT_CONDUCTIVITY_SCALAR, Global::REAL, heatCoef);
      if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
        condNL = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, condNL, refTemp, CoefXpr::OP_DIV));
      }

      // create stiffness integrator
      BaseBDBInt* stiffInt = NULL;
      if( dim_ == 2 ) {
        stiffInt = new BBInt<>(new GradientOperator<FeH1,2>(), condNL, 1.0, updatedGeo_);
      } else {
        stiffInt = new BBInt<>(new GradientOperator<FeH1,3>(), condNL, 1.0, updatedGeo_);
      }
      stiffInt->SetName("StiffnessIntegrator-NL");

      BiLinFormContext* stiffContext = new BiLinFormContext(stiffInt, STIFFNESS);
      stiffContext->SetEntities( actSDList, actSDList );
      stiffContext->SetFeFunctions( feFunc, feFunc );

      assemble_->AddBiLinearForm( stiffContext );
      bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffInt) );

      // =================================
      //  Nonlinear RHS-integrator
      // =================================
      //      LinearForm * rhsNlinForm = new KXIntegrator<Double>(stiffInt, -1.0,
      //                                                          feFunc );
      //      rhsNlinForm->SetName("RHSNonLinFormHeatStiff");
      //      LinearFormContext * rhsNlinContext =
      //          new LinearFormContext( rhsNlinForm );
      //      rhsNlinContext->SetEntities( actSDList );
      //      rhsNlinContext->SetFeFunction( feFunc );
      //      assemble_->AddLinearForm( rhsNlinContext );
    }
    else {
      // ====================================================================
      // linear stiffness integrator
      // ====================================================================
      shared_ptr<CoefFunction > curCoef = actSDMat->GetTensorCoefFnc( HEAT_CONDUCTIVITY_TENSOR, tensorType_, Global::REAL );
      if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
        curCoef = CoefFunction::Generate(mp_, Global::REAL, CoefXprTensScalOp(mp_, curCoef, refTemp, CoefXpr::OP_DIV));
      }

      // when we do optimization we wrap the original CoefFunction. Don't check for region to handle dim-1 pressure on dim elements
      if(domain->HasDesign()) {
        CoefFunctionOpt* tmpFnc = new CoefFunctionOpt(domain->GetDesign(), curCoef, HEAT_CONDUCTIVITY_TENSOR, this); // takes double and complex
        curCoef.reset(tmpFnc);
      }

      BaseBDBInt* stiffInt = NULL;
      if(isMapping){
        // define an infinite mapping layer
        PtrCoefFct curCoefScl = CoefFunction::Generate(mp_, Global::REAL, CoefXprTensScalOp(mp_, curCoef, coeffMAPScal, CoefXpr::OP_MULT_TENSOR));
        if(dim_ == 2)
          stiffInt = new BDBInt<Double>(new ScaledGradientOperator<FeH1,2,Double>(), curCoefScl, 1.0, updatedGeo_ );
        else
          stiffInt = new BDBInt<Double>(new ScaledGradientOperator<FeH1,3,Double>(), curCoefScl, 1.0, updatedGeo_ );

        stiffInt->SetBCoefFunctionOpA(coeffMAPVec);
      }else{
        // no mapping
        if(dim_ == 2)
          stiffInt = new BDBInt<>(new GradientOperator<FeH1,2>(), curCoef,1.0, updatedGeo_);
        else
          stiffInt = new BDBInt<>(new GradientOperator<FeH1,3>(), curCoef,1.0, updatedGeo_);
      }

      stiffInt->SetName("HeatConductivity");
      // the integrator has a coef function but for the optimization case the opt coef needs to know also the integrator
      if(domain->HasDesign())
        dynamic_pointer_cast<CoefFunctionOpt>(curCoef)->SetForm(stiffInt);

      BiLinFormContext* stiffIntDescr = new BiLinFormContext(stiffInt, STIFFNESS );

      //stiffIntDescr->SetPtPdes(this, this);
      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions(feFunc,feFunc);
      stiffInt->SetFeSpace( feFunctions_[HEAT_TEMPERATURE]->GetFeSpace());

      assemble_->AddBiLinearForm( stiffIntDescr );
      bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,stiffInt) );

      //      if ( nonLinTypes.Find(NLHEAT_CONDUCTIVITY) ||
      //                nonLinTypes.Find(NLHEAT_CAPACITY) != -1 ) {
      //        // === Additional RHS integrator in case of Non-linearity ===
      //        LinearForm * rhsNlinForm = new KXIntegrator<Double>(stiffInt, -1.0,
      //                                                            feFunc );
      //        rhsNlinForm->SetName("RHSNonLinFormHeatStiff-Lin");
      //        LinearFormContext * rhsNlinContext =
      //            new LinearFormContext( rhsNlinForm );
      //        rhsNlinContext->SetEntities( actSDList );
      //        rhsNlinContext->SetFeFunction( feFunc );
      //        assemble_->AddLinearForm( rhsNlinContext );
      //      }
    }

    StabilisationType stabilisation = NO_STABILISATION;
    std::string velocityId = curRegNode->Get("velocityId")->As<std::string>();
    if(velocityId != "") {
      // Add the region information
      PtrParamNode velNode = myParam_->Get("velocityList")->GetByVal("velocity","name",velocityId.c_str());
      stabilisation = BasePDE::stabilisationType.Parse(velNode->Get("stabilisation")->As<std::string>());
    }

    // ====================================================================
    // mass integrator
    // ====================================================================
    if ( nonLinTypes.Find(NLHEAT_CAPACITY) != -1 ) {
      // ====================================================================
      // nonlinear mass integrator
      // ====================================================================
      // Factor for mass matrix: density * heatCapacity
      PtrCoefFct heatCoef = this->GetCoefFct(HEAT_TEMPERATURE);
      PtrCoefFct density = NULL;
      if(nonLinTypes.Find(NLHEAT_DENSITY) != -1){
        // nonlinear (temperature-dependent) mass density
        density = actSDMat->GetScalCoefFncNonLin( DENSITY, Global::REAL, heatCoef );
      }else{
        // linear mass density
        density = actSDMat->GetScalCoefFnc( DENSITY, Global::REAL );
      }

      //BaseBOperator * bOp = new IdentityOperator<FeH1>();
      PtrCoefFct capNL = actSDMat->GetScalCoefFncNonLin( HEAT_CAPACITY, Global::REAL, heatCoef );
      PtrCoefFct nlMassCoeff = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, capNL, density, CoefXpr::OP_MULT));
      //Symmetries the forumlation by divinding with Temp_ref
      if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
        nlMassCoeff = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp( mp_, nlMassCoeff, refTemp, CoefXpr::OP_DIV));
      }

      // create stiffness integrator
      BaseBDBInt* massIntNL = NULL;
      if(dim_ == 2)
        massIntNL = new BBInt<>(new IdentityOperator<FeH1,2>(), nlMassCoeff, 1, updatedGeo_ );
      else
        massIntNL = new BBInt<>(new IdentityOperator<FeH1,3>(), nlMassCoeff, 1, updatedGeo_ );

      massIntNL->SetName("MassIntegrator-NL");

      BiLinFormContext * massNLContext = new BiLinFormContext(massIntNL, DAMPING );
      massNLContext->SetEntities( actSDList, actSDList );
      massNLContext->SetFeFunctions( feFunc, feFunc );

      assemble_->AddBiLinearForm( massNLContext );
      //bdbInts_[actRegion] = massIntNL;

    }
    else {
      // ====================================================================
      // linear mass integrator
      // ====================================================================
      // Factor for mass matrix: density * heatCapacity
      PtrCoefFct density = actSDMat->GetScalCoefFnc( DENSITY, Global::REAL );
      PtrCoefFct heatCapacity = actSDMat->GetScalCoefFnc( HEAT_CAPACITY, Global::REAL );
      PtrCoefFct massFactor = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp( mp_, density, heatCapacity, CoefXpr::OP_MULT));
      if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
        massFactor = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp( mp_, massFactor, refTemp, CoefXpr::OP_DIV));
      }

      if(domain->HasDesign()) {
        CoefFunctionOpt* tmpFnc = new CoefFunctionOpt(domain->GetDesign(), massFactor, NO_MATERIAL, this); // No single material property
        massFactor.reset(tmpFnc);
      }

      BiLinearForm *massInt = NULL;
      if(dim_ == 2)
        massInt = new BBInt<>(new IdentityOperator<FeH1,2,1,Double>(), massFactor,1.0, updatedGeo_ );
      else
        massInt = new BBInt<>(new IdentityOperator<FeH1,3,1,Double>(), massFactor,1.0, updatedGeo_ );

      massInt->SetName("MassIntegrator");
      // the integrator has a coef function but for the optimization case the opt coef needs to know also the integrator
      if(domain->HasDesign())
        dynamic_pointer_cast<CoefFunctionOpt>(massFactor)->SetForm(massInt);
      massInt->SetFeSpace( feFunc->GetFeSpace() );

      BiLinFormContext *massContext =  new BiLinFormContext(massInt, DAMPING );
      massContext->SetEntities( actSDList, actSDList );
      massContext->SetFeFunctions( feFunc,feFunc);
      assemble_->AddBiLinearForm( massContext );

      switch (stabilisation) 
      {
        case StabilisationType::SUPG:
        {
          EXCEPTION("SUPG is not validated and tested for Mass matrix! Use ArtificialDiffusion or validate the implementation below ...");
          // get material property for tau calculation (we need to get ) 
          shared_ptr<CoefFunction > curCoef = actSDMat->GetTensorCoefFnc( HEAT_CONDUCTIVITY_TENSOR, tensorType_, Global::REAL );
          if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
            curCoef = CoefFunction::Generate(mp_, Global::REAL, CoefXprTensScalOp(mp_, curCoef, refTemp, CoefXpr::OP_DIV));
          }
          // when we do optimization we wrap the original CoefFunction. Don't check for region to handle dim-1 pressure on dim elements
          if(domain->HasDesign()) {
            CoefFunctionOpt* tmpFnc = new CoefFunctionOpt(domain->GetDesign(), curCoef, HEAT_CONDUCTIVITY_TENSOR, this); // takes double and complex
            curCoef.reset(tmpFnc);
          }
          curCoef = CoefFunction::Generate(mp_, Global::REAL, CoefXprTensScalOp(mp_, curCoef, massFactor, CoefXpr::OP_DIV));
          // Calculate the stabilisation coeffitient here
          PtrCoefFct coeffUpwindingFactor;
          coeffUpwindingFactor.reset(new CoefFunctionSUPG(convecVelCoef_, curCoef, feFunc));
          // combination for material properites and tau coeffitient 
          PtrCoefFct factor = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp( mp_, coeffUpwindingFactor, massFactor, CoefXpr::OP_MULT )) ;
          // Create the integrators
          BiLinearForm *massIntSUPG = NULL;
          BaseBOperator *bOperator = NULL;
          if(dim_ == 2)
          {
            bOperator = new ConvectiveOperator<FeH1,2,1>();
            bOperator->SetCoefFunction(convecVelCoef_);
            massIntSUPG = new ABInt<>(new IdentityOperator<FeH1,2,1,Double>(), bOperator, factor,1.0, updatedGeo_ );
          }
          else
          {
            bOperator = new ConvectiveOperator<FeH1,3,1>();
            bOperator->SetCoefFunction(convecVelCoef_);
            massIntSUPG = new ABInt<>(new IdentityOperator<FeH1,3,1,Double>(), bOperator, factor,1.0, updatedGeo_ );
          }

          massIntSUPG->SetName("MassIntegratorSUPG");
          // the integrator has a coef function but for the optimization case the opt coef needs to know also the integrator
          if(domain->HasDesign())
            dynamic_pointer_cast<CoefFunctionOpt>(factor)->SetForm(massIntSUPG);
          massIntSUPG->SetFeSpace( feFunc->GetFeSpace() );

          BiLinFormContext *massContext =  new BiLinFormContext(massIntSUPG, DAMPING );
          massContext->SetEntities( actSDList, actSDList );
          massContext->SetFeFunctions( feFunc,feFunc);
          assemble_->AddBiLinearForm( massContext );
          break;
        }
        default:
        {break;}
      }
    }


    // ====================================================================
    // check for convective velocity (no infinite mapping allowed)
    // ====================================================================
    if(velocityId != "") {
      if(isMapping)
        EXCEPTION("Infinite mapping, applied to a region with defined velocity is permitted!!");

      // Get result info object for flow
      shared_ptr<ResultInfo> velInfo = GetResultInfo(MEAN_FLUIDMECH_VELOCITY);

      // Add the region information
      PtrParamNode velNode = myParam_->Get("velocityList")->GetByVal("velocity","name",velocityId.c_str());
     
      // Read velocity coefficient function for this region and add it to velocity functor
      PtrCoefFct regionMoving;
      std::set<UInt> definedDofs;
      bool coefUpdateGeo;
      //we assume that velocity is real
      ReadUserFieldValues( actSDList, velNode, velInfo->dofNames, velInfo->entryType, isComplex_, regionMoving, definedDofs, coefUpdateGeo );
      convecVelCoef_->AddRegion( actRegion, regionMoving );

      // Factor for convective matrix: density * heatCapacity
      PtrCoefFct density = actSDMat->GetScalCoefFnc( DENSITY, Global::REAL );
      PtrCoefFct heatCapacity = NULL;
      if ( nonLinTypes.Find(NLHEAT_CAPACITY) != -1 ) {
        PtrCoefFct heatCoef = this->GetCoefFct(HEAT_TEMPERATURE);
        heatCapacity = actSDMat->GetScalCoefFncNonLin( HEAT_CAPACITY, Global::REAL, heatCoef );
      }else{
        heatCapacity = actSDMat->GetScalCoefFnc( HEAT_CAPACITY, Global::REAL );
      }
      PtrCoefFct velFactor = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp( mp_, density, heatCapacity, CoefXpr::OP_MULT ) );
      if(domain->HasDesign()) {
        CoefFunctionOpt* tmpFnc = new CoefFunctionOpt(domain->GetDesign(), velFactor, NO_MATERIAL, this); // takes double and complex
        velFactor.reset(tmpFnc);
      }
      if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
        velFactor = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp( mp_, velFactor, refTemp, CoefXpr::OP_DIV));
      }

      // Create the integrators
      BaseBDBInt   *convectiveStiff = NULL;
      if( isComplex_ ) {
        if(dim_ == 2)
          convectiveStiff  = new ABInt<>(new IdentityOperator<FeH1,2,1>(),new ConvectiveOperator<FeH1,2,1,Complex>(), velFactor, 1.0, coefUpdateGeo);
        else
          convectiveStiff  = new ABInt<>(new IdentityOperator<FeH1,3,1>(),new ConvectiveOperator<FeH1,3,1,Complex>(), velFactor, 1.0, coefUpdateGeo);
      } else{
        if(dim_ == 2)
          convectiveStiff  = new ABInt<>(new IdentityOperator<FeH1,2,1>(),new ConvectiveOperator<FeH1,2,1>(),velFactor, 1.0, coefUpdateGeo);
        else
          convectiveStiff  = new ABInt<>(new IdentityOperator<FeH1,3,1>(),new ConvectiveOperator<FeH1,3,1>(),velFactor, 1.0, coefUpdateGeo);
      }

      convectiveStiff->SetBCoefFunctionOpB(convecVelCoef_);
      if ( nonLinTypes.Find(NLHEAT_CAPACITY) != -1 ) {
        convectiveStiff->SetName("ConvectiveStiffInt-NL");
      }else{
        convectiveStiff->SetName("ConvectiveStiffInt");
      }
      // the integrator has a coef function but for the optimization case the opt coef needs to know also the integrator
      if(domain->HasDesign())
        dynamic_pointer_cast<CoefFunctionOpt>(velFactor)->SetForm(convectiveStiff);

      convectiveInts_[actRegion] = convectiveStiff;

      BiLinFormContext *convectiveContextStiff =  new BiLinFormContext(convectiveStiff, STIFFNESS );
      convectiveContextStiff->SetEntities( actSDList, actSDList );
      convectiveContextStiff->SetFeFunctions( feFunctions_[HEAT_TEMPERATURE],feFunc);
      assemble_->AddBiLinearForm( convectiveContextStiff );

      PtrParamNode in = infoNode_->Get("velocity");
      in->Get("id")->SetValue(velocityId);

      switch (stabilisation)
      {
        case StabilisationType::SUPG:
        {
          EXCEPTION("SUPG is not implemented for all the terms of the algebraic system of equations! Use ArtificialDiffusion instead!");
          // IMPLEMENT check if an element is linear in a region for SUPG 
          // because if it's not linear, we would have one more term in stiffness matrix
        }
        case StabilisationType::ARTIFICIAL_DIFFUSION:
        {
          PtrCoefFct curCoef;
          if ( nonLinTypes.Find(NLHEAT_CONDUCTIVITY) != -1 ) {
            // NON-LINEAR MATERIAL
            // non-linear heat conductivity
            PtrCoefFct heatCoef = this->GetCoefFct(HEAT_TEMPERATURE);
            curCoef = actSDMat->GetScalCoefFncNonLin( HEAT_CONDUCTIVITY_SCALAR, Global::REAL, heatCoef);
            if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
              curCoef = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, curCoef, refTemp, CoefXpr::OP_DIV));
            }
            // material parameter - heat conductivity devided by capacity and density 
            curCoef = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, curCoef, velFactor, CoefXpr::OP_DIV));
          }
          else
          {
            // LINEAR MATERIAL
            // get material property for stabilisation factor calculation
            // heat conductivity
            curCoef = actSDMat->GetTensorCoefFnc( HEAT_CONDUCTIVITY_TENSOR, tensorType_, Global::REAL );
            if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
              curCoef = CoefFunction::Generate(mp_, Global::REAL, CoefXprTensScalOp(mp_, curCoef, refTemp, CoefXpr::OP_DIV));
            }
            // when we do optimization we wrap the original CoefFunction. Don't check for region to handle dim-1 pressure on dim elements
            if(domain->HasDesign()) {
              CoefFunctionOpt* tmpFnc = new CoefFunctionOpt(domain->GetDesign(), curCoef, HEAT_CONDUCTIVITY_TENSOR, this); // takes double and complex
              curCoef.reset(tmpFnc);
            }
            // material parameter - heat conductivity devided by capacity and density 
            curCoef = CoefFunction::Generate(mp_, Global::REAL, CoefXprTensScalOp(mp_, curCoef, velFactor, CoefXpr::OP_DIV));
          }
          
          // calculate stabilisation parameter
          PtrCoefFct coeffUpwindingFactor;
          coeffUpwindingFactor.reset(new CoefFunctionSUPG(convecVelCoef_, curCoef, feFunc));
          // combination for material properites and stabilisation coeffitient for bilinear terms
          PtrCoefFct factor = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp( mp_, coeffUpwindingFactor, velFactor, CoefXpr::OP_MULT )) ;
          // Create the integrators
          BaseBDBInt   *convectiveSUPG = NULL;
          BaseBOperator *bOperator = NULL;
          if( isComplex_ ) {
            if(dim_ == 2)
              bOperator = new ConvectiveOperator<FeH1,2,1,Complex>();
            else
              bOperator = new ConvectiveOperator<FeH1,3,1,Complex>();
          } else{
            if(dim_ == 2)
              bOperator = new ConvectiveOperator<FeH1,2,1>();
            else
              bOperator = new ConvectiveOperator<FeH1,3,1>();
          }
          bOperator->SetCoefFunction(convecVelCoef_);
          convectiveSUPG  = new BBInt<>(bOperator, factor, 1.0, coefUpdateGeo);
          // store the stabilisation term in the STIFFNESS matrix
          if ( nonLinTypes.Find(NLHEAT_CAPACITY) != -1 || nonLinTypes.Find(NLHEAT_CONDUCTIVITY) != -1) {
            convectiveSUPG->SetName("ConvectiveStiffSUPG-NL");
          }else{
            convectiveSUPG->SetName("ConvectiveStiffSUPG");
          }
               // the integrator has a coef function but for the optimization case the opt coef needs to know also the integrator
          if(domain->HasDesign())
            dynamic_pointer_cast<CoefFunctionOpt>(velFactor)->SetForm(convectiveSUPG);
          BiLinFormContext *convectiveContextStiffSUPG =  new BiLinFormContext(convectiveSUPG, STIFFNESS );
          convectiveContextStiffSUPG->SetEntities( actSDList, actSDList );
          convectiveContextStiffSUPG->SetFeFunctions( feFunctions_[HEAT_TEMPERATURE],feFunc);
          assemble_->AddBiLinearForm( convectiveContextStiffSUPG );
          break;
        }
        default:
        {
          break;
        }
      }
    } //end convective term
  } //end loop over materials_


  // ===============
  //  electric power density
  // ===============
  LOG_DBG(heatcondpde) << "Reading electric power densities";

  shared_ptr<BaseFeFunction> myFct = feFunctions_[HEAT_TEMPERATURE];
  StdVector<std::string> dispDofNames = myFct->GetResultInfo()->dofNames;
  StdVector<shared_ptr<EntityList> > ent;
  StdVector<PtrCoefFct > coef;
  LinearForm * lin = NULL;
  //LinearForm * linSUPG = NULL;
  StdVector<std::string> volumeRegions;

  ReadRhsExcitation( "elecPowerDensity", dispDofNames, ResultInfo::SCALAR, isComplex_, ent, coef, updatedGeo_ );
  for( UInt i = 0; i < ent.GetSize(); ++i ) {
    // check type of entitylist
    if (ent[i]->GetType() == EntityList::NODE_LIST) {
      EXCEPTION("Electric power density must be defined on elements");
    }

    if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
      RegionIdType entRegion = ent[i]->GetRegion();
      PtrCoefFct refTemp = materials_[entRegion]->GetScalCoefFnc(HEAT_REF_TEMPERATURE, Global::REAL);
      coef[i] = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp( mp_, coef[i], refTemp, CoefXpr::OP_DIV));
    }

    if( dim_ == 2) {
      if(isComplex_) {
        lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,2>(),Complex(1.0), coef[i], updatedGeo_ );
      } else {
        lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,2>(),1.0, coef[i], updatedGeo_ );
      }
    } else  {
      if(isComplex_) {
        lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,3>(), Complex(1.0), coef[i], updatedGeo_ );
      } else {
        lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,3>(), 1.0, coef[i], updatedGeo_ );
      }
    }
    lin->SetName("ElectricPowerDensityInt");
    LinearFormContext *ctx = new LinearFormContext( lin );
    ctx->SetEntities( ent[i] );
    ctx->SetFeFunction(myFct);
    assemble_->AddLinearForm(ctx);
    myFct->AddEntityList(ent[i]);

    // implementation of SUPG stabilisation for a convective term
    // it does not work because of the region implementation for curRegNodeID
    // Probably, the reason is that ent[i] is not a node, but an element
/*  
    BaseBOperator *bOperator = NULL;
    RegionIdType curRegNodeID = ent[i]->GetRegion();
    PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","id",curRegNodeID);

    StabilisationType stabilisation = NO_STABILISATION;
    std::string velocityId = curRegNode->Get("velocityId")->As<std::string>();
    if(velocityId != "") {
      // Add the region information
      PtrParamNode velNode = myParam_->Get("velocityList")->GetByVal("velocity","name",velocityId.c_str());
      stabilisation = BasePDE::stabilisationType.Parse(velNode->Get("stabilisation")->As<std::string>());
    }

    switch (stabilisation)
    {
      case StabilisationType::SUPG:
      {
        EXCEPTION("SUPG is not validated and tested for electric power input! Use ArtificialDiffusion");
        // Factor = density * heatCapacity
        PtrCoefFct velFactor = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp( mp_, actSDMat->GetScalCoefFnc( HEAT_CAPACITY, Global::REAL ), actSDMat->GetScalCoefFnc( DENSITY, Global::REAL ), CoefXpr::OP_MULT ) );
        // get tau
        shared_ptr<CoefFunction > curCoef = actSDMat->GetTensorCoefFnc( HEAT_CONDUCTIVITY_TENSOR, tensorType_, Global::REAL );
        curCoef = CoefFunction::Generate(mp_, Global::REAL, CoefXprTensScalOp(mp_, curCoef, velFactor, CoefXpr::OP_DIV));
        PtrCoefFct coeffUpwindingFactor;
        coeffUpwindingFactor.reset(new CoefFunctionSUPG(convecVelCoef_, curCoef, myFct));
        PtrCoefFct newCoef;
        newCoef = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp( mp_, coeffUpwindingFactor, coef[i], CoefXpr::OP_MULT )) ;

        if( isComplex_ )
        {
          if(dim_ == 2)
            bOperator = new ConvectiveOperator<FeH1,2,1,Complex>();
          else
            bOperator = new ConvectiveOperator<FeH1,3,1,Complex>();

          bOperator->SetCoefFunction(convecVelCoef_);
          linSUPG = new BUIntegrator<Complex> (bOperator, Complex(1.0), newCoef, updatedGeo_ );
        } else{
          if(dim_ == 2)
            bOperator = new ConvectiveOperator<FeH1,2,1>();
          else
            bOperator = new ConvectiveOperator<FeH1,3,1>();

          bOperator->SetCoefFunction(convecVelCoef_);
          linSUPG = new BUIntegrator<Double> (bOperator, 1.0, newCoef, updatedGeo_ );
        }
        linSUPG->SetName("ElectricPowerDensityIntSUPG");
        LinearFormContext *ctxSUPG = new LinearFormContext( linSUPG );
        ctxSUPG->SetEntities( ent[i] );
        ctxSUPG->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctxSUPG);
        myFct->AddEntityList(ent[i]);
        break;
      }
      default:
      {
        break; 
      }

    }
    */
  } // end loop over entities


  // ======================================================================
  // Heat flux boundary condition
  // ======================================================================
  ReadRhsExcitation( "heatFlux", dispDofNames, ResultInfo::SCALAR, isComplex_, ent, coef, updatedGeo_ );
  if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
    // Get the volumeRegions for the heatFlux, because the reference temperature of the volumeRegion is needed to make the system symmetric
    ReadVolumeRegions("heatFlux", volumeRegions);
    if( ent.GetSize() != volumeRegions.GetSize()){
      EXCEPTION("Define a volumeRegion for the heatFlux or set the attribute symmetric=\"false\" in <couplingList> -> <direct> -> <linFlowHeatDirect>");
    }
  }
  for( UInt i = 0; i < ent.GetSize(); ++i ) {
    // check type of entitylist
    if (ent[i]->GetType() == EntityList::NODE_LIST) {
      EXCEPTION("Heat flux must be defined on elements")
    }

    if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
      if(volumeRegions[i] == ""){
        EXCEPTION("Define a volumeRegion for the heatFlux or set the attribute symmetric=\"false\" in <couplingList> -> <direct> -> <linFlowHeatDirect>");
      }
      RegionIdType volRegion = ptGrid_->GetRegion().Parse(volumeRegions[i]);
      PtrCoefFct refTemp = materials_[volRegion]->GetScalCoefFnc(HEAT_REF_TEMPERATURE, Global::REAL);
      coef[i] = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp( mp_, coef[i], refTemp, CoefXpr::OP_DIV));
    }

    if( dim_ == 2) {
      if(isComplex_) {
        lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,2>(), Complex(1.0), coef[i], updatedGeo_ );
      } else {
        lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,2>(), 1.0, coef[i], updatedGeo_ );
      }
    } else  {
      if(isComplex_) {
        lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,3>(), Complex(1.0), coef[i], updatedGeo_ );
      } else {
        lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,3>(), 1.0, coef[i], updatedGeo_ );
      }
    }
    lin->SetName("HeatFluxInt");
    LinearFormContext *ctx = new LinearFormContext( lin );
    ctx->SetEntities( ent[i] );
    ctx->SetFeFunction(myFct);
    assemble_->AddLinearForm(ctx);
    myFct->AddEntityList(ent[i]);
  } // end loop over entities


  //========================================================================================
  // Heat transfer boundary condition
  //========================================================================================
  this->HeatTransferBC();

  //========================================================================================
  // Thermal Radiation boundary condition
  //========================================================================================
  this->ThermalRadiationBC();


  //     // =======================================================================
  //     // Integrators for NonConforming Interfaces
  //     // =======================================================================

  //     // Get index of LAGRANGE_MULT result, just in case... who knows...
  //     UInt lmResultIdx = 0;
  //     for(UInt i=0, n=results_.GetSize(); i<n; i++) {
  //       if(results_[i]->resultType == LAGRANGE_MULT) {
  //         lmResultIdx = i;
  //         break;
  //       }
  //     }
  //     LOG_DBG2(heatcondpde) << "NonMatching: Index of LAGRANGE_MULT result: "
  //                      << lmResultIdx;

  //     for( UInt i = 0; i < ncIFaces_.GetSize(); i++ ) {

  //       // get regionId of Lagrangian surface
  //       StdVector<std::string> keyVec, attrVec, valVec;
  //       std::string secondarySide;
  //       std::string ncIfaceName = ptGrid_->GetRegion().ToString(ncIFaces_[i]);

  //       PtrParamNode ncIfaceListNode;
  //       ncIfaceListNode = domain_->GetParamRoot()->Get("domain")->Get("ncInterfaceList");

  //       secondarySide = ncIfaceListNode->
  //           GetByVal("ncInterface", "name",
  //                    ncIfaceName)->Get("secondarySide")->As<std::string>();

  //       // Part 1: Define integrator M(Psi, Lambda) on
  //       //         non-conforming interface
  //       LOG_DBG2(heatcondpde) << "NonMatching: Defining nonconforming integrator"
  //                         << " for M on interface '"
  //                         << ptGrid_->GetRegion().ToString(ncIFaces_[i]) << "'.";
  //       shared_ptr<ElemList> actNcList( new ElemList(ptGrid_ ) );
  //       actNcList->SetRegion( ncIFaces_[i] );

  //       NonConformingInt * ncInt =
  //         new NonConformingInt( 1, isaxi_ );

  //       NcBiLinFormContext * stiffIntDescr =
  //      	  new NcBiLinFormContext( ncInt , STIFFNESS );

  //       // Force assembling of M(Psi, Lambda)^T
  //       stiffIntDescr->SetCounterPart( true );

  //       stiffIntDescr->SetPtPdes(this, this);
  //       stiffIntDescr->SetResults( results_[0], results_[lmResultIdx],
  //                                  actNcList, actNcList );

  //       assemble_->AddBiLinearForm( stiffIntDescr );


  //       // Part 2: Define integrator D(Psi, Lambda) on
  //       //         Lagrangian surface
  //       LOG_DBG2(heatcondpde) << "NonMatching: Defining mass integrator"
  //                         << " for D on interface '"
  //                         << ptGrid_->GetRegion().ToString(ncIFaces_[i]) << "'.";
  //       shared_ptr<SurfElemList> actSDList( new SurfElemList(ptGrid_ ) );
  //       actSDList->SetRegion( ptGrid_->GetRegion().Parse(secondarySide));

  //       // D(Psi, Lambda) has the form of a standard mass
  //       // integrator with factor 1.0
  //       MassInt * dMatInt = new MassInt( 1.0, 1, isaxi_ );
  //       BiLinFormContext * dMatContext =
  //         new BiLinFormContext( dMatInt, STIFFNESS );

  //       // Force assembling of D(Psi, Lambda)^T
  //       dMatContext->SetCounterPart( true );
  //       dMatContext->SetPtPdes( this, this );
  //       dMatContext->SetResults( results_[0], results_[lmResultIdx],
  //                                actSDList, actSDList );

  //       assemble_->AddBiLinearForm( dMatContext );

  //       // Give result LAGRANGE_MULT to equation numbering class
  //       eqnMap_->AddResult( *results_[lmResultIdx], actSDList );
  //     }


  // TO BE DONE

  //   try {
  //     // iterate over all parameter nodes
  //     for( UInt i = 0; i < rhsValuesNodes.GetSize(); i++ )
  //     {
  //       std::string rhsRegion(rhsValuesNodes[i]->Get("region")->As<std::string>());
  //       rhsValuesNodes[i]->GetValue("isharmonic", isharmonic, ParamNode::EX);
  //       //rhsFileId = rhsValuesNodes[i]->Get("inputId")->As<std::string>();

  //       if ( isharmonic ) {
  //         Info->PrintF( pdename_,
  //             "Use same RHS source vector for region '%s' for all timesteps.\n\n",
  //             rhsRegion.c_str() );
  //       } else {
  //         Info->PrintF( pdename_,
  //             "Use RHS source vector for region '%s' from acoustic results.\n\n",
  //             rhsRegion.c_str() );

  //       }

  //       // get material density
  //       BaseMaterial * actMat = materials_[ ptGrid_->GetRegion().Parse(rhsRegion) ];
  //       actMat->GetScalar(density,DENSITY,Global::REAL);

  //       linAcouPowerSourceInt* sourceRHSInt = new linAcouPowerSourceInt( isaxi_,
  //           isharmonic, rhsFileId, rhsRegion, density );

  //       LinearFormContext* sourceRHSContext = new LinearFormContext( sourceRHSInt );
  //       sourceRHSContext->SetPtPde( this );

  //       shared_ptr<ElemList> rhsElemList( new ElemList(ptGrid_ ) );
  //       rhsElemList->SetRegion( ptGrid_->GetRegion().Parse(rhsRegion) );
  //       sourceRHSContext->SetResult( results_[0], rhsElemList );
  //       assemble_->AddLinearForm( sourceRHSContext );
  //     }
  //   }
  //   catch(Exception & ex)
  //   {
  //     RETHROW_EXCEPTION(ex, "Could not assemble RHS source integrator"
  //         <<" in HeatPDE" );
  //   }

}


void HeatPDE::DefineNcIntegrators() {
  StdVector< NcInterfaceInfo >::iterator ncIt = ncInterfaces_.Begin(),
                                         endIt = ncInterfaces_.End();
  for ( ; ncIt != endIt; ++ncIt ) {
    switch (ncIt->type) {
    case NC_MORTAR:
      if (dim_ == 2)
        DefineMortarCoupling<2,1>(HEAT_TEMPERATURE, *ncIt);
      else
        DefineMortarCoupling<3,1>(HEAT_TEMPERATURE, *ncIt);
      break;
    case NC_NITSCHE:
      if (dim_ == 2)
        DefineNitscheCoupling<2,1>(HEAT_TEMPERATURE, *ncIt );
      else
        DefineNitscheCoupling<3,1>(HEAT_TEMPERATURE, *ncIt );
      break;
    default:
      EXCEPTION("Unknown type of ncInterface");
      break;
    }
  }
}


void HeatPDE::HeatTransferBC(){
  // here we cannot simply call ReadRhsExcitation, because it only looks for "value" and "name"
  // but we have two values (heat transfer coefficient and bulk temperature), so we have to do
  // the whole thing here

  LOG_TRACE(heatcondpde) << "Defining heat transfer BC for thermal PDE";

  // Get FESpace and FeFunction
  shared_ptr<BaseFeFunction> myFct = feFunctions_[HEAT_TEMPERATURE];
  shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();

  StdVector<shared_ptr<EntityList> > ent;
  StdVector<PtrCoefFct > coef;
  LinearForm * lin = NULL;
  BiLinearForm * heatTransInt = NULL;

  bool coefUpdateGeo = true;
  std::string elemName = "heatTransport";
  ParamNodeList elems = myParam_->Get( "bcsAndLoads" )->GetList( "heatTransport" );
  ent.Resize(elems.GetSize());
  coef.Resize(elems.GetSize());
  std::string entName;

  for( UInt i = 0; i < elems.GetSize(); ++i ) {
        PtrParamNode xml = elems[i];
        try {
            // get entity list, depending on type
            entName = xml->Get("name")->As<std::string>();
            // check if we only have surface elements in the given region
            EntityList::ListType listType = EntityList::ELEM_LIST;
            if( ptGrid_->GetEntityDim( entName ) == ptGrid_->GetDim() - 1) {
              listType = EntityList::SURF_ELEM_LIST;
            }else{
              EXCEPTION("Heat transfer BC can only be applied on surface regions");
            }
            switch( ptGrid_->GetEntityType(entName) ) {
              case EntityList::NAMED_NODES:
                ent[i] = ptGrid_->GetEntityList( EntityList::NODE_LIST, entName);
                break;
              case EntityList::REGION:
              case EntityList::NAMED_ELEMS:
                ent[i] = ptGrid_->GetEntityList( listType, entName );
                break;
              case EntityList::NO_TYPE:
                EXCEPTION("No entities with name '" << entName << "' known");
                break;
            }
        } catch (Exception& e) {
          RETHROW_EXCEPTION(e, pdename_ << ": Could not read definition for '" << elemName
                            << "' on entities '" << entName <<"'");
        }

        // check just defined type of entitylist
        if (ent[i]->GetType() == EntityList::NODE_LIST) {
          EXCEPTION("Heat flux must be defined on elements")
        }


        std::string T_b = xml->Get("bulkTemperature")->As<std::string>();
        std::string alpha = xml->Get("heatTransferCoefficient")->As<std::string>();
        PtrCoefFct alphaCoef = CoefFunction::Generate( mp_, Global::REAL, alpha);
        PtrCoefFct bulkCoef = CoefFunction::Generate( mp_, Global::REAL, T_b);
        std::string volRegName = xml->Get("volumeRegion")->As<std::string>();
        RegionIdType volRegion = ptGrid_->GetRegion().Parse(volRegName);

        // Get the reference Temperatur to make the system symmetric
        PtrCoefFct refTemp = nullptr;
        if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
          refTemp = materials_[volRegion]->GetScalCoefFnc(HEAT_REF_TEMPERATURE, Global::REAL);
        }

        // check if volume region has mapping or PML...should already be catched
        if( dampingList_[volRegion] == MAPPING ){
          EXCEPTION("Infinite mapping and heat transfer BC not yet verified and tested, therefore it's disabled");
        }


        //========================================================================================
        // First part of heat transfer boundary condition \int_{\Omega} \alpha T' T dS
        //========================================================================================
        // factor for mass matrix: \frac{\alpha}{\lambda}
        // \alpha... heat transfer coefficient
        PtrCoefFct factor1 = CoefFunction::Generate( mp_, Global::REAL, "1.0");

        PtrCoefFct coeffMass = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_, factor1, alphaCoef, CoefXpr::OP_MULT ) );
        if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
          coeffMass = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp( mp_, coeffMass, refTemp, CoefXpr::OP_DIV));
        }

        if (coeffMass->IsComplex()) {
          EXCEPTION("HeatPDE: Complex coefficient for mass matrix not implemented");
        } else {
            if( dim_ == 2 ) {
              heatTransInt = new BBInt<>(new IdentityOperator<FeH1,2,1>(), coeffMass, 1.0, updatedGeo_ );
            } else {
              heatTransInt = new BBInt<>(new IdentityOperator<FeH1,3,1>(), coeffMass, 1.0, updatedGeo_ );
            }
        }
        heatTransInt->SetName("heatTransportPart1Integrator");
        heatTransInt->SetFeSpace(feFunctions_[HEAT_TEMPERATURE]->GetFeSpace());

        BiLinFormContext *heatTransContext = new BiLinFormContext(heatTransInt, STIFFNESS);
        heatTransContext->SetEntities( ent[i], ent[i] );
        heatTransContext->SetFeFunctions( myFct , myFct);
        feFunctions_[HEAT_TEMPERATURE]->AddEntityList( ent[i] );
        assemble_->AddBiLinearForm( heatTransContext );


        // =======================================================================================
        //  Second Part of heat transfer BC \int_{\Omega} \alpha T' T_bulk dS
        // =======================================================================================
        PtrCoefFct factor2 = CoefFunction::Generate( mp_, Global::REAL, "1.0");
        PtrCoefFct bulkTimesAlpha = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, bulkCoef, alphaCoef, CoefXpr::OP_MULT ));
        PtrCoefFct coeffRHS = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_, factor2, bulkTimesAlpha, CoefXpr::OP_MULT ) );
        if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
          coeffRHS = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp( mp_, coeffRHS, refTemp, CoefXpr::OP_DIV));
        }
        if(isComplex_) {
          lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1>(), Complex(1.0), coeffRHS, coefUpdateGeo);
        } else  {
          lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1>(), 1.0, coeffRHS, coefUpdateGeo);
        }
        // obtain coordinate system and set it at coefficient function
        std::string coordSysId = "default";
        xml->GetValue("coordSysId", coordSysId, ParamNode::PASS);
        if( coordSysId != "default" ) {
          coeffRHS->SetCoordinateSystem( domain_->GetCoordSystem(coordSysId) );
        }

        lin->SetName("HeatTransPart2Int");
        lin->SetFeSpace(feFunctions_[HEAT_TEMPERATURE]->GetFeSpace());

        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);

  }
}

void HeatPDE::ThermalRadiationBC(){
  // incorporate thermal radiation bc
  LOG_TRACE(heatcondpde) << "Defining thermal radiation BC for thermal PDE";

  // Get FESpace and FeFunction
  shared_ptr<BaseFeFunction> myFct = feFunctions_[HEAT_TEMPERATURE];
  shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();

  StdVector<shared_ptr<EntityList> > ent;
  StdVector<PtrCoefFct > coef;
  LinearForm * lin = NULL;
  BiLinearForm * thermRadInt = NULL;

  bool coefUpdateGeo = true;
  std::string elemName = "thermalRadiation";
  ParamNodeList elems = myParam_->Get( "bcsAndLoads" )->GetList( "thermalRadiation" );
  ent.Resize(elems.GetSize());
  coef.Resize(elems.GetSize());
  std::string entName;

  for (UInt i = 0; i < elems.GetSize(); ++i) {
	  PtrParamNode xml = elems[i];
	  try {
		  // get entity list, depending on type
		  entName = xml->Get("name")->As<std::string>();
		  // check if we only have surface elements in the given region
		  EntityList::ListType listType = EntityList::ELEM_LIST;
		  if( ptGrid_->GetEntityDim( entName ) == ptGrid_->GetDim() - 1) {
			listType = EntityList::SURF_ELEM_LIST;
		  } else {
			EXCEPTION("Thermal Radiation BC can only be applied on surface regions");
		  }
		  switch( ptGrid_->GetEntityType(entName) ) {
		    case EntityList::NAMED_NODES:
			  ent[i] = ptGrid_->GetEntityList( EntityList::NODE_LIST, entName);
			  break;
			case EntityList::REGION:
			case EntityList::NAMED_ELEMS:
			  ent[i] = ptGrid_->GetEntityList( listType, entName );
			  break;
			case EntityList::NO_TYPE:
			  EXCEPTION("No entities with name '" << entName << "' known");
			  break;
		  }
	  } catch (Exception &e) {
		  RETHROW_EXCEPTION(e, pdename_ << ": Could not read definition for '" << elemName
		                              << "' on entities '" << entName <<"'");
	  }

	  // check just defined type of entitylist
	  if (ent[i]->GetType() == EntityList::NODE_LIST) {
		EXCEPTION("Thermal Radiation BC must be defined on elements")
	  }

	  std::string T0 = xml->Get("bulkTemperature")->As<std::string>();
	  std::string epsilon = xml->Get("emissivity")->As<std::string>();
	  std::string sigma = "5.670374419e-8"; //stefan-boltzmann-constant

	  PtrCoefFct epsilonCoef = CoefFunction::Generate(mp_, Global::REAL, epsilon);
	  PtrCoefFct T0Coef = CoefFunction::Generate(mp_, Global::REAL, T0);
	  PtrCoefFct sigmaCoef = CoefFunction::Generate(mp_, Global::REAL, sigma);

    if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
      EXCEPTION("Thermal Radiation is not implemented for the symmetric formulation of the LinFLow-Heat-Coupling.")
    }

	  //========================================================================================
	  // First part of thermal radiation boundary condition  4 * \epsilon \sigma * (T_{k-1})^3 \int_{\Gamma} T' T_k dS
	  //========================================================================================
	  // \epsilon ... emissivity
	  // \sigma ... stefan-boltzmann constant

	  // get current temperature T_{k-1}
	  PtrCoefFct currT = this->GetCoefFct(HEAT_TEMPERATURE);

	  currT->SetSolDependent();

	  PtrCoefFct factor3 = CoefFunction::Generate(mp_, Global::REAL, "3.0");
	  PtrCoefFct factor4 = CoefFunction::Generate(mp_, Global::REAL, "4.0");
	  // T_{k-1}^3
	  PtrCoefFct TPow3 = CoefFunction::Generate(mp_, Global::REAL,
	  				CoefXprBinOp(mp_, currT, factor3, CoefXpr::OP_POW));
	  // \epsilon * \sigma
	  PtrCoefFct SigmaTimesEpsilon = CoefFunction::Generate(mp_, Global::REAL,
	  				CoefXprBinOp(mp_, sigmaCoef, epsilonCoef, CoefXpr::OP_MULT));
	  // 4 * \epsilon * \sigma
	  PtrCoefFct FourTimesSigmaTimesEpsilon = CoefFunction::Generate(mp_, Global::REAL,
	  						CoefXprBinOp(mp_, SigmaTimesEpsilon, factor4, CoefXpr::OP_MULT));
	  // coefficient function of LHS : 4 * \epsilon \sigma * (T_{k-1})^3
	  PtrCoefFct coeffLHS = CoefFunction::Generate(mp_, Global::REAL,
	  				CoefXprBinOp(mp_, TPow3, FourTimesSigmaTimesEpsilon, CoefXpr::OP_MULT));

	  if (dim_ == 2) {
		thermRadInt = new BBInt<>(new IdentityOperator<FeH1, 2, 1>(),
		coeffLHS, 1.0, updatedGeo_);
	  } else {
		thermRadInt = new BBInt<>(new IdentityOperator<FeH1, 3, 1>(),
		coeffLHS, 1.0, updatedGeo_);
	  }
	  thermRadInt->SetName("ThermalRadiationLHSNonLinInt");
	  thermRadInt->SetFeSpace(feFunctions_[HEAT_TEMPERATURE]->GetFeSpace());

	  BiLinFormContext *thermRadContext = new BiLinFormContext(thermRadInt,
	  				STIFFNESS);
	  thermRadContext->SetEntities(ent[i], ent[i]);
	  thermRadContext->SetFeFunctions(myFct, myFct);
	  feFunctions_[HEAT_TEMPERATURE]->AddEntityList(ent[i]);
	  assemble_->AddBiLinearForm(thermRadContext);

	  // =======================================================================================
	  //  Second Part of heat transfer BC \int_{\Gamma} \epsilon * \sigma  (3 *T_{k-1}^4 + T0^4) T' dS
	  // =======================================================================================

	  // T_{k-1}^4
	  PtrCoefFct TPow4 = CoefFunction::Generate(mp_, Global::REAL,
				CoefXprBinOp(mp_, currT, factor4, CoefXpr::OP_POW));
	  // T0^4
	  PtrCoefFct T0Pow4 = CoefFunction::Generate(mp_, Global::REAL,
				CoefXprBinOp(mp_, T0Coef, factor4, CoefXpr::OP_POW));
	  // 3 * T_{k-1}^4
	  PtrCoefFct ThreeTimesTPow4 = CoefFunction::Generate(mp_, Global::REAL,
				CoefXprBinOp(mp_, TPow4, factor3, CoefXpr::OP_MULT));
	  // 3 * T_{k-1}^4 + T0^4
	  PtrCoefFct ThreeTimesTPow4PlusT0Pow4 = CoefFunction::Generate(mp_, Global::REAL,
				CoefXprBinOp(mp_, ThreeTimesTPow4, T0Pow4, CoefXpr::OP_ADD));
	  // coefficient function of RHS : \epsilon * \sigma  (3 *T_{k-1}^4 + T0^4)
	  PtrCoefFct coeffRHS = CoefFunction::Generate(mp_, Global::REAL,
				CoefXprBinOp(mp_, SigmaTimesEpsilon, ThreeTimesTPow4PlusT0Pow4, CoefXpr::OP_MULT));
	  if (isComplex_) {
	    lin = new BUIntegrator<Complex>(new IdentityOperator<FeH1>(),
					Complex(1.0), coeffRHS, coefUpdateGeo);
	  } else {
	  	lin = new BUIntegrator<Double>(new IdentityOperator<FeH1>(), 1.0,
					coeffRHS, coefUpdateGeo);
	  }
	  // obtain coordinate system and set it at coefficient function
	  std::string coordSysId = "default";
	  xml->GetValue("coordSysId", coordSysId, ParamNode::PASS);
	  if (coordSysId != "default") {
	  	coeffRHS->SetCoordinateSystem(domain_->GetCoordSystem(coordSysId));
	  }
      lin->SetName("ThermalRadiationRHSNonLinInt");
	  lin->SetFeSpace(feFunctions_[HEAT_TEMPERATURE]->GetFeSpace());
	  lin->SetSolDependent();
	  LinearFormContext *ctx = new LinearFormContext(lin);
	  ctx->SetEntities(ent[i]);
	  ctx->SetFeFunction(myFct);
	  assemble_->AddLinearForm(ctx);
  }
}


void HeatPDE::DefineRhsLoadIntegrators() {
  
  LOG_TRACE(heatcondpde) << "Defining rhs load integrators for thermal PDE";

  // Get FESpace and FeFunction of electric potential
  shared_ptr<BaseFeFunction> myFct = feFunctions_[HEAT_TEMPERATURE];
  shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();

  StdVector<shared_ptr<EntityList> > ent;
  StdVector<PtrCoefFct > coef;
  LinearForm * lin = NULL;
  StdVector<std::string> dofNames;
  StdVector<std::string> volumeRegions;
    

  bool coefUpdateGeo = true;
  // =====================
  //  HEAT SOURCE DENSITY
  // =====================
  LOG_DBG(heatcondpde) << "Reading heat source density";

  ReadRhsExcitation( "heatSourceDensity", dofNames, ResultInfo::SCALAR, isComplex_, ent, coef, coefUpdateGeo );
  if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
    // Get the volumeRegions for the heatSourceDensity, because the reference temperature of the volumeRegion is needed to make the system symmetric
    ReadVolumeRegions("heatSourceDensity", volumeRegions);
    if( ent.GetSize() != volumeRegions.GetSize()){
      EXCEPTION("Define a volumeRegion for the heatSourceDensity or set the attribute symmetric=\"false\" in <couplingList> -> <direct> -> <linFlowHeatDirect>");
    }
  }
  
  for( UInt i = 0; i < ent.GetSize(); ++i ) {
    // check type of entitylist
    if (ent[i]->GetType() == EntityList::NODE_LIST) {
      EXCEPTION("Heat source density must be defined on elements")
    }
    EntityIterator it = ent[i]->GetIterator();
    it.Begin();

    if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
      if(volumeRegions[i] == ""){
        EXCEPTION("Define a volumeRegion for the heatSourceDensity or set the attribute symmetric=\"false\" in <couplingList> -> <direct> -> <linFlowHeatDirect>");
      }
      RegionIdType volRegion = ptGrid_->GetRegion().Parse(volumeRegions[i]);
      PtrCoefFct refTemp = materials_[volRegion]->GetScalCoefFnc(HEAT_REF_TEMPERATURE, Global::REAL);
      coef[i] = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp( mp_, coef[i], refTemp, CoefXpr::OP_DIV));
    }

    if(isComplex_) {
      //pure complex case (harmonic simulation)
      lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1>(), Complex(1.0), coef[i], coefUpdateGeo);
    } else  {
    	if(coef[i]->IsComplex()){
    		// rhs excitation comes from a harmonic simulation (actually a real result with zero imaginary part)
    		lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1>(), Complex(1.0), coef[i], coefUpdateGeo, true,
    		      	  	  	  	  	  	  	  	   coef[i]->IsComplex() );
    	}else{
    		lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1>(), 1.0, coef[i], coefUpdateGeo, true,
    		      	  	  	  	  	  	  	  	   coef[i]->IsComplex() );
    	}

    }
    lin->SetName("HeatSourceDensityInt");
    LinearFormContext *ctx = new LinearFormContext( lin );
    ctx->SetEntities( ent[i] );
    ctx->SetFeFunction(myFct);
    assemble_->AddLinearForm(ctx);
  } // end loop over entities

  // ========================
  //  HEAT SOURCE(nodal)
  // ========================
  LOG_DBG(heatcondpde) << "Reading heat source values";

  coefUpdateGeo = false;

  ReadRhsExcitation( "heatSource", dofNames, ResultInfo::SCALAR, isComplex_, ent, coef, coefUpdateGeo);
  if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
    // Get the volumeRegions for the heatSource, because the reference temperature of the volumeRegion is needed to make the system symmetric
    ReadVolumeRegions("heatSource", volumeRegions);
    if( ent.GetSize() != volumeRegions.GetSize()){
      EXCEPTION("Define a volumeRegion for the heatSource or set the attribute symmetric=\"false\" in <couplingList> -> <direct> -> <linFlowHeatDirect>");
    }
  }

  for( UInt i = 0; i < ent.GetSize(); ++i ) {
    // assume that we have elem list due to specification of a region instead of named nodes in xml file
    if (ent[i]->GetType() != EntityList::NODE_LIST && ent[i]->GetType() != EntityList::ELEM_LIST) {
      EXCEPTION("Heat source must be defined on nodes!")
    }

    unsigned int numNodes = ent[i]->GetSize();
    const std::string numNodes_str = std::to_string(numNodes);
    if(numNodes > 1 && coef[i]->DoNormalize())
      coef[i] = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, coef[i], numNodes_str, CoefXpr::OP_DIV));

    if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
      if(volumeRegions[i] == ""){
        EXCEPTION("Define a volumeRegion for the heatSource or set the attribute symmetric=\"false\" in <couplingList> -> <direct> -> <linFlowHeatDirect>");
      }
      RegionIdType volRegion = ptGrid_->GetRegion().Parse(volumeRegions[i]);
      PtrCoefFct refTemp = materials_[volRegion]->GetScalCoefFnc(HEAT_REF_TEMPERATURE, Global::REAL);
      coef[i] = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp( mp_, coef[i], refTemp, CoefXpr::OP_DIV));
    }

    lin = new SingleEntryInt(coef[i]);
    lin->SetName("NodalHeatInt");
    if(nonLinTypes_.find("source") != nonLinTypes_.end()) {
      lin->SetName("NonLinNodalHeatInt");
      lin->SetSolDependent();
    }
    LinearFormContext *ctx = new LinearFormContext( lin );
    ctx->SetEntities( ent[i] );
    ctx->SetFeFunction(myFct);
    assemble_->AddLinearForm(ctx);
  } // end loop over entities

  // ========================
  //  DESIGN DEPENDENT HEAT SOURCE
  // ========================
  LOG_DBG(heatcondpde) << "Reading heat source values (design dependent)";

  ReadRhsExcitation( "designDependentHeatSource", dofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo);
  if(GetParamNode()->Has("bcsAndLoads/designDependentHeatSource"))
    interfaceDrivenHeatSource_ = true;
  for( UInt i = 0; i < ent.GetSize(); ++i ) {
    // assume that we have elem list due to specification of a region instead of named nodes in xml file
    if(ent[i]->GetType() != EntityList::NODE_LIST && ent[i]->GetType() != EntityList::ELEM_LIST) {
      EXCEPTION("Design dependent heat source must be defined on nodes!")
    }

    if(domain->HasDesign())
    {
      CoefFunctionOpt* tmpFnc = new CoefFunctionOpt(domain->GetDesign(), coef[i], NO_MATERIAL, this); // takes double and complex
      coef[i].reset(tmpFnc);
    }

    if (isLinFlowPDECoupled_ && isCouplingFormulationSymmetric_) {
      EXCEPTION("Design dependent heat sources are not implememtated for a symmetric coupling with the LinFlowPDE.");
    }

    lin = new SingleEntryInt(coef[i]);
    lin->SetName("DesignDepHeatInt");

    LinearFormContext *ctx = new LinearFormContext( lin );
    ctx->SetEntities(ent[i]);
    ctx->SetFeFunction(myFct);
    assemble_->AddLinearForm(ctx);
  }

}

// similar to MechPDE, except that we have only 3 "test strains" in 3D (2 in 2D)
void HeatPDE::DefineTestStrainIntegrator(const TestStrain test, StdVector<LinearFormContext*>* linForms)
{
  LOG_DBG(heatcondpde) << "DTTGI: test=" << test << " lf=" << (linForms ? linForms->ToString() : "null");

  shared_ptr<BaseFeFunction> myFct = feFunctions_[HEAT_TEMPERATURE];

  // to generate the vector values coef function for the test strain we need scalar const coef functions for zero and one
  PtrCoefFct one  = CoefFunction::Generate(mp_, Global::REAL, "1.0");
  PtrCoefFct zero = CoefFunction::Generate(mp_, Global::REAL, "0.0");

  StdVector<PtrCoefFct> tempGrad(dim_ == 2 ? 2 : 3);
  tempGrad.Init(zero);
  tempGrad[test] = one;
  LOG_DBG(heatcondpde) << "DTTGI: idx=" << test << " -> one";

  std::map<RegionIdType, BaseMaterial*>::iterator it;
  for(it = materials_.begin(); it != materials_.end(); it++)
  {
    // Set current region and material
    RegionIdType actRegion = it->first;

    shared_ptr<EntityList> actSDList = ptGrid_->GetEntityList( EntityList::ELEM_LIST, ptGrid_->GetRegionName(actRegion));

    std::multimap<RegionIdType, BaseBDBInt*>::iterator stiffIt = bdbInts_.begin();
    BaseBDBInt* bdb = nullptr;
    if( bdbInts_.count(actRegion)==1 ) {
      for(; stiffIt != bdbInts_.end(); ++stiffIt ) {
        RegionIdType region = stiffIt->first;
        if(region==actRegion) {
          bdb = stiffIt->second;
          break;
        }
      }
    } else {
      EXCEPTION("Implementation error: Multiple bdbInts_ defined, can't return single coefFunction. You probably used bdbInts.insert() multiple times per region and did not use SetIntegratorName().")
    }

    PtrCoefFct curCoef = bdb->GetCoef();
    PtrCoefFct ttg = CoefFunction::Generate(mp_, Global::REAL, tempGrad);

    LinearForm* lin = NULL;
    if(dim_ == 3)
      lin = new BDUIntegrator<GradientOperator<FeH1,3>, double>(1.0, ttg, curCoef, false); // no updateGeo
    else
      lin = new BDUIntegrator<GradientOperator<FeH1,2>, double>(1.0, ttg, curCoef, false); // no updateGeo

    lin->SetName("TestTempGradInt");
    LinearFormContext* ctx = new LinearFormContext(lin);
    ctx->SetEntities(actSDList);
    ctx->SetFeFunction(myFct);

    if(linForms != NULL)
      linForms->Push_back(ctx);
    else
      assemble_->AddLinearForm(ctx);
  }
}


void HeatPDE::DefineSolveStep() {

  solveStep_ = new StdSolveStep(*this);

}

// ======================================================
// TIME STEPPING SECTION
// ======================================================

void HeatPDE::InitTimeStepping() {

  // Until now no effective mass formulation in the trapezoidal
  //  integration scheme is implemented!
  Double gamma = 0.5; 
  GLMScheme * scheme = new Trapezoidal(gamma);

  if ( nonLinTotalFormulation_ ) {
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0)); //, nlType) );
    feFunctions_[HEAT_TEMPERATURE]->SetTimeScheme(myScheme);
  }
  else {
    TimeSchemeGLM::NonLinType nlType = (nonLin_)? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType) );
    feFunctions_[HEAT_TEMPERATURE]->SetTimeScheme(myScheme);
  }
}


void HeatPDE::DefinePrimaryResults() {

  // === TEMPERATURE ===
  shared_ptr<ResultInfo> res1( new ResultInfo);
  res1->resultType = HEAT_TEMPERATURE;
  res1->dofNames = "";
  res1->unit = MapSolTypeToUnit(HEAT_TEMPERATURE);
  res1->definedOn = ResultInfo::NODE;
  res1->entryType = ResultInfo::SCALAR;
  feFunctions_[HEAT_TEMPERATURE]->SetResultInfo(res1);
  results_.Push_back( res1 );
  availResults_.insert( res1 );
  res1->SetFeFunction(feFunctions_[HEAT_TEMPERATURE]);
  DefineFieldResult( feFunctions_[HEAT_TEMPERATURE], res1 );

  shared_ptr<ResultInfo> res2( new ResultInfo);
  res2->resultType = HEAT_MEAN_TEMPERATURE;
  res2->dofNames = "";
  res2->unit = MapSolTypeToUnit(HEAT_MEAN_TEMPERATURE);
  res2->definedOn = ResultInfo::NODE;
  res2->entryType = ResultInfo::SCALAR;
  results_.Push_back( res2 );
  availResults_.insert( res2 );
  PtrCoefFct tmpReal = NULL;
  tmpReal.reset(new CoefFunctionComplexToReal<Double>(feFunctions_[HEAT_TEMPERATURE], ptGrid_));
  DefineFieldResult( tmpReal, res2 );

  //  DefineFieldResult( convecVelCoef_, velocity );
  //  res2->SetFeFunction(feFunctions_[HEAT_TEMPERATURE]);
  //  DefineFieldResult( feFunctions_[HEAT_TEMPERATURE], res2 );



  // -----------------------------------
  //  Define xml-names of Dirichlet BCs
  // -----------------------------------
  idbcSolNameMap_[HEAT_TEMPERATURE] = "temperature";
  

  //creates the convective velocity
  StdVector<std::string> vecDofNames;
  if( ptGrid_->GetDim() == 3 ) {
    vecDofNames = "x", "y", "z";
  } else {
    if( ptGrid_->IsAxi() ) {
      vecDofNames = "r", "z";
    } else {
      vecDofNames = "x", "y";
    }
  }

  //// === VELOCITY ===
  shared_ptr<ResultInfo> velocity( new ResultInfo);
  velocity->resultType = MEAN_FLUIDMECH_VELOCITY;
  velocity->dofNames = vecDofNames;
  velocity->unit = MapSolTypeToUnit(MEAN_FLUIDMECH_VELOCITY);

  velocity->definedOn = ResultInfo::NODE;
  velocity->entryType = ResultInfo::VECTOR;

  convecVelCoef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
  DefineFieldResult( convecVelCoef_, velocity );

  results_.Push_back( velocity );
  availResults_.insert( velocity );

}

void HeatPDE::DefinePostProcResults() {
  shared_ptr<BaseFeFunction> feFct = feFunctions_[HEAT_TEMPERATURE];

  if ( analysistype_ != STATIC ) {
    // === TEMPERATURE D1===
    // first time derivative of Temperature 
    // ${\frac{\partial T} {\partial t}}$ or $j \omega T$
    shared_ptr<ResultInfo> heatD1( new ResultInfo);
    heatD1->resultType = HEAT_TEMPERATURE_D1;

    heatD1->dofNames = "";
    heatD1->unit = MapSolTypeToUnit(HEAT_TEMPERATURE_D1);
    heatD1->definedOn = ResultInfo::NODE;
    heatD1->entryType = ResultInfo::SCALAR;
    availResults_.insert( heatD1 );
    DefineTimeDerivResult( HEAT_TEMPERATURE_D1, 1, HEAT_TEMPERATURE );
  }
  
  // === TEMPERATURE RHS ===
  shared_ptr<ResultInfo> rhs ( new ResultInfo );
  rhs->resultType = HEAT_RHS_LOAD;
  rhs->dofNames = "";
  rhs->unit = "J";
  rhs->definedOn = ResultInfo::NODE;
  rhs->entryType = ResultInfo::SCALAR;
  rhsFeFunctions_[HEAT_TEMPERATURE]->SetResultInfo(rhs);
  DefineFieldResult( rhsFeFunctions_[HEAT_TEMPERATURE], rhs );

  // === HEAT FLUX DENSITY ===
  // Heat Flux Density $\bm{q}  = -k \nabla T$
  shared_ptr<ResultInfo> fluxDens ( new ResultInfo );
  fluxDens->resultType = HEAT_FLUX_DENSITY;
  fluxDens->SetVectorDOFs(dim_, isaxi_);
  fluxDens->unit = MapSolTypeToUnit(HEAT_FLUX_DENSITY);
  fluxDens->definedOn = ResultInfo::ELEMENT;
  fluxDens->entryType = ResultInfo::VECTOR;
  shared_ptr<CoefFunctionFormBased> fluxDensFunc;
  if( isComplex_ ) {
    fluxDensFunc.reset(new CoefFunctionFlux<Complex>(feFct, fluxDens, Complex(-1.0)));
  } else {
    fluxDensFunc.reset(new CoefFunctionFlux<Double>(feFct, fluxDens, -1.0));
  }
  DefineFieldResult( fluxDensFunc, fluxDens );
  stiffFormCoefs_.insert(fluxDensFunc);

  // === HEAT FLUX INTENSITY ===
  // Heat Flux Density $q_n = \bm{q} \cdot \bm{n}$
  shared_ptr<ResultInfo> fluxNormal ( new ResultInfo );
  fluxNormal->resultType = HEAT_FLUX_INTENSITY;
  fluxNormal->unit = MapSolTypeToUnit(HEAT_FLUX_INTENSITY);
  fluxNormal->dofNames = "";
  fluxNormal->definedOn = ResultInfo::SURF_ELEM;
  fluxNormal->entryType = ResultInfo::SCALAR;
  shared_ptr<CoefFunctionSurf> fluxNormalFunc;
  fluxNormalFunc.reset(new CoefFunctionSurf(true, 1.0, fluxNormal));
  DefineFieldResult( fluxNormalFunc, fluxNormal );
  surfCoefFcts_[fluxNormalFunc] = fluxDensFunc;

  // === HEAT FLUX ===
  // Heat Flux $\dot{Q} = \int_{\Gamma} \bm{q} \cdot \bm{n} \mathrm{d} \Gamma$
  shared_ptr<ResultInfo> flux ( new ResultInfo );
  flux->resultType = HEAT_FLUX;
  flux->unit = MapSolTypeToUnit(HEAT_FLUX);
  flux->dofNames = "";
  flux->definedOn = ResultInfo::SURF_REGION;
  flux->entryType = ResultInfo::SCALAR;
  shared_ptr<ResultFunctor> fluxFunc;
  availResults_.insert( flux );
  if( isComplex_ ) {
    fluxFunc.reset(new ResultFunctorIntegrate<Complex>(fluxNormalFunc, feFct, flux) );
  } else {
    fluxFunc.reset(new ResultFunctorIntegrate<Double>(fluxNormalFunc, feFct, flux) );
  }
  resultFunctors_[HEAT_FLUX] = fluxFunc;

  // === HEAT_CONDUCTIVITY_TENSOR_HOM ====
  // cannot call it HEAT_CONDUCTIVITY_TENSOR because there exists already a material type with this name
  shared_ptr<ResultInfo> conduct_tensor(new ResultInfo);
  conduct_tensor->resultType = HEAT_CONDUCTIVITY_TENSOR_HOM;
  conduct_tensor->dofNames = "e11", "e12", "e13", "e22", "e23", "e33";
  conduct_tensor->unit = MapSolTypeToUnit(HEAT_CONDUCTIVITY_TENSOR_HOM);
  conduct_tensor->entryType = ResultInfo::TENSOR;
  conduct_tensor->definedOn = ResultInfo::ELEMENT;
  shared_ptr<CoefFunctionFormBased> conduct_coef;
  conduct_coef.reset(new CoefFunctionHomogenization<double, App::HEAT>(feFct));
  DefineFieldResult(conduct_coef, conduct_tensor);
  stiffFormCoefs_.insert(conduct_coef); // will define the forms

  // optimization results are provided in DesignSpace::ExtractResults()
  // copied from MechPDE
  // === MECH_PSEUDO_DENISTY ===
  shared_ptr<ResultInfo> mpd(new ResultInfo);
  mpd->resultType = MECH_PSEUDO_DENSITY;
  mpd->entryType = ResultInfo::SCALAR;
  mpd->definedOn = ResultInfo::ELEMENT;
  mpd->dofNames = MapSolTypeToUnit(MECH_PSEUDO_DENSITY);
  mpd->fromOptimization = true;
  DefineFieldResult(shared_ptr<FeFunction<double> >(new FeFunction<double>(NULL)), mpd); // the fe-function is only a dummy

  // === PHYSICAL_PSEUDO_DENISTY ===
  shared_ptr<ResultInfo> ppd(new ResultInfo);
  ppd->resultType = PHYSICAL_PSEUDO_DENSITY;
  ppd->entryType = ResultInfo::SCALAR;
  ppd->definedOn = ResultInfo::ELEMENT;
  ppd->dofNames = MapSolTypeToUnit(PHYSICAL_PSEUDO_DENSITY);
  ppd->fromOptimization = true;
  DefineFieldResult(shared_ptr<FeFunction<double> >(new FeFunction<double>(NULL)), ppd);
}

} // end of namespace CoupledField
