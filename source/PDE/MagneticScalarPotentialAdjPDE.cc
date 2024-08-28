#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <set>

#include "MagneticScalarPotentialAdjPDE.hh"

#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionMapping.hh"
#include "Domain/CoefFunction/CoefFunctionDiagTensorFromScalar.hh"
#include "Utils/StdVector.hh"
#include "Driver/Assemble.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/FeFunctions.hh"
#include "Domain/CoefFunction/CoefXpr.hh"

// new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"
#include "Driver/SolveSteps/SolveStepEB.hh"
#include "CoupledPDE/IterCoupledPDE.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"

namespace CoupledField
{

  DEFINE_LOG(magscaladjpde, "magneticScalarPotentialAdjPDE")

  // ***************
  //   Constructor
  // ***************
  MagneticScalarPotentialAdjPDE::MagneticScalarPotentialAdjPDE(Grid *aptgrid, PtrParamNode paramNode,
                                                         PtrParamNode infoNode,
                                                         shared_ptr<SimState> simState, Domain *domain)
      : SinglePDE(aptgrid, paramNode, infoNode, simState, domain)
  {
    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_ = "magneticScalarPotentialAdjPDE";
    pdematerialclass_ = ELECTROMAGNETIC;

    nonLin_ = false;
    nonLinMaterial_ = false;

    //! Always use updated Lagrangian formulation
    updatedGeo_ = true;

    // this is the default value
    modelName_ = "nonlinearCurve";

    // Be aware that the hysteresis via CoefFunctionMaterialModel does not use the isHysteresis_ flag!
    matModelCoef_.reset(new CoefFunctionMaterialModel<Complex>());
    // init the nlFluxCoef
    nlFluxCoef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_, true));
  }

  void MagneticScalarPotentialAdjPDE::DefineIntegrators()
  {

    RegionIdType actRegion;
    BaseMaterial *actSDMat = NULL;

    perm_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, false));

    // get FEFunctions and space
    shared_ptr<BaseFeFunction> feFunc_reduced = feFunctions_[MAG_POTENTIAL_ADJ];
    shared_ptr<FeSpace> feSpace_reduced = feFunc_reduced->GetFeSpace();

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial *>::iterator it;

    //coefFunction holding the magnetic intensity H (result of previous sequence step!!)
    PtrCoefFct magFieldCoef = this->GetCoefFct(MAG_FIELD_INTENSITY);

    // Init material model for hysteretic transient analysis
    if (((analysistype_ == STATIC) || (analysistype_ == TRANSIENT)) && nonLin_ && (modelName_ != "nonlinearCurve"))
    {
      matModelCoef_->Init(magFieldCoef, modelName_, dim_);
    }

    // currently we are only allowed to have one hysteresis region
    bool moreThan1HystRegion = false;

    // iterate over the region (or materials)
    for (UInt iRegion = 0; iRegion < regions_.GetSize(); iRegion++) {
      // set current region and material
      actRegion = regions_[iRegion];
      actSDMat = materials_[actRegion];
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];

      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList(new ElemList(ptGrid_)); 
      actSDList->SetRegion(actRegion);

      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region", "name", regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();

      // --- Set the approximation for the current region ---
      feSpace_reduced->SetRegionApproximation(actRegion, polyId, integId);
      feFunc_reduced->AddEntityList(actSDList);

      // ====================================================================
      // stiffness integrator - NONLINEAR
      // ====================================================================
      BaseBDBInt *stiffInt = NULL;
      PtrCoefFct muNL = NULL;
      PtrCoefFct flux = NULL;
      if (nonLinTypes.Find(PERMEABILITY) != -1) {
        // ====================================================================
        // stiffness integrator - NONLINEAR
        // ====================================================================

        // hysteretic - Energy Based
        if (modelName_ == "JilesAthertonModel") {
          EXCEPTION("Jiles-Atherton model not implemented for MagneticScalarPotentialPDE");
        }
        else if (modelName_ == "EBHysteresisModel") {
          if(moreThan1HystRegion){
            EXCEPTION("Currently only ONE hysteretic region is allowed!");
          }
          moreThan1HystRegion = true;
          
          std::map<std::string, double> ParameterMap;
          actSDMat->GetScalar(ParameterMap["Ps"], MAG_PS_EB, Global::REAL);
          actSDMat->GetScalar(ParameterMap["A"], MAG_A_EB, Global::REAL);
          actSDMat->GetScalar(ParameterMap["mu0"], MAG_MU0_EB, Global::REAL);
          actSDMat->GetScalar(ParameterMap["numS"], MAG_NUMS_EB, Global::REAL);
          actSDMat->GetScalar(ParameterMap["chi_factor"], MAG_CHI_FACTOR_EB, Global::REAL);
          ParameterMap["isMH"] = 0;
          matModelCoef_->InitModel(ParameterMap, actSDList);

          muNL = matModelCoef_; //actSDMat->GetTensorCoefFncModel(matModelCoef_);
          flux = matModelCoef_; //actSDMat->GetVectorCoefFncModel(matModelCoef_);
          nlFluxCoef_->AddRegion(actRegion, matModelCoef_);

          if (dim_ == 2) {
            stiffInt = new BDBInt<>(new GradientOperator<FeH1, 2>(), muNL, 1.0, updatedGeo_);
          }
          else {
            stiffInt = new BDBInt<>(new GradientOperator<FeH1, 3>(), muNL, 1.0, updatedGeo_);
          }
        }
        else {
          // classical nonlinearity with analytic prescription
          PtrCoefFct muNl = actSDMat->GetScalCoefFncNonLin( MAG_PERMEABILITY_SCALAR, Global::REAL, magFieldCoef);
          
          if (dim_ == 2) {
            stiffInt = new BBInt<>(new GradientOperator<FeH1, 2>(), muNl, 1.0, updatedGeo_);
          }
          else {
            stiffInt = new BBInt<>(new GradientOperator<FeH1, 3>(), muNl, 1.0, updatedGeo_);
          }
          perm_->AddRegion(actRegion, muNl);
        }

        stiffInt->SetName("StiffnessIntegratorHysteresis");
      }
      else {
        // ====================================================================
        // stiffness integrator - LINEAR
        // ====================================================================
        //get magnetic permeability
        PtrCoefFct coef, mu; 
        StdVector<NonLinType> matDepenTypes = regionMatDepTypes_[actRegion]; // material dependency
        if ( matDepenTypes.Find(PERMEABILITY) != -1 ) {
          // purely dependency of magnetic permeability of magnetic intensity (computed by forward problem)            
          StdVector<std::string> dispDofNames; // = feFunc_reduced->GetResultInfo()->dofNames;
          shared_ptr<EntityList> ent = ptGrid_->GetEntityList( EntityList::ELEM_LIST, regionName );

          //get coeff-Fnc to evaluate the temperature
          ReadMaterialDependency( "permeability", dispDofNames, ResultInfo::VECTOR, isComplex_,
                                  ent, coef, updatedGeo_ );
          //coef-Fnc for magnetic permeability
          mu = actSDMat->GetScalCoefFncNonLin( MAG_PERMEABILITY_SCALAR, Global::REAL, coef);     
          Hsmap_[actRegion] = coef;    
        }
        else {
          //get coeFunction for magnetic field intensity of the PDE
          Hsmap_[actRegion] = iterCplPde_->GetCouplingCoefFct(MAG_FIELD_INTENSITY, actSDList, "magneticScalarPotential", updatedGeo_);

          //get the reluctivity          
          mu = actSDMat->GetScalCoefFnc(MAG_PERMEABILITY_SCALAR, Global::REAL);
        }

        perm_->AddRegion(actRegion, mu);

        if (dim_ == 2) {
          stiffInt = new BBInt<>(new GradientOperator<FeH1, 2>(), perm_, 1.0, updatedGeo_);
        }
        else {
          stiffInt = new BBInt<>(new GradientOperator<FeH1, 3>(), perm_, 1.0, updatedGeo_);
        }
        stiffInt->SetName("StiffnessIntegrator");
      }

      BiLinFormContext *stiffIntDescr = new BiLinFormContext(stiffInt, STIFFNESS);

      stiffIntDescr->SetEntities(actSDList, actSDList);
      stiffIntDescr->SetFeFunctions(feFunc_reduced, feFunc_reduced);
      stiffInt->SetFeSpace(feSpace_reduced);

      assemble_->AddBiLinearForm(stiffIntDescr);
      bdbInts_.insert(std::pair<RegionIdType, BaseBDBInt *>(actRegion, stiffInt));
    }
  }

 void MagneticScalarPotentialAdjPDE::DefineRhsLoadIntegrators() {
    // get FEFunctions and space
    shared_ptr<BaseFeFunction> feFunc_reduced = feFunctions_[MAG_POTENTIAL_ADJ];
    shared_ptr<FeSpace> feSpace_reduced = feFunc_reduced->GetFeSpace();

    RegionIdType actRegion;
    bool coefUpdateGeo = true;

    // ========================================
    //  Current Vector Potential (T or Hs)
    // ========================================
    StdVector<shared_ptr<EntityList>> ent;
    StdVector<PtrCoefFct> coef;
    LinearForm* lin1 = NULL;
    LinearForm* lin2 = NULL;

    StdVector<std::string> dofNames(dim_);
    dofNames[0] = "Hx";
    dofNames[1] = "Hy";
    if ( dim_ == 3)
      dofNames[2] = "Hz";

    std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );

    //read the measured field intensities in the sensor positions and define RHS integrator
    ReadRhsExcitation("measuredFieldIntensity", dofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo);
    for (UInt reg = 0; reg < ent.GetSize(); ++reg) {
      // set current region and material
      actRegion = regions_[reg];

      UInt entDim = ent[reg]->GetGrid()->GetEntityDim(ptGrid_->GetRegionName(ent[reg]->GetRegion()));
      std::string regName = ptGrid_->GetRegionName(ent[reg]->GetRegion());
      std::string entName = ent[reg]->GetName();
      LOG_DBG(magscaladjpde) << "Entity with name " << entName << " is in region with name "
                          << entName << "\n and has the spatial dimension " << entDim;

      // check type of entitylist
      if (ent[reg]->GetType() == EntityList::NODE_LIST ||
          ent[reg]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EXCEPTION("Magnetic field intensity must be defined on elements")
      }

      //check for entity name
      if (regName.compare(entName) != 0) {
        EXCEPTION("There seems to be an error with region and entity names")
      }
     
      // //compute mu * Hmeas
      // actSDMat = materials_[actRegion];            
      // PtrCoefFct mu = actSDMat->GetScalCoefFnc(MAG_PERMEABILITY_SCALAR, Global::REAL);
      // PtrCoefFct muHmeas = CoefFunction::Generate( mp_,  Global::REAL,
      //                                              CoefXprBinOp(mp_, perm_, coef[reg], CoefXpr::OP_MULT ) );                                          
      if ( dim_ == 2 ) {
        lin1 = new BUIntegrator<Double>(new GradientOperator<FeH1, 2>(), 1.0, coef[reg], volRegions, coefUpdateGeo);        
      }
      else {
        lin1 = new BUIntegrator<Double>(new GradientOperator<FeH1, 3>(), 1.0, coef[reg], volRegions, coefUpdateGeo); 
      }      
      lin1->SetName("MeasuredMagFieldIntensityInt");
      lin1->SetNormalizeToVol();
      LinearFormContext *ctx = new LinearFormContext(lin1);
      ctx->SetEntities(ent[reg]);
      ctx->SetFeFunction(feFunc_reduced);
      assemble_->AddLinearForm(ctx);
    } // end loop over entities


    //read magnetic field intensity from forward computation 
    ReadRhsExcitation("fieldIntensity", dofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo);

    //here, we have to add the magnetic source field Hs from the current loaded coil
    //please note: we have do adapt the vector Hs and set all components to zero., which 
    //where not measured by the sensors; e.g., when the sensor just measures the x-component 
    //then we have to set the y- and z-component to zero, which we do by the scalVec!
    for (UInt i = 0; i < ent.GetSize(); ++i) {
      // set current region and material
      actRegion = regions_[i];
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];

      UInt entDim = ent[i]->GetGrid()->GetEntityDim(ptGrid_->GetRegionName(ent[i]->GetRegion()));
      std::string regName = ptGrid_->GetRegionName(ent[i]->GetRegion());
      std::string entName = ent[i]->GetName();
      LOG_DBG(magscaladjpde) << "Entity with name " << entName << " is in region with name "
                          << entName << "\n and has the spatial dimension " << entDim;

      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
          ent[i]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EXCEPTION("Magnetic field intensity must be defined on elements")
      }

      //check for entity name
      if (regName.compare(entName) != 0) {
        EXCEPTION("There seems to be an error with region and entity names")
      }

      // Here we store the H-field of the forward simulation
      //Hsmap_[ent[i]->GetRegion()] = coef[i];

      //PtrCoefFct coefB = Hsmap_[ent[i]->GetRegion()]; 
      Vector<Double> scalVec(dim_);
      scalVec.Init(0.0);
      std::istringstream iss(coef[i]->GetDofNames());
      do {
        std::string sub;
        iss >> sub;
        if( sub == "all" ) {
          // add all dofs to the definedDofs
          for( UInt k = 0; k<dim_; k++ ) {
            scalVec[k] = 1.0;
          }
          break;
        } 
        else {
            if ( sub == "Hx" ) 
              scalVec[0] = 1.0;
            else if ( sub == "Hy")
              scalVec[1] = 1.0;
            else if ( sub == "Hz" && dim_ == 3 )
              scalVec[2] = 1.0;
        }
      } while (iss);

      // actSDMat = materials_[actRegion];      
      // //we assembleHs permeability times magnetic field intensity of forward simulation
      // PtrCoefFct mu = actSDMat->GetScalCoefFnc(MAG_PERMEABILITY_SCALAR, Global::REAL);
      // PtrCoefFct muHs = CoefFunction::Generate( mp_, Global::REAL,
      //                                           CoefXprBinOp(mp_, mu, coef[i], CoefXpr::OP_MULT ) );

      if ( dim_ == 2 ) {
        lin2 = new BUIntegrator<Double>(new GradientOperator<FeH1, 2>(), -1.0, coef[i], volRegions, coefUpdateGeo);        
      }
      else {
        lin2 = new BUIntegrator<Double>(new GradientOperator<FeH1, 3>(), -1.0, coef[i], volRegions, coefUpdateGeo); 
      }
      lin2->SetScalVec(scalVec);
      lin2->SetNormalizeToVol();
      lin2->SetName("SourceMagFieldIntensityInt");
      LinearFormContext *ctx = new LinearFormContext(lin2);
      ctx->SetEntities(ent[i]);
      ctx->SetFeFunction(feFunc_reduced);
      assemble_->AddLinearForm(ctx);
    } // end loop over entities
 
  }


  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void MagneticScalarPotentialAdjPDE::InitNonLin() {
    SinglePDE::InitNonLin();
  }

  void MagneticScalarPotentialAdjPDE::DefineSolveStep()
  {
    if(nonLin_ && (modelName_ == "EBHysteresisModel")){
      solveStep_ = new SolveStepEB(*this);
    }else{
      solveStep_ = new StdSolveStep(*this);
    }
  }

  void MagneticScalarPotentialAdjPDE::InitTimeStepping()
  {
    Double gamma = 1.0;
    GLMScheme *scheme = new Trapezoidal(gamma);

    TimeSchemeGLM::NonLinType nlType = (nonLin_ || isHysteresis_) ? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType));
    feFunctions_[MAG_POTENTIAL_ADJ]->SetTimeScheme(myScheme);
  }

  void MagneticScalarPotentialAdjPDE::DefinePrimaryResults()
  {

    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL_ADJ];

    // Magnetic Potential
    shared_ptr<ResultInfo> res1(new ResultInfo);
    res1->resultType = MAG_POTENTIAL_ADJ;
    res1->dofNames = "";
    res1->unit = "A";
    res1->definedOn = ResultInfo::NODE;
    res1->entryType = ResultInfo::SCALAR;

    feFunctions_[MAG_POTENTIAL_ADJ]->SetResultInfo(res1);
    res1->SetFeFunction(feFunctions_[MAG_POTENTIAL_ADJ]);
    results_.Push_back(res1);
    availResults_.insert(res1);
    DefineFieldResult(feFunctions_[MAG_POTENTIAL_ADJ], res1);

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MAG_POTENTIAL_ADJ] = "ground";
    idbcSolNameMap_[MAG_POTENTIAL_ADJ] = "potential";
  }

  void MagneticScalarPotentialAdjPDE::DefinePostProcResults()
  {

    StdVector<std::string> vecDofNames;
    if (ptGrid_->GetDim() == 3) {
      vecDofNames = "x", "y", "z";
    } else {
      if (ptGrid_->IsAxi()) {
        vecDofNames = "r", "z";
      }
      else {
        vecDofNames = "x", "y";
      }
    }

    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL_ADJ];

    // === MAGNETIC RHS ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = MAG_RHS_LOAD_ADJ;
    rhs->dofNames = "";
    rhs->unit = "";
    rhs->entryType = ResultInfo::SCALAR;
    rhs->definedOn = ResultInfo::NODE;
    rhsFeFunctions_[MAG_POTENTIAL_ADJ]->SetResultInfo(rhs);
    DefineFieldResult( rhsFeFunctions_[MAG_POTENTIAL_ADJ], rhs );

    // === Minus GRADIENT PHI of ADJOINT ===
    shared_ptr<ResultInfo> hf(new ResultInfo);
    hf->resultType = MAG_POTENTIAL_GRAD_ADJ;
    hf->dofNames = vecDofNames;
    hf->unit = "A/m";
    hf->definedOn = ResultInfo::ELEMENT;
    hf->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionFormBased> hFunc;
    hFunc.reset(new CoefFunctionBOp<Double>(feFct, hf, 1.0));
    DefineFieldResult(hFunc, hf);
    stiffFormCoefs_.insert(hFunc); 

    // === Compute the volume
    shared_ptr<ResultInfo> vol(new ResultInfo);
    vol->resultType = VOLUME;
    vol->dofNames = "";
    vol->unit = "m^3";
    vol->entryType = ResultInfo::SCALAR;
    vol->definedOn = ResultInfo::REGION;
    shared_ptr<CoefFunction> coefVol = CoefFunction::Generate( mp_,  Global::REAL, "1.0");
    shared_ptr<ResultFunctor> volFunctor;
    volFunctor.reset(new ResultFunctorIntegrate<Double>(coefVol, feFct, vol));
    resultFunctors_[VOLUME] = volFunctor;
    availResults_.insert(vol);

    // === compute averaged magnetic field intensity of forward simulation
    shared_ptr<ResultInfo> averagedH(new ResultInfo);
    averagedH->resultType = MAG_FIELD_INTENSITY;
    averagedH->dofNames = vecDofNames;
    averagedH->unit = "A/m";
    averagedH->entryType = ResultInfo::VECTOR;
    averagedH->definedOn = ResultInfo::ELEMENT;
    // The computation is defined in FinalizePostProcResults()
    shared_ptr<CoefFunctionMulti> averagedHfnc(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
    DefineFieldResult(averagedHfnc, averagedH);    

    //functor for averaged h-field of forward simulation
    shared_ptr<ResultInfo> resH1(new ResultInfo());
    resH1->resultType = MAG_AVERAGED_FIELD_INTENSITY;
    resH1->dofNames = vecDofNames; //hField->GetDofNames();
    resH1->unit = "A/m";
    resH1->entryType = ResultInfo::VECTOR;
    resH1->definedOn = ResultInfo::REGION;   
    availResults_.insert(resH1);                             
    shared_ptr<ResultFunctor> averagedHfield;
    averagedHfield.reset(new ResultFunctorIntegrate<Double>(averagedHfnc, feFct, resH1));
    dynamic_pointer_cast< ResultFunctorIntegrate<Double> >(averagedHfield)->SetAveraged(true);
    resultFunctors_[MAG_AVERAGED_FIELD_INTENSITY] = averagedHfield;   
  }

  
  void MagneticScalarPotentialAdjPDE::FinalizePostProcResults() {
    StdVector<std::string> vecDofNames;
    if (ptGrid_->GetDim() == 3) {
      vecDofNames = "x", "y", "z";
    } else {
      if (ptGrid_->IsAxi()) {
        vecDofNames = "r", "z";
      }
      else {
        vecDofNames = "x", "y";
      }
    }

    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();

    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL_ADJ];

    // === Pure helper result ===
    shared_ptr<ResultInfo> ef(new ResultInfo);
    ef->resultType = OPT_RESULT_1; 
    ef->dofNames = "";
    ef->unit = "";
    ef->definedOn = ResultInfo::ELEMENT;
    ef->entryType = ResultInfo::SCALAR;   
    shared_ptr<CoefFunctionMulti> mult(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
    DefineFieldResult(mult, ef);

    // === Gradient via adjoint ===
    shared_ptr<ResultInfo> gradAdjParam;
    gradAdjParam.reset(new ResultInfo);
    gradAdjParam->resultType = MAG_GRAD_ADJ_PARAM;
    gradAdjParam->dofNames = "";
    gradAdjParam->unit = "";
    gradAdjParam->entryType = ResultInfo::SCALAR;
    gradAdjParam->definedOn = ResultInfo::REGION;    
    availResults_.insert(gradAdjParam);                            
    shared_ptr<ResultFunctor> gradRegion;
    gradRegion.reset(new ResultFunctorIntegrate<Double>(mult, feFct, gradAdjParam));
    resultFunctors_[MAG_GRAD_ADJ_PARAM] = gradRegion;

    // ============ MAG_GRAD_ADJ_PARAM & avaeraged flux density from forward simulation ============
    // define functor to compute gradient part: 
    shared_ptr<CoefFunctionMulti> scalMult = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[OPT_RESULT_1]);
    shared_ptr<CoefFunctionMulti> hField = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_FIELD_INTENSITY]);
    
    StdVector<RegionIdType>::iterator regIt = regions_.Begin();
    regIt = regions_.Begin();
    for (; regIt != regions_.End(); ++regIt) {
      // ========= H field =============
      PtrCoefFct Hforward = NULL;
      PtrCoefFct mult = NULL;
      RegionIdType actRegion = *regIt;
      if(Hsmap_.find(actRegion) != Hsmap_.end()){
        // H from previous (forward) simulation
        Hforward = Hsmap_[actRegion];
        //result grad(p), where p is the adjoint solution!
        PtrCoefFct gradPot = this->GetCoefFct(MAG_POTENTIAL_GRAD_ADJ);
        mult = CoefFunction::Generate( mp_, Global::REAL,
                                       CoefXprBinOp(mp_, Hforward, gradPot, CoefXpr::OP_MULT) );    
        scalMult->AddRegion(actRegion, mult);

        //H from previous (forward) simulation
        hField->AddRegion(actRegion, Hforward);
      }
    } // loop over regions

  }    
  // create the FEspace (H1 in this case)
  std::map<SolutionType, shared_ptr<FeSpace>>
  MagneticScalarPotentialAdjPDE::CreateFeSpaces(const std::string &formulation, PtrParamNode infoNode)
  {
    std::map<SolutionType, shared_ptr<FeSpace>> crSpaces;
    if (formulation == "default" || formulation == "H1")
    {
      PtrParamNode potSpaceNode_reduced = infoNode->Get("magPotentialAdj");
      crSpaces[MAG_POTENTIAL_ADJ] = FeSpace::CreateInstance(myParam_, potSpaceNode_reduced, FeSpace::H1, ptGrid_);
      crSpaces[MAG_POTENTIAL_ADJ]->Init(solStrat_);
    }
    else
    {
      EXCEPTION("The formulation " << formulation << "of adjoint magnetic potential PDE is not known!");
    }
    return crSpaces;
  }
}