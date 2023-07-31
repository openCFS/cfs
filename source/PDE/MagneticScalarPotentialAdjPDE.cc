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
    pdename_ = "magnetic";
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
          std::cout << "ReadMaterialDependency OK, now define coefFunction mu" << std::endl;
          //coef-Fnc for magnetic permeability
          mu = actSDMat->GetScalCoefFncNonLin( MAG_PERMEABILITY_SCALAR, Global::REAL, coef);
          std::cout << "coef-Fnc for magnetic permeability OK" << std::endl;           
        }
        else {
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

    BaseMaterial *actSDMat = NULL;
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
      LOG_DBG(magscaladjpde) << "Entity with name " << entName
                          << " is in region with name "
                          << entName
                          << "\n and has the spatial dimension "
                          << entDim;

      // check type of entitylist
      if (ent[reg]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Magnetic field intensity must be defined on elements")
      }
      if (regName.compare(entName) != 0) {
        EXCEPTION("There seems to be an error with region and entity names")
      }

      // Here we store the H-field of the previous (forward) simulation
      //std::cout << "measuredFieldIntensity: Add coef to Hsmap_ for region: " << regName << std::endl;
      Hsmap_[ent[reg]->GetRegion()] = coef[reg];
      
      actSDMat = materials_[actRegion];      
      PtrCoefFct magFieldCoef = this->GetCoefFct(MAG_FIELD_INTENSITY);   
      PtrCoefFct mu = actSDMat->GetScalCoefFnc(MAG_PERMEABILITY_SCALAR, Global::REAL);
      PtrCoefFct muHmeas = CoefFunction::Generate( mp_,  Global::REAL,
                                                   CoefXprBinOp(mp_, perm_, coef[reg], CoefXpr::OP_MULT ) );                                          
      if ( dim_ == 2 ) {
        lin1 = new BUIntegrator<Double>(new GradientOperator<FeH1, 2>(), -1.0, muHmeas, volRegions, coefUpdateGeo);        
      }
      else {
        lin1 = new BUIntegrator<Double>(new GradientOperator<FeH1, 3>(), -1.0, muHmeas, volRegions, coefUpdateGeo); 
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
      LOG_DBG(magscaladjpde) << "Entity with name " << entName
                          << " is in region with name "
                          << entName
                          << "\n and has the spatial dimension "
                          << entDim;

      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Magnetic field intensity must be defined on elements")
      }
      if (regName.compare(entName) != 0) {
        EXCEPTION("There seems to be an error with region and entity names")
      }

      // Here we store the H-field of the previous (forward) simulation
      //std::cout << "fieldIntensity: Add coef to Hsmap_ for region: " << regName << std::endl;
      Hsmap_[ent[i]->GetRegion()] = coef[i];

      //check dofs
      // PtrCoefFct coefOne  = CoefFunction::Generate(mp_, Global::REAL, "1.0");
      // PtrCoefFct coefZero = CoefFunction::Generate(mp_, Global::REAL, "0.0");
      // StdVector<PtrCoefFct> tensorComponents = StdVector<PtrCoefFct>(dim_*dim_);
      // tensorComponents.Init(coefZero);
      // tensorComponents[0] = coefOne;
      // tensorComponents[4] = coefOne;
      // tensorComponents[8] = coefOne;

      // StdVector<PtrCoefFct> vecComponents = StdVector<PtrCoefFct>(dim_);
      // vecComponents.Init(coefOne);

      // //Initialite wit zero
      // for (UInt k=0; k<dim_; k++)
      //   diagValues[k] = coefOne;

      // for (UInt k=dim_; k<2*dim_; k++)
      //   diagValues[k] = coefZero;

      // std::istringstream iss(coef[i]->GetDofNames());
      // do {
      //   std::string sub;
      //   iss >> sub;
      //   if( sub == "all" ) {
      //     // add all dofs to the definedDofs
      //     for( UInt k = 0; k<dim_; k++ ) {
      //       diagValues[k] = coefOne;
      //     }
      //     break;
      //   } 
      //   else {
      //       if ( sub == "Hx" ) 
      //         diagValues[0] = coefOne;
      //       else if ( sub == "Hy")
      //         diagValues[1] = coefOne;
      //       else if ( sub == "Hz" && dim_ == 3 )
      //         diagValues[2] = coefOne;
      //   }
      // } while (iss);

      //now we can set the scaling tensor (coefFunction)
      // PtrCoefFct measuredHcomp = CoefFunction::Generate(mp_,Global::REAL,dim_,dim_,tensorComponents);
      // PtrCoefFct vecHcomp = CoefFunction::Generate(mp_,Global::REAL,vecComponents);

      //shared_ptr<CoefFunction> measuredHcomp (new CoefFunctionDiagTensorFromScalar(diagValues));

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

      actSDMat = materials_[actRegion];      
      PtrCoefFct magFieldCoef = this->GetCoefFct(MAG_FIELD_INTENSITY);   
      //linear region and we assembleHs permeability times magnetic field intensity    
      PtrCoefFct mu = actSDMat->GetScalCoefFnc(MAG_PERMEABILITY_SCALAR, Global::REAL);
      PtrCoefFct muHs = CoefFunction::Generate( mp_,  Global::REAL,
                                                CoefXprBinOp(mp_, perm_, coef[i], CoefXpr::OP_MULT ) );
      // PtrCoefFct muHsMess = CoefFunction::Generate( mp_,  Global::REAL,
      //                                               CoefXprBinOp(mp_, measuredHcomp, vecHcomp, CoefXpr::OP_MULT ) );
      if ( dim_ == 2 ) {
        lin2 = new BUIntegrator<Double>(new GradientOperator<FeH1, 2>(), 1.0, muHs, volRegions, coefUpdateGeo);        
      }
      else {
        lin2 = new BUIntegrator<Double>(new GradientOperator<FeH1, 3>(), 1.0, muHs, volRegions, coefUpdateGeo); 
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
    res1->resultType = MAG_POTENTIAL;
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

    // === Pure helper result ===
    shared_ptr<ResultInfo> ef(new ResultInfo);
    ef->resultType = RHS_PSEUDO_DENSITY;
    ef->dofNames = "";
    ef->unit = "";
    ef->definedOn = ResultInfo::ELEMENT;
    ef->entryType = ResultInfo::SCALAR;
    // The computation of the gradient of the design parameters is defined in FinalizePostProcResults()
    shared_ptr<CoefFunctionMulti> mult(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
    DefineFieldResult(mult, ef);

    // === MAGNETIC RHS ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = MAG_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = "";
    rhs->entryType = ResultInfo::VECTOR;
    rhs->definedOn = ResultInfo::NODE;
    rhsFeFunctions_[MAG_POTENTIAL_ADJ]->SetResultInfo(rhs);
    DefineFieldResult( rhsFeFunctions_[MAG_POTENTIAL_ADJ], rhs );

    // === - GRADIENT PHI (also helper result)===
    shared_ptr<ResultInfo> hf(new ResultInfo);
    hf->resultType = MAG_POTENTIAL_GRAD;
    hf->dofNames = vecDofNames;
    hf->unit = "A/m";
    hf->definedOn = ResultInfo::ELEMENT;
    hf->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionFormBased> hFunc;
    hFunc.reset(new CoefFunctionBOp<Double>(feFct, hf, 1.0));
    DefineFieldResult(hFunc, hf);
    stiffFormCoefs_.insert(hFunc); 
    availResults_.insert(hf);  

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
    averagedH->definedOn = ResultInfo::REGION;
    // The computation is defined in FinalizePostProcResults()
    shared_ptr<CoefFunctionMulti> averagedHfnc(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
    DefineFieldResult(averagedHfnc, averagedH);       
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

    // ============ MAG_GRAD_ADJ_PARAM ============
    // define functor to compute gradient part: 

    shared_ptr<CoefFunctionMulti> scalMult = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[RHS_PSEUDO_DENSITY]);
    shared_ptr<CoefFunctionMulti> hField = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_FIELD_INTENSITY]);
    
    StdVector<RegionIdType>::iterator regIt = regions_.Begin();
    regIt = regions_.Begin();
    for (; regIt != regions_.End(); ++regIt) {
      // ========= H field =============
      PtrCoefFct hs, mult;
      if(Hsmap_.find(*regIt) != Hsmap_.end()){
        // H from previous (forward) simulation
        hs = Hsmap_[*regIt];
        // result MAG_POTENTIAL_DIV is grad(p), where p is the adjoint solution!
        PtrCoefFct gradPot = this->GetCoefFct(MAG_POTENTIAL_GRAD);
        mult = CoefFunction::Generate( mp_, Global::REAL,
                               CoefXprBinOp(mp_, gradPot, hs, CoefXpr::OP_MULT) );        
        scalMult->AddRegion(*regIt, mult);

        //H from previous (forward) simulation
        hField->AddRegion(*regIt, hs);
      }
    } // loop over regions

    // === Gradient via adjoint ===
    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL_ADJ];
    shared_ptr<ResultInfo> gradAdjParam;
    gradAdjParam.reset(new ResultInfo);
    gradAdjParam->resultType = MAG_GRAD_ADJ_PARAM;
    gradAdjParam->dofNames = "";
    gradAdjParam->unit = "";
    gradAdjParam->entryType = ResultInfo::SCALAR;
    gradAdjParam->definedOn = ResultInfo::REGION;    
    availResults_.insert(gradAdjParam);                            
    shared_ptr<ResultFunctor> gradRegion;
    gradRegion.reset(new ResultFunctorIntegrate<Double>(scalMult, feFct, gradAdjParam));
    resultFunctors_[MAG_GRAD_ADJ_PARAM] = gradRegion;

    //averaged hfield
    shared_ptr<ResultInfo> resH;
    resH.reset(new ResultInfo);
    resH->resultType = MAG_AVERAGED_FIELD_INTENSITY;
    resH->dofNames = vecDofNames; //hField->GetDofNames();
    resH->unit = "A/m";
    resH->entryType = ResultInfo::VECTOR;
    resH->definedOn = ResultInfo::REGION;   
    availResults_.insert(resH);                             
    shared_ptr<ResultFunctor> averagedH;
    averagedH.reset(new ResultFunctorIntegrate<Double>(hField, feFct, resH));
    dynamic_pointer_cast< ResultFunctorIntegrate<Double> >(averagedH)->SetAveraged(true);
    resultFunctors_[MAG_AVERAGED_FIELD_INTENSITY] = averagedH;
  }    

  // create the FEspace (H1 in this case)
  std::map<SolutionType, shared_ptr<FeSpace>>
  MagneticScalarPotentialAdjPDE::CreateFeSpaces(const std::string &formulation, PtrParamNode infoNode)
  {
    std::map<SolutionType, shared_ptr<FeSpace>> crSpaces;
    if (formulation == "default" || formulation == "H1")
    {
      PtrParamNode potSpaceNode_reduced = infoNode->Get("magPotential");
      crSpaces[MAG_POTENTIAL_ADJ] = FeSpace::CreateInstance(myParam_, potSpaceNode_reduced, FeSpace::H1, ptGrid_);
      crSpaces[MAG_POTENTIAL_ADJ]->Init(solStrat_);
    }
    else
    {
      EXCEPTION("The formulation " << formulation << "of magnetic PDE is not known!");
    }
    return crSpaces;
  }
}