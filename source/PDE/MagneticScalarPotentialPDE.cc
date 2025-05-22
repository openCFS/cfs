#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <set>

#include "MagneticScalarPotentialPDE.hh"

#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionMapping.hh"
#include "Utils/StdVector.hh"
#include "Driver/Assemble.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/FeFunctions.hh"

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

  DEFINE_LOG(magscalpde, "magneticScalarPotentialPDE")

  // ***************
  //   Constructor
  // ***************
  MagneticScalarPotentialPDE::MagneticScalarPotentialPDE(Grid *aptgrid, PtrParamNode paramNode,
                                                         PtrParamNode infoNode,
                                                         shared_ptr<SimState> simState, Domain *domain)
      : SinglePDE(aptgrid, paramNode, infoNode, simState, domain)
  {
    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_ = "magneticScalarPotential";
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

  void MagneticScalarPotentialPDE::DefineIntegrators()
  {

    RegionIdType actRegion;
    BaseMaterial *actSDMat = NULL;

    perm_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, false));

    // get FEFunctions and space
    shared_ptr<BaseFeFunction> feFunc_reduced = feFunctions_[MAG_POTENTIAL];
    shared_ptr<FeSpace> feSpace_reduced = feFunc_reduced->GetFeSpace();

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial *>::iterator it;

    PtrCoefFct magFieldCoef = this->GetCoefFct(MAG_FIELD_INTENSITY);

    // currently we are only allowed to have one hysteresis region
    bool moreThan1HystRegion = false;

    // iterate over the region (or materials)
    for (UInt iRegion = 0; iRegion < regions_.GetSize(); iRegion++)
    {
      // set current region and material
      actRegion = regions_[iRegion];
      actSDMat = materials_[actRegion];
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];
      
      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList(new ElemList(ptGrid_)); // changed newSDList to actSDList
      /* actSDList = newSDList; */
      actSDList->SetRegion(actRegion);

      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region", "name", regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();

      // --- Set the approximation for the current region ---
      feSpace_reduced->SetRegionApproximation(actRegion, polyId, integId);
      feFunc_reduced->AddEntityList(actSDList);

      // ====================================================================
      // STIFFNESS INTEGRATOR (start)
      // ====================================================================
      BaseBDBInt *stiffInt = NULL;
      PtrCoefFct muNL = NULL;
      PtrCoefFct flux = NULL;
      if (nonLinTypes.Find(PERMEABILITY) != -1)
      {
        // NONLINEAR CASE ( \int gradPhi dBdH gradPhi' )
        if (modelName_ == "JilesAthertonModel")
        {
          EXCEPTION("Jiles-Atherton model not implemented for MagneticScalarPotentialPDE");
        }
        else if (modelName_ == "EBHysteresisModel")
        {
          // if(moreThan1HystRegion){
          //   EXCEPTION("Currently only ONE hysteretic region is allowed!");
          // }
    
          
          // Init material model for hysteretic transient analysis
          if (((analysistype_ == STATIC) || (analysistype_ == TRANSIENT)) && nonLin_ && (modelName_ != "nonlinearCurve"))
          {
            matModelCoef_->Init(magFieldCoef, modelName_, dim_);
            matModelCoefm_[actRegion].reset(new CoefFunctionMaterialModel<Complex>()); // = matModelCoef_;
            matModelCoefm_[actRegion]->Init(magFieldCoef, modelName_, dim_); // = matModelCoef_;
          }

          //TODO kroppert: I do not like that, we should pass the actSDMat instead of these ParameterMap-stuff
          moreThan1HystRegion = true;
          std::map<std::string, double> ParameterMap;
          actSDMat->GetScalar(ParameterMap["Ps"], MAG_PS_EB, Global::REAL);
          actSDMat->GetScalar(ParameterMap["A"], MAG_A_EB, Global::REAL);
          actSDMat->GetScalar(ParameterMap["mu0"], MAG_MU0_EB, Global::REAL);
          actSDMat->GetScalar(ParameterMap["numS"], MAG_NUMS_EB, Global::REAL);
          actSDMat->GetScalar(ParameterMap["chi_factor"], MAG_CHI_FACTOR_EB, Global::REAL);
          actSDMat->GetScalar(ParameterMap["jacobian_method"], MAG_JACOBIAN_METHOD_EB, Global::REAL);
          ParameterMap["isMH"] = 0;
          
          //matModelCoef_->InitModel(ParameterMap, actSDList);          
          matModelCoefm_[actRegion]->InitModel(ParameterMap, actSDList);

          if(actSDMat->GetAnhystMagModel() == "multiscale_anhysteresis"){
            // check if we have a dependency of the anhysteretic curve with mechanical stress
            PtrCoefFct stresscoef;
            //get coeff-Fnc to evaluate the temperature
            StdVector<std::string> dispDofNames = feFunc_reduced->GetResultInfo()->dofNames;
            shared_ptr<EntityList> ent = ptGrid_->GetEntityList( EntityList::ELEM_LIST, regionName );
            ReadMaterialDependency( "permeability", dispDofNames, ResultInfo::VECTOR, false,
                                    ent, stresscoef, updatedGeo_ );
            matModelCoef_->RegisterStressDependence(stresscoef);
          }

          //muNL = matModelCoef_;          
          //nlFluxCoef_->AddRegion(actRegion, matModelCoef_);
          muNL = matModelCoefm_[actRegion];
          nlFluxCoefm_[actRegion].reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_, true)); //= nlFluxCoef_;
          nlFluxCoefm_[actRegion]->AddRegion(actRegion, matModelCoefm_[actRegion]);
          //nlFluxCoef_->AddRegion(actRegion, matModelCoefm_[actRegion]);

          if (dim_ == 2)
          {
            stiffInt = new BDBInt<>(new GradientOperator<FeH1, 2>(), muNL, 1.0, updatedGeo_);
          }
          else
          {
            stiffInt = new BDBInt<>(new GradientOperator<FeH1, 3>(), muNL, 1.0, updatedGeo_);
          }
        }
        stiffInt->SetName("(db/dh grad phi, grad phi')");
      }
      else
      {
        // LINEAR CASE ( \int gradPhi \mu gradPhi' )
        PtrCoefFct mu = actSDMat->GetScalCoefFnc(MAG_PERMEABILITY_SCALAR, Global::REAL);
        perm_->AddRegion(actRegion, mu);
        if (dim_ == 2)
        {
          stiffInt = new BBInt<>(new GradientOperator<FeH1, 2>(), perm_, 1.0, updatedGeo_);
        }
        else
        {
          stiffInt = new BBInt<>(new GradientOperator<FeH1, 3>(), perm_, 1.0, updatedGeo_);
        }
        stiffInt->SetName("(mu grad phi, grad phi')");
      }
      BiLinFormContext *stiffIntDescr = new BiLinFormContext(stiffInt, STIFFNESS);
      stiffIntDescr->SetEntities(actSDList, actSDList);
      stiffIntDescr->SetFeFunctions(feFunc_reduced, feFunc_reduced);
      stiffInt->SetFeSpace(feSpace_reduced);
      assemble_->AddBiLinearForm(stiffIntDescr);
      bdbInts_.insert(std::pair<RegionIdType, BaseBDBInt *>(actRegion, stiffInt));
      // ====================================================================
      // STIFFNESS INTEGRATOR (end)
      // ====================================================================
    }
  }

  void MagneticScalarPotentialPDE::DefineRhsLoadIntegrators()
  {

    // get FEFunctions and space
    shared_ptr<BaseFeFunction> feFunc_reduced = feFunctions_[MAG_POTENTIAL];
    shared_ptr<FeSpace> feSpace_reduced = feFunc_reduced->GetFeSpace();

    LinearForm * lin = NULL;
    RegionIdType actRegion;

    bool coefUpdateGeo = true;
    bool isHystereticMat = false;


    // iterate over the region (or materials)
    for (UInt iRegion = 0; iRegion < regions_.GetSize(); iRegion++)
    {
      // set current region and material
      actRegion = regions_[iRegion];
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];

      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList(new ElemList(ptGrid_));
      actSDList->SetRegion(actRegion);

      PtrCoefFct fluxDensityNL = NULL;
      fluxDensityNL = matModelCoef_;

      // ===============================================================================================
      // NONLINEAR CASE AND NONLINEAR REGION: \int B(H) \gradPhi' (start)
      // ===============================================================================================
      if (nonLinTypes.Find(PERMEABILITY) != -1) {
        if (modelName_ == "JilesAthertonModel")
        {
          EXCEPTION("Jiles-Atherton model not implemented for MagneticScalarPotentialPDE");
        }
        else if (modelName_ == "EBHysteresisModel")
        {
          isHystereticMat = true;
          if( dim_ == 2 ) {
            lin = new BUIntegrator<Double>( new GradientOperator<FeH1,2> (),
                    (1.0), fluxDensityNL, coefUpdateGeo, false);
          } else {
            lin = new BUIntegrator<Double>( new GradientOperator<FeH1,3> (),
                    (1.0), fluxDensityNL, coefUpdateGeo, false);
          }
        }
      }else{
        EXCEPTION("MagneticScalarPotentialPDE: no valid material model provided");
      }

    lin->SetName("(B,grad phi'): nonlinear problem; nonlinear subregion RHS");
    lin->SetSolDependent();
    LinearFormContext *ctx = new LinearFormContext( lin );
    ctx->SetEntities( actSDList );
    ctx->SetFeFunction(feFunc_reduced);
    assemble_->AddLinearForm(ctx);
    }
    // ===============================================================================================
    // NONLINEAR CASE AND NONLINEAR REGION: \int B(H) \gradPhi' (end)
    // ===============================================================================================


    StdVector<shared_ptr<EntityList>> ent;
    StdVector<PtrCoefFct> coef;
    StdVector<PtrCoefFct> remanent_flux_density;
    StdVector<std::string> vecDofNames = feFunc_reduced->GetResultInfo()->dofNames;
    LinearForm *lin2 = NULL;
    LinearForm *lin_br = NULL; // (b_r,gradN), describes linear permanent magnets in the linear/nonlinear case
    StdVector<std::string> dofNames;
    std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );
    ReadRhsExcitation("fieldIntensity", dofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo);
    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
      UInt entDim = ent[i]->GetGrid()->GetEntityDim(ptGrid_->GetRegionName(ent[i]->GetRegion()));
      std::string regName = ptGrid_->GetRegionName(ent[i]->GetRegion());
      std::string entName = ent[i]->GetName();
      LOG_DBG(magscalpde) << "Entity with name " << entName
                          << " is in region with name "
                          << entName
                          << "\n and has the spatial dimension "
                          << entDim;

      actRegion = ptGrid_->GetRegionId(regName);
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];
    

      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST)
      {
        EXCEPTION("Magnetic field intensity must be defined on elements")
      }
      if (regName.compare(entName) != 0)
      {
        EXCEPTION("There seems to be an error with region and entity names")
      }
      
      // ===============================================================================================
      // RHS FOR NONLINEAR CASE BUT LINEAR SUBREGION OR ENTIRELY LINEAR (start)
      // ===============================================================================================
      Hsmap_[ent[i]->GetRegion()] = coef[i]; // Here we store the Hs field for every region to have it ready for postprocessing
      if ( nonLinTypes.Find(PERMEABILITY) == -1){ 
        if (isHystereticMat) // NONLINEAR CASE BUT LINEAR REGION: \int B gradPhi'
        {
          if(dim_ == 2){
            lin2 = new BUIntegrator<Double>(new GradientOperator<FeH1, 2>(), 1.0, GetCoefFct( MAG_FLUX_DENSITY ), volRegions, coefUpdateGeo);
          }else{
            lin2 = new BUIntegrator<Double>(new GradientOperator<FeH1, 3>(), 1.0, GetCoefFct( MAG_FLUX_DENSITY ), volRegions, coefUpdateGeo);
          }
          lin2->SetName("(B,grad phi'): nonlinear problem; linear subregion RHS");
          LinearFormContext *ctx = new LinearFormContext(lin2);
          ctx->SetEntities(ent[i]);
          ctx->SetFeFunction(feFunc_reduced);
          assemble_->AddLinearForm(ctx);
        }
        else // LINEAR CASE : \int \mu Hs gradPhi'
        {
          PtrCoefFct mu_times_Hs = CoefFunction::Generate(mp_, Global::REAL, CoefXprVecScalOp(mp_, coef[i], perm_, CoefXpr::OP_MULT));  
          if(dim_ == 2){
            lin2 = new BUIntegrator<Double>(new GradientOperator<FeH1, 2>(), 1.0, mu_times_Hs, volRegions, coefUpdateGeo);
          }else{
            lin2 = new BUIntegrator<Double>(new GradientOperator<FeH1, 3>(), 1.0, mu_times_Hs, volRegions, coefUpdateGeo);
          }

          lin2->SetName("(mu Hs,grad phi'): linear RHS");
          LinearFormContext *ctx = new LinearFormContext(lin2);
          ctx->SetEntities(ent[i]);
          ctx->SetFeFunction(feFunc_reduced);
          assemble_->AddLinearForm(ctx); 
        }
      }
    // ===============================================================================================
    // RHS FOR NONLINEAR CASE BUT LINEAR SUBREGION OR ENTIRELY LINEAR (end)
    // ===============================================================================================
    } // end loop over entities

    // ===============================================================================================
    // ADDING REMANENCE FLUX DENSITY TO RHS: lin_br: (b_r,gradN) (start)
    // The b_r has to be defined in the input .xml file for the simulation of permanent magnets
    // ===============================================================================================
    ReadRhsExcitation("fluxDensity", vecDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );

    PtrCoefFct br;
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
          ent[i]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EXCEPTION("Prescribed magnetic flux density can only be specified in a volume!")
      }

      //save the remanent magnetic flux of the permanent magnet
      BremMap_[ent[i]->GetRegion()] = coef[i]; 

      if (isComplex_) {
        if( dim_ == 2) {
          lin = new BUIntegrator<Complex>( new GradientOperator<FeH1,2,1,Complex>(), Complex(1.0),coef[i],  coefUpdateGeo, true);
        } else {
          // 3D case
          lin = new BUIntegrator<Complex> ( new GradientOperator<FeH1,3,1,Complex>(), Complex(1.0),coef[i],  coefUpdateGeo, true);
        }
      } else {
        if( dim_ == 2) {
          lin = new BUIntegrator<Double>( new GradientOperator<FeH1,2>(), 1.0, coef[i],  coefUpdateGeo, true);
        } else {
          // 3D case
          lin = new BUIntegrator<Double>( new GradientOperator<FeH1,3>(), 1.0, coef[i],  coefUpdateGeo, true);         
        }
      }
      
      lin_br->SetName("fluxIntegrator: (b_r,gradN): linear RHS, remanent flux density");
      LinearFormContext *ctx = new LinearFormContext( lin_br );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(feFunc_reduced);
      assemble_->AddLinearForm(ctx);
      feFunc_reduced->AddEntityList(ent[i]);
    }
    // ===============================================================================================
    //  ADDING REMANENCE FLUX DENSITY TO RHS: lin_br: (b_r,gradN) (end)
    // ===============================================================================================
  }

  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void MagneticScalarPotentialPDE::InitNonLin() {
    SinglePDE::InitNonLin();
  }

  void MagneticScalarPotentialPDE::DefineSolveStep()
  {
    if(nonLin_ && (modelName_ == "EBHysteresisModel")){
      solveStep_ = new SolveStepEB(*this);
    }else{
      solveStep_ = new StdSolveStep(*this);
    }
  }

  void MagneticScalarPotentialPDE::InitTimeStepping()
  {
    Double gamma = 1.0;
    GLMScheme *scheme = new Trapezoidal(gamma);

    TimeSchemeGLM::NonLinType nlType = (nonLin_ || isHysteresis_) ? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType));
    feFunctions_[MAG_POTENTIAL]->SetTimeScheme(myScheme);
  }

  void MagneticScalarPotentialPDE::DefinePrimaryResults()
  {

    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];

    // Magnetic Potential
    shared_ptr<ResultInfo> res1(new ResultInfo);
    res1->resultType = MAG_POTENTIAL;
    res1->dofNames = "";
    res1->unit = "A";
    res1->definedOn = ResultInfo::NODE;
    res1->entryType = ResultInfo::SCALAR;
    feFunctions_[MAG_POTENTIAL]->SetResultInfo(res1);
    res1->SetFeFunction(feFunctions_[MAG_POTENTIAL]);
    results_.Push_back(res1);
    availResults_.insert(res1);
    DefineFieldResult(feFunctions_[MAG_POTENTIAL], res1);

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MAG_POTENTIAL] = "ground";
    idbcSolNameMap_[MAG_POTENTIAL] = "potential";
  }

  void MagneticScalarPotentialPDE::DefinePostProcResults()
  {

    StdVector<std::string> vecDofNames;
    if (ptGrid_->GetDim() == 3)
    {
      vecDofNames = "x", "y", "z";
    }
    else
    {
      if (ptGrid_->IsAxi())
      {
        vecDofNames = "r", "z";
      }
      else
      {
        vecDofNames = "x", "y";
      }
    }

    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];

    // === MAGNETIC FIELD INTENSITY ===
    shared_ptr<ResultInfo> ef(new ResultInfo);
    ef->resultType = MAG_FIELD_INTENSITY;
    ef->dofNames = vecDofNames;
    ef->unit = "A/m";
    ef->definedOn = ResultInfo::ELEMENT;
    ef->entryType = ResultInfo::VECTOR;
    ef->SetFeFunction(feFunctions_[MAG_POTENTIAL]);
    // The recipe about how to actually evaluate H, is defined in FinalizePostProcResults()
    shared_ptr<CoefFunctionMulti> hIntensFunc(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
    DefineFieldResult(hIntensFunc, ef);

    // === PERMEABILITY  ===
    shared_ptr<ResultInfo> perm(new ResultInfo);
    perm->resultType = MAG_ELEM_PERMEABILITY;
    perm->SetFeFunction(feFunctions_[MAG_POTENTIAL]);
    shared_ptr<CoefFunctionFormBased> perm_coef;
    shared_ptr<CoefFunctionMulti> permFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, false));
    if(nonLin_ && (modelName_ == "EBHysteresisModel")){
      if(dim_ == 2){
        perm->dofNames = "xx", "yy","xy";
      } else {
        perm->dofNames = "xx", "yy","zz","yz","xz","xy";
      }
      perm->unit = "Vs/Am";
      perm->definedOn = ResultInfo::ELEMENT;
      perm->entryType = ResultInfo::TENSOR;
      perm_coef.reset(new CoefFunctionHomogenization<double, App::MAG>(feFct));
      DefineFieldResult(perm_coef , perm);
      stiffFormCoefs_.insert(perm_coef);
    }else{
      perm->dofNames = "";
      perm->unit = "Vs/Am";
      perm->definedOn = ResultInfo::ELEMENT;
      perm->entryType = ResultInfo::SCALAR;
      DefineFieldResult( perm_, perm );
    }

    // =====================================================
    // MAG_FLUX_DENSTIY (START)
    // =====================================================
    shared_ptr<ResultInfo> bf(new ResultInfo);
    bf->resultType = MAG_FLUX_DENSITY;
    bf->SetFeFunction(feFunctions_[MAG_POTENTIAL]);
    bf->dofNames = vecDofNames;
    bf->unit = "T";
    bf->definedOn = ResultInfo::ELEMENT;
    bf->entryType = ResultInfo::VECTOR;
    // The recipe about how to actually evaluate B, is defined in FinalizePostProcResults()
    shared_ptr<CoefFunctionMulti> bFunc(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
    DefineFieldResult(bFunc, bf);
    // =====================================================
    // MAG_FLUX_DENSTIY (END)
    // =====================================================

    // =====================================================
    // GRADIENT PHI (helper result) (START)
    // =====================================================
    shared_ptr<ResultInfo> hf(new ResultInfo);
    hf->resultType = MAG_POTENTIAL_DIV;
    hf->SetFeFunction(feFunctions_[MAG_POTENTIAL]);
    hf->dofNames = vecDofNames;
    hf->unit = "A/m";
    hf->definedOn = ResultInfo::ELEMENT;
    hf->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionFormBased> hFunc;
    hFunc.reset(new CoefFunctionBOp<Double>(feFct, hf, 1.0));
    DefineFieldResult(hFunc, hf);
    stiffFormCoefs_.insert(hFunc);
    // =====================================================
    // GRADIENT PHI (helper result) (END)
    // =====================================================
  }


  // **********************************************************
  // recipe to do the postprocessing for certain quantities:
  // MAG_FLUX_DENSTIY
  // MAG_FIELD_INTENSITY
  // **********************************************************
  void MagneticScalarPotentialPDE::FinalizePostProcResults()
  {
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;

    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();

    // ============ MAGNETIC FIELD INTENSITY ============
    // Assemble coefficient function for
    // H = Hs - \grad \Phi
    // whereas Hs is stored in a CoefFunction and typically read in from a previous sequence
    // step as the solution of a curl curl type problem (magnetostatics)

    // Hs
    shared_ptr<CoefFunctionMulti> hIntensCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_FIELD_INTENSITY]);
    shared_ptr<CoefFunctionMulti> bCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_FLUX_DENSITY]);
    //shared_ptr<CoefFunctionMulti> mu = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_ELEM_PERMEABILITY]);
    // We only need this to check which kind of nonlinearity is specified in the different regions
    StdVector<RegionIdType>::iterator regIt = regions_.Begin();
    regIt = regions_.Begin();

    for (; regIt != regions_.End(); ++regIt){
      // =====================================================
      // MAG_FIELD_INTENSITY (START)
      // h = hs - gradPhi
      // =====================================================
      PtrCoefFct hs, magIntens;
      if(Hsmap_.find(*regIt) != Hsmap_.end()){
        // Hs is prescribed in the region
        hs = Hsmap_[*regIt];
        // result MAG_POTENTIAL_DIV is (grad(Phi))!!
        magIntens = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, hs, GetCoefFct(MAG_POTENTIAL_DIV), CoefXpr::OP_SUB));        
        hIntensCoef->AddRegion(*regIt, magIntens);
      }else{
        // Hs is NOT prescribed in the region
        hIntensCoef->AddRegion(*regIt, GetCoefFct(MAG_POTENTIAL_DIV));
      }
      // =====================================================
      // MAG_FIELD_INTENSITY (END)
      // =====================================================


      // =====================================================
      // MAG_FLUX_DENSTIY (START)
      // b = mu*h + b_r (linear case)
      // b = b(h) + b_r (nonlinear/hysteresis case)
      // =====================================================
      // Just to find out which linear/nonlinear type is defined in this region
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[*regIt];
      if( nonLinTypes.Find(PERMEABILITY) != -1 ){
        // hysteretic case
        bCoef->AddRegion(*regIt, nlFluxCoefm_[*regIt]);
      }else{
        // classical nonlinear case and linear case
        PtrCoefFct b = CoefFunction::Generate( mp_, part, CoefXprVecScalOp(mp_, hIntensCoef, perm_, CoefXpr::OP_MULT));

        //check for remanent flux density of a permanent magnet
        PtrCoefFct btotal;
        if(BremMap_.find(*regIt) != BremMap_.end()){          
          PtrCoefFct br = BremMap_[*regIt];
          btotal = b = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, b, br, CoefXpr::OP_ADD));
        } else {
          btotal = b;          
        }
        bCoef->AddRegion(*regIt, b);
      }
      // =====================================================
      // MAG_FLUX_DENSTIY (END)
      // =====================================================    
    } // loop over regions
  }

  // create the FEspace (H1 in this case)
  std::map<SolutionType, shared_ptr<FeSpace>>
  MagneticScalarPotentialPDE::CreateFeSpaces(const std::string &formulation, PtrParamNode infoNode)
  {
    std::map<SolutionType, shared_ptr<FeSpace>> crSpaces;
    if (formulation == "default" || formulation == "H1")
    {
      PtrParamNode potSpaceNode_reduced = infoNode->Get("magPotential");
      crSpaces[MAG_POTENTIAL] = FeSpace::CreateInstance(myParam_, potSpaceNode_reduced, FeSpace::H1, ptGrid_);
      crSpaces[MAG_POTENTIAL]->Init(solStrat_);
    }
    else
    {
      EXCEPTION("The formulation " << formulation << "of magnetic PDE is not known!");
    }
    return crSpaces;
  }
}
