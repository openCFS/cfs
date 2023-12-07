#include <fstream>

#include "MagEdgeHPDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/Coil.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"

#include "Driver/Assemble.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "FeBasis/HCurl/FeSpaceHCurlHi.hh"
#include "FeBasis/HCurl/HCurlElems.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "Domain/CoefFunction/CoefFunctionHarmBalance.hh"
#include "Domain/CoefFunction/CoefFunctionExpression.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"



// forms
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/BiLinForms/BiLinWrappedLinForm.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/Operators/CurlOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/ConvectiveOperator.hh"
#include "Forms/Operators/GradientOperator.hh"

// solve steps
#include "Driver/SolveSteps/SolveStepEB.hh"

//time stepping
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

#include "Driver/MultiHarmonicDriver.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
namespace CoupledField {

// declare class specific logging stream
DEFINE_LOG(magEdgeHPde, "magEdgeHPde")


  // **************
  //  Constructor
  // **************
  MagEdgeHPDE::MagEdgeHPDE( Grid * aptgrid, PtrParamNode paramNode,
                          PtrParamNode infoNode,
                          shared_ptr<SimState> simState, Domain* domain )
    :MagBasePDE( aptgrid, paramNode, infoNode, simState, domain ) {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magneticHEdge";
    pdematerialclass_ = ELECTROMAGNETIC;
    formulation_ = MagBasePDE::EDGE;

    nonLin_ = false;
    nonLinMaterial_ = false;

    //! Always use updated Lagrangian formulation
    updatedGeo_        = true;

    // default is false
    useGradFields_ = paramNode->Get("useGradientFields")->As<bool>();

    // this is the default value
    modelName_ = "nonlinearCurve";

    // Be aware that the hysteresis via CoefFunctionMaterialModel does not use the isHysteresis_ flag!
    matModelCoef_.reset(new CoefFunctionMaterialModel<Complex>());
    // init the nlFluxCoef
    nlFluxCoef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_, true));
    nlScalCoef_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_, true));

    if ( !is3d_ )
      EXCEPTION("MagEdgeHPDE is just implemented for 3D setups!");
  }


  // **********************************************************
  //  Destructor
  // **********************************************************
  MagEdgeHPDE::~MagEdgeHPDE(){
  }


  // **********************************************************
  //  Initialize Nonlinearities
  // **********************************************************
  void MagEdgeHPDE::InitNonLin(){
    SinglePDE::InitNonLin();
  }


  // **********************************************************
  // Define SolveStep depending on the material law
  // nonlinear/hysteresis: SolveStepEB
  // linear:               StdSolveStep
  // **********************************************************
  void MagEdgeHPDE::DefineSolveStep(){
    if(nonLin_ && (modelName_ == "EBHysteresisModel")){
      solveStep_ = new SolveStepEB(*this);
    } 
    else{
      solveStep_ = new StdSolveStep(*this);
    }
  }


  // **********************************************************
  // INIT. TIME STEPPING (not a real time stepping, since the
  // resulting FEM system is not an ODE)
  // **********************************************************
  void MagEdgeHPDE::InitTimeStepping(){
    Double gamma = 1.0;
    GLMScheme *scheme = new Trapezoidal(gamma);

    TimeSchemeGLM::NonLinType nlType = (nonLin_ || isHysteresis_) ? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType));
    feFunctions_[MAG_FIELD_INTENSITY]->SetTimeScheme(myScheme);
  }


  // **********************************************************
  //  Definition of Integrators
  // **********************************************************
  void MagEdgeHPDE::DefineIntegrators(){
    this->DefineStandardIntegrators();
    DefineCoilIntegrators(); // in MagBasePDE
  }


  // **********************************************************
  // just a helper function for DefineCoilIntegrators()
  // **********************************************************
  LinearForm* MagEdgeHPDE::GetCurrentDensityInt( Double factor, PtrCoefFct coef, std::string coilId){
    LinearForm * ret = NULL;
    ret = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3, Double>(), factor, coef, updatedGeo_);
    return ret;
  }


  // *************************************************************
  // Definition of linearforms that involve coils
  // linear:    (rho_art j,curlN)
  // nonlinear: - (rho_art j,curlN); negative due to Newton scheme
  // *************************************************************
  void MagEdgeHPDE::DefineCoilIntegrators(){

    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;

    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_FIELD_INTENSITY];
    shared_ptr<FeSpace> feSpace = feFunc->GetFeSpace();

    std::map<Coil::IdType, shared_ptr<Coil>>::iterator coilIt;
    coilIt = coils_.begin();
    for (; coilIt != coils_.end(); coilIt++){
      Coil &actCoil = *(coilIt->second);
      // run over all parts
      std::map<RegionIdType, shared_ptr<Coil::Part>>::iterator partIt;
      partIt = actCoil.parts_.begin();

      if (actCoil.sourceType_ != Coil::CURRENT){
        EXCEPTION("Only current excitation is implemented in H-formulation!!");
      } 
      else{

        for (partIt = actCoil.parts_.begin(); partIt != actCoil.parts_.end(); partIt++){
          Coil::Part &actPart = *(partIt->second);
          RegionIdType actRegion = partIt->first;
          shared_ptr<ElemList> actSDList(new ElemList(ptGrid_));
          actSDList->SetRegion(actRegion);
          LinearForm* curInt = NULL;

          std::map<UInt, PtrCoefFct> jFct;
          CoefXprVecScalOp iVec = CoefXprVecScalOp(mp_, actPart.jUnitVec, actCoil.srcVal_, CoefXpr::OP_MULT);
          PtrCoefFct iFct = CoefFunction::Generate(mp_, part, iVec);

          // ============= IMPORTANT (START) =======================================
          // Here the value for the current density vector j is determined.
          // The linearform of the formulation is as follows
          //         (rho_art j, curl(h')) ,
          // where rho_art is an artificial resistivity that occurs due to the
          // penalty method to end up in this formulation.
          // A common way to calculate this coefficient is
          //         rho_art = 1/(n^0.01)*mu_rmax*mu*10^(n/2) ,
          // where mu is the permeability in this region and n is a penalty parameter
          // that should be chosen rather small. However, since coils are usually 
          // in regions where mu ~ mu0 this factor is harcoded and only n can be 
          // altered.  Anyway, this should work also for the nonlinear hysteretic case,
          // since the coil region is usually linear material.
          // (Honestly, I did not know how to implement it differently)
          // ============= IMPORTANT (END) ======================================= 
          Double penaltyParameter = 1e-6; // default value for n
          myParam_->GetValue("penaltyFunctionParameter", penaltyParameter, ParamNode::PASS); // get the penaltyParameter from the .xml file and set the default value
          std::cout << "penalty j curlh" << penaltyParameter << std::endl;
          //CoefXprVecScalOp jVec = CoefXprVecScalOp(mp_, iFct, boost::lexical_cast<std::string>(actPart.wireCrossSect*(1/((1/(std::pow(penaltyParameter,0.01)))*1*1.256637061e-06*std::pow(10,penaltyParameter/2)))), CoefXpr::OP_DIV);
          CoefXprVecScalOp jVec = CoefXprVecScalOp(mp_, iFct, boost::lexical_cast<std::string>(actPart.wireCrossSect), CoefXpr::OP_DIV);
          PtrCoefFct jVec_coefFct = CoefFunction::Generate(mp_, part, jVec);
          CoefXprVecScalOp rho_art_times_jVec = CoefXprVecScalOp(mp_, jVec_coefFct, boost::lexical_cast<std::string>(1.256637061e-06*std::pow(10,penaltyParameter/2)), CoefXpr::OP_MULT);
          jFct[0] = CoefFunction::Generate(mp_, part, rho_art_times_jVec);
          coilCurrentDens_[actRegion] = jFct[0];
          curInt = GetCurrentDensityInt( 1.0, jFct[0] );
          curInt->SetName("(rho_art j,curlN): CoilIntegrator");
          LinearFormContext * coilContext = new LinearFormContext( curInt );
          coilContext->SetEntities( actSDList );
          coilContext->SetFeFunction( feFunc );
          assemble_->AddLinearForm( coilContext );
        } // loop: parts
      }
    } // loop: coils
  }


  // **********************************************************
  // Definition of all bilinearforms
  // linear:    (rho_art curlN, curlN) + (mu N, N)
  // nonlinear: (rho_art curlN, curlN) + (db/dh N, N)
  // **********************************************************
  void MagEdgeHPDE::DefineStandardIntegrators(){

    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

    // global storage for permeability makes the postprocessing easier
    perm_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, false));

    // get the penaltyParameter from the .xml file and set the default value
    Double penaltyParameter = 1e-6;
    myParam_->GetValue("penaltyFunctionParameter", penaltyParameter, ParamNode::PASS);
    std::cout << "penalty bilinear" << penaltyParameter << std::endl;

    std::cout << penaltyParameter << std::endl;

    // get FEFunction and space
    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_FIELD_INTENSITY];
    shared_ptr<FeSpace> feSpace = feFunc->GetFeSpace();

    PtrCoefFct magFieldCoef = this->GetCoefFct(MAG_FIELD_INTENSITY);
    // Init material model for hysteretic transient analysis
    if (((analysistype_ == STATIC) || (analysistype_ == TRANSIENT)) && nonLin_ && (modelName_ != "nonlinearCurve")){
      matModelCoef_->Init(magFieldCoef, modelName_, dim_);
    }

    // currently we are only allowed to have one hysteresis region
    bool moreThan1HystRegion = false;
    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion ++){
        // set current region and materials
        actRegion = regions_[iRegion];
        actMat    = materials_[actRegion];
        StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];

        // get current region name
        std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

        // create new entity list
        shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
        actSDList->SetRegion( actRegion );

        // Set the FE-Ansatz for the current region
        PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
        std::string polyId = curRegNode->Get("polyId")->As<std::string>();
        std::string integId = curRegNode->Get("integId")->As<std::string>();

        // pass entitylist ot fespace / fefunction
        feSpace->SetRegionApproximation(actRegion,polyId,integId);
        feFunc->AddEntityList( actSDList );

        // ====================================================================
        // MASS INTEGRATOR [START] (c N, N)
        // The coefficient c of this integrator depends on the material case:
        // ~ linear:     use just the value of the permeability of this region.
        // ~ nonlinear:  use the derivative of the material law b = b(h) as the
        //             value of the coefficient, since this is the Jacobian J 
        //             of the Newton-Raphson scheme J = db/dh.
        // ~ hysteresis: same as in the nonlinear case. However typically the 
        //             Jacobian can not be determined analytically, but by 
        //             a numerical method like finite differences or Broyden.
        // So here the only thing that is "if cased", is the determination of
        // coefficient and the integrator that accepts also tensors as 
        // material coefficent. The integrator (structure, (c N, N)) itself
        // stays the same!
        // ====================================================================
        PtrCoefFct mu = NULL;
        BaseBDBInt *massInt = NULL;
        if (nonLinTypes.Find(PERMEABILITY) != -1){ // NONLINEAR/HYSTERESIS CASE
          if (modelName_ == "JilesAthertonModel"){ // just because, nobody wonders that it does not work with this model
            EXCEPTION("Jiles-Atherton model not implemented for MagEdgeHPDE");
          }
          else if (modelName_ == "EBHysteresisModel"){ // (this is the model we want to use!)
            if(moreThan1HystRegion){ // has to be fixed!
              EXCEPTION("Currently only ONE hysteretic region is allowed!");
            }
            // get nonlinear/hysteretic material parameter
            std::map<std::string, double> ParameterMap;
            actMat->GetScalar(ParameterMap["Ps"], MAG_PS_EB, Global::REAL);
            actMat->GetScalar(ParameterMap["A"], MAG_A_EB, Global::REAL);
            actMat->GetScalar(ParameterMap["mu0"], MAG_MU0_EB, Global::REAL);
            actMat->GetScalar(ParameterMap["numS"], MAG_NUMS_EB, Global::REAL);
            actMat->GetScalar(ParameterMap["chi_factor"], MAG_CHI_FACTOR_EB, Global::REAL);
            ParameterMap["isMH"] = 0;
            matModelCoef_->InitModel(ParameterMap, actSDList);
            // evaluate dbdh = mu
            mu = matModelCoef_;
            nlFluxCoef_->AddRegion(actRegion, matModelCoef_);
            nlScalCoef_->AddRegion(actRegion, matModelCoef_);
            massInt = new BDBInt<>(new IdentityOperator<FeHCurl,3,1,Double>(),mu,1.0,updatedGeo_);
            massInt->SetName("(db/dh N,N): IdentityIntegrator");
          }
        }
        else{ // LINEAR CASE (difference is described above)
          mu = actMat->GetScalCoefFnc(MAG_PERMEABILITY_SCALAR, Global::REAL);
          perm_->AddRegion(actRegion, mu);
          massInt = new BBInt<>(new IdentityOperator<FeHCurl,3,1,Double>(),mu,1.0,updatedGeo_);
          massInt->SetName("(mu N,N): IdentityIntegrator");
        }
        BiLinFormContext * massContext = new BiLinFormContext(massInt, STIFFNESS );
        massContext->SetEntities( actSDList, actSDList );
        massContext->SetFeFunctions( feFunc, feFunc );
        assemble_->AddBiLinearForm( massContext );
        massInts_[actRegion] = massInt; // insert mass integrator to list of defined mass integrators
        // ====================================================================
        // MASS INTEGRATOR [END]
        // ====================================================================

        // ====================================================================
        // STIFFNESS INTEGRATOR [START] (rho_art curlN, curlN)
        // The coefficent rho_art depends on the material case:
        // ~ linear:      rho_art = (1/(std::pow(penaltyParameter,0.01)))*1*mu_regularize*std::pow(10,penaltyParameter/2)
        // ~ nonlinear:   rho_art = (1/penaltyParameter^0.01)*||mu_regularize||_inf*10^(penaltyParameter/2)
        //                since mu in the nonlinear case is a tensor the inf-norm is a good choice for mu_regularize
        // ~ hysteresis: same as in the nonlinear case
        // ====================================================================
        // get material coefficient (in this case: artificial resistivity that regularizes the PDE and is defined as 1/sigma_artificial = mu*10^(n/2)
        // where n is an integer number)
        Double mu_regularize;
        materials_[actRegion]->GetScalar( mu_regularize, MAG_PERMEABILITY_SCALAR, Global::REAL );
        PtrCoefFct rho_art;
        PtrCoefFct rho_art_nl;
        BaseBDBInt* curlcurl = NULL;

        if (nonLinTypes.Find(PERMEABILITY) != -1){ // NONLINEAR/HYSTERESIS CASE
          if (modelName_ == "JilesAthertonModel"){ // just because, nobody wonders that it does not work with this model
            EXCEPTION("Jiles-Atherton model not implemented for MagEdgeHPDE");
          }
          else if (modelName_ == "EBHysteresisModel"){ // (this is the model we want to use!)
            if(moreThan1HystRegion){ // has to be fixed!
              EXCEPTION("Currently only ONE hysteretic region is allowed!");
            }
            moreThan1HystRegion = true;
            rho_art_nl = nlScalCoef_;
            curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), rho_art_nl,1.0, updatedGeo_);
            curlcurl->SetName("(rho_art_nl curlN,CurlN): CurlCurlIntegrator");
          }
        } else{ // LINEAR CASE
          //rho_art = CoefFunction::Generate(mp_, Global::REAL,lexical_cast<std::string>((1/(std::pow(penaltyParameter,0.01)))*1*mu_regularize*std::pow(10,penaltyParameter/2)));
          rho_art = CoefFunction::Generate(mp_, Global::REAL,lexical_cast<std::string>(mu_regularize*std::pow(10,penaltyParameter/2)));
          curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), rho_art,1.0, updatedGeo_);
          curlcurl->SetName("(rho_art curlN,CurlN): CurlCurlIntegrator");
        }

        BiLinFormContext * stiffContext = new BiLinFormContext(curlcurl, STIFFNESS );
        stiffContext->SetEntities( actSDList, actSDList );
        stiffContext->SetFeFunctions( feFunc, feFunc );
        assemble_->AddBiLinearForm( stiffContext );
        // ====================================================================
        // STIFFNESS INTEGRATOR [END]
        // ====================================================================

        // ====================================================================
        // AUXILIARY CURL-CURL INTEGRATOR [START]
        // ====================================================================
        BaseBDBInt* AUX_curlcurl = NULL;
        PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
        AUX_curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), constOne,1.0, updatedGeo_);
        AUX_curlcurl->SetName("(curlN,CurlN): Auxiliary,CurlCurlIntegrator (just for postprocessing: get curlh from h)");
        BiLinFormContext * AUX_stiffContext = new BiLinFormContext(AUX_curlcurl, STIFFNESS );
        AUX_stiffContext->SetEntities( actSDList, actSDList );
        AUX_stiffContext->SetFeFunctions( feFunc, feFunc );
        // Important: Add bdb-integrator to global list, as we need them later
        // for calculation of postprocessing results
        // - curlh
        bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,AUX_curlcurl));
        // ====================================================================
        // AUXILIARY CURL-CURL INTEGRATOR [END]
        // ====================================================================
        }// end for regions
    }


  // **********************************************************
  // Definition of all linearforms that do NOT involve coils
  // linear:    nothing
  // nonlinear: (rho_art curlh, curlN) + (b(h), N)
  // **********************************************************
  void MagEdgeHPDE::DefineRhsLoadIntegrators(){
    // ============= IMPORTANT (START) ================================
    // all linearforms defined here only occur in the nonlinear case.
    // In the linear case only (rho_art j,curlN) is needed and that is
    // defined in DefineCoilIntegrators()! For this reason an "if" is
    // wrapped around all statements, such that this method is empty {}
    // in the linear case.
    // ============= IMPORTANT (END) ==================================
    if(nonLin_ && (modelName_ == "EBHysteresisModel"))
    {
      // get FEFunctions and space
      shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_FIELD_INTENSITY];
      shared_ptr<FeSpace> feSpace_reduced = feFunc->GetFeSpace();

      LinearForm * lin1 = NULL; // (rho_art curlh,curlN), only in the nonlinear case necessary
      LinearForm * lin2 = NULL; // (b(h),N), only in the nonlinear case necessary
      RegionIdType actRegion;

      bool coefUpdateGeo = true;
      bool isHystereticMat = false;

      Double penaltyParameter = 1e-6;
      myParam_->GetValue("penaltyFunctionParameter", penaltyParameter, ParamNode::PASS);

      // iterate over the region (or materials)
      for (UInt iRegion = 0; iRegion < regions_.GetSize(); iRegion++)
      {
        // set current region and material
        actRegion = regions_[iRegion];
        StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];
        std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );

        // Get current region name
        std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

        // create new entity list
        shared_ptr<ElemList> actSDList(new ElemList(ptGrid_));
        actSDList->SetRegion(actRegion);

        PtrCoefFct fluxDensityNL = NULL;
        fluxDensityNL = matModelCoef_;

        // ===============================================================================================
        // lin1: (rho_art curlh,curlN) [START]
        // curlh:   is obtained by the mag. field intensity h from the last Newton iteration.
        // rho_art: regularization parameter that depends on the current scalar permeability mu
        //    -    nonlinear subregion: uses nlScalCoef_ to get the norm of dbdh
        //    -    linear subregion:    uses linear case to get rho_art
        // ===============================================================================================
        // generate the coefFct that is the multiplication of the curlh and rho_art (rho_art*curlh)
        Double mu_regularize;
        materials_[actRegion]->GetScalar( mu_regularize, MAG_PERMEABILITY_SCALAR, Global::REAL );
        PtrCoefFct rho_art;
        PtrCoefFct rho_art_nl;
        PtrCoefFct rho_times_curlh;
        if (nonLinTypes.Find(PERMEABILITY) != -1){ // NONLINEAR CASE, NONLINEAR SUBREGION
          CoefXprVecScalOp temp = CoefXprVecScalOp(mp_, GetCoefFct( MAG_FIELD_INTENSITY_CURL ), nlScalCoef_, CoefXpr::OP_MULT);
          PtrCoefFct rho_times_curlh = CoefFunction::Generate(mp_, Global::REAL, temp);
          lin1 = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3,Double>(),-1.0, rho_times_curlh, volRegions, coefUpdateGeo);
          lin1->SetName("(rho_art_nl curlh,curlN): residual");
        } else{ // NONLINEAR CASE, LINEAR SUBREGION
          rho_art = CoefFunction::Generate(mp_, Global::REAL,lexical_cast<std::string>(mu_regularize*std::pow(10,penaltyParameter/2)));
          CoefXprVecScalOp temp = CoefXprVecScalOp(mp_, GetCoefFct( MAG_FIELD_INTENSITY_CURL ), rho_art, CoefXpr::OP_MULT);
          PtrCoefFct rho_times_curlh = CoefFunction::Generate(mp_, Global::REAL, temp);
          lin1 = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3,Double>(),-1.0, rho_times_curlh, volRegions, coefUpdateGeo);
          lin1->SetName("(rho_art curlh,curlN): residual");
        }
        LinearFormContext *ctx = new LinearFormContext(lin1);
        ctx->SetEntities(actSDList);
        ctx->SetFeFunction(feFunc);
        assemble_->AddLinearForm(ctx);
        // ===============================================================================================
        // lin1: (rho_art curlh,curlN) [END]
        // ===============================================================================================

        // ===============================================================================================
        // lin2: (b(h),N) [START]
        // the b(h) is obtained by the used material model where the input is the mag. field intensity h
        // from the last Newton iteration.
        // ===============================================================================================
        if (nonLinTypes.Find(PERMEABILITY) != -1 && modelName_ != "nonlinearCurve"){
          if (modelName_ == "JilesAthertonModel")
          {
            EXCEPTION("Jiles-Atherton model not implemented for MagneticScalarPotentialPDE");
          }
          else if (modelName_ == "EBHysteresisModel") // NONLINEAR CASE, AND NONLINEAR SUBREGION
          {
            isHystereticMat = true;
            lin2 = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3>(),(-1.0), fluxDensityNL, coefUpdateGeo, false);
          }
          lin2->SetName("(b(h),N): residual");
          lin2->SetSolDependent();
          LinearFormContext *ctx = new LinearFormContext( lin2 );
          ctx->SetEntities( actSDList );
          ctx->SetFeFunction(feFunc);
          assemble_->AddLinearForm(ctx);
        } else { // NONLINEAR CASE, BUT LINEAR SUBREGION
          lin2 = new BUIntegrator<Double>(new IdentityOperator<FeHCurl, 3>(),(-1.0), GetCoefFct( MAG_FLUX_DENSITY ), volRegions, coefUpdateGeo);
          lin2->SetName("(b(h)_linear,N): residual");
          LinearFormContext *ctx = new LinearFormContext(lin2);
          ctx->SetEntities(actSDList);
          ctx->SetFeFunction(feFunc);
          assemble_->AddLinearForm(ctx);
        }
        // ===============================================================================================
        // lin2: (b(h),N) [END]
        // ===============================================================================================
      }// end loop over entities
    } 
  }


  // **********************************************************
  // definition of the primary result 
  // magnetic field intensity h
  // **********************************************************
  void MagEdgeHPDE::DefinePrimaryResults() {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    // =====================================================
    // MAG_FIELD_INTENSITY (START)
    // =====================================================
    shared_ptr<ResultInfo> potInfo(new ResultInfo);
    potInfo->resultType = MAG_FIELD_INTENSITY;
    potInfo->SetFeFunction(feFunctions_[MAG_FIELD_INTENSITY]);
    potInfo->dofNames = vecComponents;
    potInfo->unit = "A/m";
    potInfo->definedOn = ResultInfo::ELEMENT;
    potInfo->entryType = ResultInfo::VECTOR;
    feFunctions_[MAG_FIELD_INTENSITY]->SetResultInfo(potInfo);
    DefineFieldResult( feFunctions_[MAG_FIELD_INTENSITY], potInfo );
    // =====================================================
    // MAG_FIELD_INTENSITY (END)
    // =====================================================

  }


  // **********************************************************
  // definition of postprocessed results 
  // - magnetic flux density:  b = mu h
  // - source current density: j = curlh
  // **********************************************************
  void MagEdgeHPDE::DefinePostProcResults() {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";
    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_FIELD_INTENSITY];

    // =====================================================
    // MAG_FLUX_DENSTIY (START)
    // =====================================================
    shared_ptr<ResultInfo> magFlux(new ResultInfo);
    magFlux->resultType = MAG_FLUX_DENSITY;
    magFlux->SetVectorDOFs(dim_, isaxi_);
    magFlux->dofNames = vecComponents;
    magFlux->unit = "Vs/m^2";
    magFlux->definedOn = ResultInfo::ELEMENT;
    magFlux->entryType = ResultInfo::VECTOR;
    // The recipe about how to actually evaluate B, is defined in FinalizePostProcResults()
    shared_ptr<CoefFunctionMulti> bFunc(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
    DefineFieldResult(bFunc, magFlux);
    // =====================================================
    // MAG_FLUX_DENSTIY (END)
    // =====================================================

    // =====================================================
    // MAG_FIELD_INTENSITY_CURL (START)
    // =====================================================
    shared_ptr<ResultInfo> magFieldIntensCurl(new ResultInfo);
    magFieldIntensCurl->resultType = MAG_FIELD_INTENSITY_CURL;
    magFieldIntensCurl->SetVectorDOFs(dim_, isaxi_);
    magFieldIntensCurl->dofNames = vecComponents;
    magFieldIntensCurl->unit = "A/m^2";
    magFieldIntensCurl->definedOn = ResultInfo::ELEMENT;
    magFieldIntensCurl->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionFormBased> curlhFunc;
    curlhFunc.reset(new CoefFunctionBOp<Double>(feFunc, magFieldIntensCurl, 1.0));
    DefineFieldResult(curlhFunc, magFieldIntensCurl);
    stiffFormCoefs_.insert(curlhFunc);
    // =====================================================
    // MAG_FIELD_INTENSITY_CURL (END)
    // =====================================================

  }

  // **********************************************************
  // recipe to do the postprocessing for certain quantities:
  // MAG_FLUX_DENSTIY
  // MAG_FIELD_INTENSITY_CURL
  // **********************************************************
  void MagEdgeHPDE::FinalizePostProcResults() {
    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();

    shared_ptr<CoefFunctionMulti> bCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_FLUX_DENSITY]);

    StdVector<RegionIdType>::iterator regIt = regions_.Begin();
    regIt = regions_.Begin();

    for (; regIt != regions_.End(); ++regIt){
      // =====================================================
      // MAG_FLUX_DENSTIY (START)
      // =====================================================
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[*regIt];  // Just to find out which linear/nonlinear type is defined in this region
      if( nonLinTypes.Find(PERMEABILITY) != -1 && modelName_ != "nonlinearCurve" ){
        // hysteretic/nonlinear case
        bCoef->AddRegion(*regIt, nlFluxCoef_);
      }
      else{
        // classical nonlinear case and linear case
        PtrCoefFct b = CoefFunction::Generate( mp_, Global::REAL, CoefXprVecScalOp(mp_, GetCoefFct( MAG_FIELD_INTENSITY ), perm_, CoefXpr::OP_MULT));
        bCoef->AddRegion(*regIt, b);
      }
      // =====================================================
      // MAG_FLUX_DENSTIY (END)
      // =====================================================
    }
  }




  std::map<SolutionType, shared_ptr<FeSpace> >
  MagEdgeHPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode ) {
      //create grid based approximation Hcurl elements and standard Gauss integration
      std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
      if(formulation == "default" || formulation == "H_CURL"){
      PtrParamNode potSpaceNode = infoNode->Get("magFieldIntensity");
      crSpaces[MAG_FIELD_INTENSITY] =
          FeSpace::CreateInstance(myParam_, potSpaceNode, FeSpace::HCURL, ptGrid_ );
      crSpaces[MAG_FIELD_INTENSITY]->Init(solStrat_);
      }else{
      EXCEPTION("The formulation " << formulation
                  << "of magnetic edge PDE is not known!");
      }
      return crSpaces;
  }


} // end of namespace
