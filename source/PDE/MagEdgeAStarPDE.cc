#include <fstream>

#include "MagEdgeAStarPDE.hh"

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
  DEFINE_LOG(magEdgeAStarPde, "magEdgeAStarPde")


  // **************
  //  Constructor
  // **************
  MagEdgeAStarPDE::MagEdgeAStarPDE( Grid * aptgrid, PtrParamNode paramNode,
                          PtrParamNode infoNode,
                          shared_ptr<SimState> simState, Domain* domain )
    :MagBasePDE( aptgrid, paramNode, infoNode, simState, domain ) {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magneticAStarEdge";
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
    matModelCoefm_.clear();
    // init the nlFieldIntensCoefm_
    nlFieldIntensCoefm_.clear();
  }


  // **********************************************************
  //  Destructor
  // **********************************************************
  MagEdgeAStarPDE::~MagEdgeAStarPDE(){
  }


  // **********************************************************
  //  Initialize Nonlinearities
  // **********************************************************
  void MagEdgeAStarPDE::InitNonLin(){
    SinglePDE::InitNonLin();
  }


  // **********************************************************
  // Define SolveStep depending on the material law
  // nonlinear/hysteresis: SolveStepEB
  // linear:               StdSolveStep
  // **********************************************************
  void MagEdgeAStarPDE::DefineSolveStep(){


    UInt is_pseudo_time_stepping = 1; // per default pseudo time stepping is used
    myParam_->GetValue("isPseudoTimeStepping", is_pseudo_time_stepping, ParamNode::PASS);
    if(nonLin_ && (modelName_ == "invEBHysteresisModel")){
      solveStep_ = new SolveStepEB(*this, is_pseudo_time_stepping);
    } 
    else{
      solveStep_ = new StdSolveStep(*this);
    }
  }


  // **********************************************************
  // INIT. TIME STEPPING (not a real time stepping, since the
  // resulting FEM system is not an ODE)
  // **********************************************************
  void MagEdgeAStarPDE::InitTimeStepping(){

    Double gamma = 1.0;
    GLMScheme *scheme = new Trapezoidal(gamma);

    TimeSchemeGLM::NonLinType nlType = (nonLin_ || isHysteresis_) ? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType));
    feFunctions_[MAG_POTENTIAL]->SetTimeScheme(myScheme);
  }


  // **********************************************************
  //  Definition of Integrators
  // **********************************************************
  void MagEdgeAStarPDE::DefineIntegrators(){
    this->DefineStandardIntegrators();
    DefineCoilIntegrators(); // in MagBasePDE
  }


  // **********************************************************
  // just a helper function for DefineCoilIntegrators()
  // **********************************************************
  LinearForm* MagEdgeAStarPDE::GetCurrentDensityInt( Double factor, PtrCoefFct coef, std::string coilId){
    LinearForm * ret = NULL;
    if (dim_ == 3) {
      ret = new BUIntegrator<Double>(new IdentityOperator<FeHCurl, 3, 1, Double>(), factor, coef, updatedGeo_);
    } else {
      EXCEPTION("AStar is only available in 3D!")
    }
    return ret;
  }


  // *************************************************************
  // Definition of linearforms that involve coils
  // linear:    (p j,curlN)
  // nonlinear: - (p j,curlN); negative due to Newton scheme
  // *************************************************************
  void MagEdgeAStarPDE::DefineCoilIntegrators(){

    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;

    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_POTENTIAL];
    shared_ptr<FeSpace> feSpace = feFunc->GetFeSpace();

    std::map<Coil::IdType, shared_ptr<Coil>>::iterator coilIt;
    coilIt = coils_.begin();
    for (; coilIt != coils_.end(); coilIt++){
      Coil &actCoil = *(coilIt->second);
      // run over all parts
      std::map<RegionIdType, shared_ptr<Coil::Part>>::iterator partIt;
      partIt = actCoil.parts_.begin();

      if (actCoil.sourceType_ != Coil::CURRENT){
        EXCEPTION("Only current excitation is implemented in AStar-formulation!!");
      } 
      else{

        for (partIt = actCoil.parts_.begin(); partIt != actCoil.parts_.end(); partIt++){
          Coil::Part &actPart = *(partIt->second);
          RegionIdType actRegion = partIt->first;
          BaseMaterial * actSDMat = NULL;
          actSDMat    = materials_[actRegion];
          shared_ptr<ElemList> actSDList(new ElemList(ptGrid_));
          actSDList->SetRegion(actRegion);
          LinearForm* curInt = NULL;

          std::map<UInt, PtrCoefFct> jFct;
          CoefXprVecScalOp iVec = CoefXprVecScalOp(mp_, actPart.jUnitVec, actCoil.srcVal_, CoefXpr::OP_MULT);
          PtrCoefFct iFct = CoefFunction::Generate(mp_, part, iVec);

          // ============= IMPORTANT (START) =======================================
          // linear: (J,A')
          // nonlinear: (J,A')
          // Since both material cases are the same there is no difference.
          // Also for real time stepping there is no difference.
          // ============= IMPORTANT (END) ======================================= 
          CoefXprVecScalOp jVec = CoefXprVecScalOp(mp_, iFct, boost::lexical_cast<std::string>(actPart.wireCrossSect), CoefXpr::OP_DIV);
          jFct[0] = CoefFunction::Generate(mp_, part, jVec);
          coilCurrentDens_[actRegion] = jFct[0];
          curInt = GetCurrentDensityInt( 1.0, jFct[0] );
          curInt->SetName("(J,A'): CoilIntegrator");
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
  // linear:       (nu curlA, curlA') + (sigma A, A')
  // nonlinear: (dH/dB curlA, curlA') + (sigma A, A')
  // **********************************************************
  void MagEdgeAStarPDE::DefineStandardIntegrators(){

    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;

    // global storage for reluctivity makes the postprocessing easier
    reluc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, false));
    conductivity_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, false));
    
    // get the beta from the .xml file and set the default value
    Double beta = 1e-6;
    myParam_->GetValue("penaltyFunctionParameter", beta, ParamNode::PASS);

    // Get the flag for pseudo time stepping
    // 1 ... pseudo time stepping ( simulate more time steps, but no time derivatives are involved )
    // 0 ... real time stepping (time derivatives are involved) 
    UInt is_pseudo_time_stepping = 1; // per default pseudo time stepping is used
    myParam_->GetValue("isPseudoTimeStepping", is_pseudo_time_stepping, ParamNode::PASS);

    // Get dt
    Double dt = 1e-3;
    myParam_->GetValue("delta_t", dt, ParamNode::PASS);

    // get FEFunction and space
    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_POTENTIAL];
    shared_ptr<FeSpace> feSpace = feFunc->GetFeSpace();

    PtrCoefFct magFluxCoef = this->GetCoefFct(MAG_FLUX_DENSITY);

    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion ++){
      // set current region and materials
      actRegion = regions_[iRegion];
      actSDMat    = materials_[actRegion];
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

      // Init material model for hysteretic transient analysis
      if (((analysistype_ == STATIC) || (analysistype_ == TRANSIENT)) && nonLin_ && (modelName_ != "nonlinearCurve")){
        matModelCoefm_[actRegion].reset(new CoefFunctionMaterialModel<Complex>());
        matModelCoefm_[actRegion]->Init(magFluxCoef, modelName_, dim_); 
      }

      // ===================================================================
      // PSEUDO TIME-STEPPING [START]
      // ===================================================================
      if (is_pseudo_time_stepping == 1){
        // ====================================================================
        // STIFFNESS INTEGRATOR [START] (c curlA, curlA')
        // The coefficient c of this integrator depends on the material case:
        // ~ linear:     use just the value of the reluctivity of this region.
        // ~ nonlinear:  use the derivative of the material law b = b(h) as the
        //             value of the coefficient, since this is the Jacobian J 
        //             of the Newton-Raphson scheme J = db/dh.
        // ~ hysteresis: same as in the nonlinear case. However, typically the 
        //             Jacobian can not be determined analytically, but by 
        //             a numerical method like finite differences or Broyden.
        // So here the only thing that is "if cased", is the determination of
        // coefficient and the integrator that accepts also tensors as 
        // material coefficent. The integrator (structure, (c N, N)) itself
        // stays the same!
        // ====================================================================
        PtrCoefFct nuNL = NULL;
        PtrCoefFct nu = NULL;
        PtrCoefFct mu = NULL;
        BaseBDBInt *stiffInt = NULL;
        if (nonLinTypes.Find(PERMEABILITY) != -1){ // NONLINEAR/HYSTERESIS CASE
          if (modelName_ == "JilesAthertonModel"){ // just because, nobody wonders that it does not work with this model
            EXCEPTION("Jiles-Atherton model not implemented for MagEdgeAStarPDE");
          } else if (modelName_ == "invEBHysteresisModel") { // (this is the model we want to use!)
            // init. inverse EB Material Model
            std::map<std::string, double> ParameterMap;
            std::map<std::string, std::string> StringParameterMap;
            if(actSDMat->GetAnhystMagModel() == "analytic_anhysteresis"){
              actSDMat->GetString(StringParameterMap["anhyst_type"], MAG_ANHYST_TYPE_INVEB);
              if(actSDMat->GetAnhystFormula() == "tan"){
                actSDMat->GetString(StringParameterMap["anhyst_formula"], MAG_ANHYST_FORMULA_INVEB);
                actSDMat->GetScalar(ParameterMap["Js"], MAG_JS_INVEB, Global::REAL);
                actSDMat->GetScalar(ParameterMap["A"], MAG_A_INVEB, Global::REAL);
              } else if(actSDMat->GetAnhystFormula() == "atan"){
                actSDMat->GetString(StringParameterMap["anhyst_formula"], MAG_ANHYST_FORMULA_INVEB);
                actSDMat->GetScalar(ParameterMap["Ms"], MAG_MS_INVEB, Global::REAL);
                actSDMat->GetScalar(ParameterMap["A"], MAG_A_INVEB, Global::REAL);
              } else if(actSDMat->GetAnhystFormula() == "pacejka"){
                actSDMat->GetString(StringParameterMap["anhyst_formula"], MAG_ANHYST_FORMULA_INVEB);
                actSDMat->GetScalar(ParameterMap["Ms"], MAG_MS_INVEB, Global::REAL);
                actSDMat->GetScalar(ParameterMap["pa"], MAG_PA_INVEB, Global::REAL);
                actSDMat->GetScalar(ParameterMap["pb"], MAG_PB_INVEB, Global::REAL);
                actSDMat->GetScalar(ParameterMap["pc"], MAG_PC_INVEB, Global::REAL);
              } else if (actSDMat->GetAnhystFormula() == "brauer") {
                actSDMat->GetString(StringParameterMap["anhyst_formula"], MAG_ANHYST_FORMULA_INVEB);
                actSDMat->GetScalar(ParameterMap["p_0"], MAG_P0_INVEB, Global::REAL);
                actSDMat->GetScalar(ParameterMap["p_1"], MAG_P1_INVEB, Global::REAL);
                actSDMat->GetScalar(ParameterMap["p_2"], MAG_P2_INVEB, Global::REAL);
              } else if (actSDMat->GetAnhystFormula() == "lookuptable") {
                actSDMat->GetString(StringParameterMap["anhyst_formula"], MAG_ANHYST_FORMULA_INVEB);
                actSDMat->GetString(StringParameterMap["lookup_table_file"], MAG_LOOKUP_TABLE_FILE_INVEB);
              }
            }
            actSDMat->GetString(StringParameterMap["weights_file_path"], MAG_WEIGHTS_FILE_PATH_EB);
            actSDMat->GetString(StringParameterMap["pinning_forces_weights_file"], MAG_PINNING_FORCES_WEIGHTS_INVEB);
            actSDMat->GetScalar(ParameterMap["approx_type"], MAG_APPROX_TYPE, Global::REAL);
            ParameterMap["isMH"] = 0;
            actSDMat->GetScalar(ParameterMap["jacobian_method"], MAG_JACOBIAN_METHOD_INVEB, Global::REAL);
            matModelCoefm_[actRegion]->InitModel(ParameterMap, StringParameterMap, actSDList);
            nuNL = matModelCoefm_[actRegion];

            nlFieldIntensCoefm_[actRegion].reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_, true)); 
            nlFieldIntensCoefm_[actRegion]->AddRegion(actRegion, matModelCoefm_[actRegion]);
            if (dim_ == 3) {
              stiffInt = new BDBInt<>(new CurlOperator<FeHCurl, 3, Double>(), nuNL, 1.0, updatedGeo_);
            } else {
              EXCEPTION("AStar is only available in 3D!");
            }
            stiffInt->SetName("(dh/db curlA, curlA'): CurlIntegrator");
          }
        } else { // LINEAR CASE (difference is described above)
          nu = actSDMat->GetScalCoefFnc(MAG_RELUCTIVITY_SCALAR, Global::REAL);
          reluc_->AddRegion(actRegion, nu);
          if (dim_ == 3) {
            stiffInt = new BBInt<>(new CurlOperator<FeHCurl, 3, Double>(), nu, 1.0, updatedGeo_);
          } else {
            EXCEPTION("AStar is only available in 3D!");
          }
          stiffInt->SetName("(nu curlA, curlA'): CurlIntegrator");
        }
        BiLinFormContext * stiffContext = new BiLinFormContext(stiffInt, STIFFNESS );
        stiffContext->SetEntities( actSDList, actSDList );
        stiffContext->SetFeFunctions( feFunc, feFunc );
        assemble_->AddBiLinearForm( stiffContext );
        bdbInts_.insert(std::pair<RegionIdType, BaseBDBInt *>(actRegion, stiffInt));
        // ====================================================================
        // STIFFNESS INTEGRATOR [END]
        // ====================================================================

        // ====================================================================
        // MASS INTEGRATOR [START] (epsilon N, N)
        // In the pseudo time stepping case the \epsilon is not a real conductivity,
        // but an artificial one that regularizes the stiffness matrix. In fact,
        // it is a regularization parameter.
        // ====================================================================
        Double mu_regularize;
        PtrCoefFct epsilon;
        BaseBDBInt* massInt = NULL;
        // get material coefficient (in this case: artificial conductivity that regularizes the PDE and is defined 
        // as epsilon = regularization parameter*mu;
        materials_[actRegion]->GetScalar( mu_regularize, MAG_PERMEABILITY_SCALAR, Global::REAL );
        epsilon = CoefFunction::Generate(mp_, Global::REAL,lexical_cast<std::string>(mu_regularize*beta));
        if (dim_ == 3) {
          //massInt = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<FeHCurl, 3, Double>(),epsilon, 1.0);
          massInt = new BBIntMassEdge<>(new IdentityOperator<FeHCurl, 3, 1, Double>(),epsilon, 1.0);
        } else {
          EXCEPTION("AStar is only available in 3D!");
        }
        massInt->SetName("(epsilon A, A): ScaledByEdgeIdentityOperator");
        BiLinFormContext * massContext = new BiLinFormContext(massInt, STIFFNESS );
        massContext->SetEntities( actSDList, actSDList );
        massContext->SetFeFunctions( feFunc, feFunc );
        assemble_->AddBiLinearForm( massContext );
        massInts_[actRegion] = massInt;
        // ====================================================================
        // MASS INTEGRATOR [END]
        // ====================================================================
      }
      // ===================================================================
      // PSEUDO TIME-STEPPING [END]
      // ===================================================================

      // ===================================================================
      // REAL TIME-STEPPING [START]
      // ===================================================================
      if (is_pseudo_time_stepping == 0){
          EXCEPTION("REAL TIME STEPPING NOT IMPLENTED YET!")
          // FOR TIME STEPPING
      }
      // ===================================================================
      // REAL TIME-STEPPING [END]
      // ===================================================================
    }// end for regions
  }

  // **********************************************************
  // Definition of all linearforms that do NOT involve coils
  // linear:    nothing
  // nonlinear: (p curlh, curlN) + (b(h), N)
  // **********************************************************
  void MagEdgeAStarPDE::DefineRhsLoadIntegrators(){
    // ============= IMPORTANT (START) ================================
    // all linearforms defined here only occur in the nonlinear case.
    // In the linear case only (J,A') is needed and that is
    // defined in DefineCoilIntegrators()! For this reason an "if" is
    // wrapped around all statements, such that this method is empty {}
    // in the linear case.
    // ============= IMPORTANT (END) ==================================
    
    // get FEFunctions and space
    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_POTENTIAL];
    shared_ptr<FeSpace> feSpace_reduced = feFunc->GetFeSpace();

    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    StdVector<std::string> vecDofNames = feFunc->GetResultInfo()->dofNames;
    LinearForm * lin1_pts = NULL; // (h(b),curlA'), only in the nonlinear case necessary (PSEUDO TIME STEPPING (pts))
    LinearForm * lin2_pts = NULL; // (epsilon A,A'), only in the nonlinear case necessary (PSEUDO TIME STEPPING (pts))
    LinearForm * lin1_rts = NULL; // (rho curlh,curlN), only in the nonlinear case necessary (REAL TIME STEPPING (rts))
    LinearForm * lin2_rts = NULL; // ((1/delta_t)*b(h)_n,N), only in the nonlinear case necessary (REAL TIME STEPPING (rts))
    LinearForm * lin3_rts = NULL; // ((1/delta_t)*b(h)_{n-1},N), only in the nonlinear case necessary (REAL TIME STEPPING (rts))
    LinearForm * lin_br = NULL; // (b_r,N), describes linear permanent magnets in the linear/nonlinear case
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;

    bool coefUpdateGeo = true;
    bool isHystereticMat = false;

    Double beta = 1e-6;
    myParam_->GetValue("penaltyFunctionParameter", beta, ParamNode::PASS);

    // Get dt
    Double dt = 1e-3;
    myParam_->GetValue("delta_t", dt, ParamNode::PASS);

    // Get the flag for pseudo time stepping
    // 1 ... pseudo time stepping ( simulate more time steps, but no time derivatives are involved )
    // 0 ... real time stepping (time derivatives are involved) 
    UInt is_pseudo_time_stepping = 1; // per default pseudo time stepping is used
    myParam_->GetValue("isPseudoTimeStepping", is_pseudo_time_stepping, ParamNode::PASS);

    // iterate over the region (or materials)
    for (UInt iRegion = 0; iRegion < regions_.GetSize(); iRegion++) {
      // set current region and materials
      actRegion = regions_[iRegion];
      actSDMat    = materials_[actRegion];

      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];
      std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );

      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList(new ElemList(ptGrid_));
      actSDList->SetRegion(actRegion);

      // ===================================================================
      // PSEUDO TIME-STEPPING [START]
      // ===================================================================
      if (is_pseudo_time_stepping == 1){
        if(nonLin_ && (modelName_ == "invEBHysteresisModel")) {
          // ===============================================================================================
          // lin1_pts: (h(b),curlA') [START]
          // h(b):   is obtained by the mag. field intensity h from the last Newton iteration.
          // ===============================================================================================
          PtrCoefFct fieldIntensityNL = NULL;
          fieldIntensityNL = matModelCoefm_[actRegion];
          if (nonLinTypes.Find(PERMEABILITY) != -1) { // NONLINEAR CASE, NONLINEAR SUBREGION
            if (dim_ == 3) {
              lin1_pts = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3,Double>(),-1.0, fieldIntensityNL, volRegions, coefUpdateGeo);
            } else {
              EXCEPTION("AStar is only available in 3D!");
            }
            lin1_pts->SetName("(h(b), curlA'): residual");
            LinearFormContext *ctx = new LinearFormContext(lin1_pts);
            ctx->SetEntities(actSDList);
            ctx->SetFeFunction(feFunc);
            assemble_->AddLinearForm(ctx);
          } else { // NONLINEAR CASE, LINEAR SUBREGION
            if (dim_ == 3) {
              lin1_pts = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3,Double>(),-1.0, GetCoefFct( MAG_FIELD_INTENSITY ), volRegions, coefUpdateGeo);
            } else {
              EXCEPTION("AStar is only available in 3D!");
            }
            lin1_pts->SetName("(h, curlA'): residual");
            LinearFormContext *ctx = new LinearFormContext(lin1_pts);
            ctx->SetEntities(actSDList);
            ctx->SetFeFunction(feFunc);
            assemble_->AddLinearForm(ctx);
          }
          // ===============================================================================================
          // lin1_pts: (h(b),curlA') [END]
          // ===============================================================================================

          // ===============================================================================================
          // lin2_pts: (epsilon A, A') [START]
          // the b(h) is obtained by the used material model where the input is the mag. field intensity h
          // from the last Newton iteration.
          // ===============================================================================================
          Double mu_regularize;
          PtrCoefFct epsilon;
          BaseBDBInt* massInt = NULL;
          // get material coefficient (in this case: artificial conductivity that regularizes the PDE and is defined 
          // as epsilon = regularization parameter*mu;
          materials_[actRegion]->GetScalar( mu_regularize, MAG_PERMEABILITY_SCALAR, Global::REAL );
          epsilon = CoefFunction::Generate(mp_, Global::REAL,lexical_cast<std::string>(mu_regularize*beta));
          CoefXprVecScalOp temp = CoefXprVecScalOp(mp_, GetCoefFct( MAG_POTENTIAL ), epsilon, CoefXpr::OP_MULT);
          PtrCoefFct epsilon_times_A = CoefFunction::Generate(mp_, Global::REAL, temp);
          if (dim_ == 3) {
            lin2_pts = new BUIntegrator<Double>(new IdentityOperator<FeHCurl, 3>(),(-1.0), epsilon_times_A, volRegions, coefUpdateGeo); 
          } else {
            EXCEPTION("AStar is only available in 3D!");
          }
          lin2_pts->SetName("(epsilon A, A'): residual");
          LinearFormContext *ctx = new LinearFormContext(lin2_pts);
          ctx->SetEntities(actSDList);
          ctx->SetFeFunction(feFunc);
          assemble_->AddLinearForm(ctx);
          // ===============================================================================================
          // lin2_pts: (b(h),N) [END]
          // ===============================================================================================
        }
      } 
      // ===================================================================
      // PSEUDO TIME-STEPPING [END]
      // ===================================================================

      // ===================================================================
      // REAL TIME-STEPPING [START]
      // ===================================================================     
      if (is_pseudo_time_stepping == 0){
          EXCEPTION("REAL TIME STEPPING NOT IMPLENTED YET!")
          // FOR TIME STEPPING
      } 
      // ===================================================================
      // REAL TIME-STEPPING [END]
      // =================================================================== 
    }// end loop over entities
  }


  // **********************************************************
  // definition of the primary result 
  // magnetic field intensity h
  // **********************************************************
  void MagEdgeAStarPDE::DefinePrimaryResults() {

    StdVector<std::string> vecComponents;
    if (dim_ == 3) {
      vecComponents = "x", "y", "z";
    } else {
      EXCEPTION("AStar is only available in 3D!");
    }
    
    // =====================================================
    // MAG_POTENTIAL (START)
    // =====================================================
    shared_ptr<ResultInfo> potInfo(new ResultInfo);
    potInfo->resultType = MAG_POTENTIAL;
    potInfo->dofNames = vecComponents;
    potInfo->unit = "Vs/m";
    potInfo->definedOn = ResultInfo::ELEMENT;
    potInfo->entryType = ResultInfo::VECTOR;
    feFunctions_[MAG_POTENTIAL]->SetResultInfo(potInfo);
    DefineFieldResult(feFunctions_[MAG_POTENTIAL], potInfo);
    // =====================================================
    // MAG_POTENTIAL (END)
    // =====================================================

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MAG_POTENTIAL] = "fluxParallel";
    idbcSolNameMap_[MAG_POTENTIAL] = "potential";
  }


  // **********************************************************
  // definition of postprocessed results 
  // - magnetic flux density:  b = curl a
  // - source current density: h = nu b
  // **********************************************************
  void MagEdgeAStarPDE::DefinePostProcResults() {

    StdVector<std::string> vecComponents;
    if (dim_ == 3) {
      vecComponents = "x", "y", "z";
    } else {
      EXCEPTION("AStar is only available in 3D!");
    }
    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_POTENTIAL];

    // =====================================================
    // MAG_FLUX_DENSTIY (START)
    // =====================================================
    shared_ptr<ResultInfo> magFluxDensity(new ResultInfo);
    magFluxDensity->resultType = MAG_FLUX_DENSITY;
    magFluxDensity->SetVectorDOFs(dim_, isaxi_);
    magFluxDensity->dofNames = vecComponents;
    magFluxDensity->unit = "T";
    magFluxDensity->definedOn = ResultInfo::ELEMENT;
    magFluxDensity->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionFormBased> bFunc;
    bFunc.reset(new CoefFunctionBOp<Double>(feFunc, magFluxDensity, 1.0));
    DefineFieldResult(bFunc, magFluxDensity);
    stiffFormCoefs_.insert(bFunc);
    // =====================================================
    // MAG_FLUX_DENSTIY (END)
    // =====================================================

    // =====================================================
    // MAG_FIELD_INTENSITY (START)
    // =====================================================
    shared_ptr<ResultInfo> magIntens(new ResultInfo);
    magIntens->resultType = MAG_FIELD_INTENSITY;
    magIntens->SetVectorDOFs(dim_, isaxi_);
    magIntens->dofNames = vecComponents;
    magIntens->unit = "A/m";
    magIntens->definedOn = ResultInfo::ELEMENT;
    magIntens->entryType = ResultInfo::VECTOR;
    // The recipe how to actually evaluate H, is defined in FinalizePostProcResults()
    shared_ptr<CoefFunctionMulti> hFunc(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
    DefineFieldResult(hFunc, magIntens);     
    // =====================================================
    // MAG_FIELD_INTENSITY (END)
    // =====================================================
  }

  // **********************************************************
  // recipe to do the postprocessing for certain quantities:
  // MAG_FIELD_INTENSITY
  // **********************************************************
  void MagEdgeAStarPDE::FinalizePostProcResults() {
    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();

    shared_ptr<CoefFunctionMulti> hCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_FIELD_INTENSITY]);
    shared_ptr<CoefFunctionMulti> eCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[ELEC_FIELD_INTENSITY]);
    shared_ptr<CoefFunctionMulti> iCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_TOTAL_CURRENT_DENSITY]);
    StdVector<RegionIdType>::iterator regIt = regions_.Begin();
    regIt = regions_.Begin();

    // =====================================================
    // MAG_FIELD_INTENSITY (START) 
    // H = nu*B (linear case)
    // H = H(B) (nonlinear/hysteresis case)
    // =====================================================
    for (; regIt != regions_.End(); ++regIt){
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[*regIt];  // Just to find out which linear/nonlinear type is defined in this region
      PtrCoefFct h;
      if( nonLinTypes.Find(PERMEABILITY) != -1 && modelName_ != "nonlinearCurve" ){
        std::cout << "here FinalizePostProcResults - nonlin\n";
        h = nlFieldIntensCoefm_[*regIt];
        hCoef->AddRegion(*regIt, h);
      } else {
        std::cout << "here FinalizePostProcResults - lin\n";
        h = CoefFunction::Generate( mp_, Global::REAL, CoefXprVecScalOp(mp_, GetCoefFct( MAG_FLUX_DENSITY ), reluc_, CoefXpr::OP_MULT));
        hCoef->AddRegion(*regIt, h);
      }
    }
    // =====================================================
    // MAG_FIELD_INTENSITY (END) 
    // =====================================================
  }

  std::map<SolutionType, shared_ptr<FeSpace> >
  MagEdgeAStarPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode ) {
      //create grid based approximation Hcurl elements and standard Gauss integration
      std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
      if(formulation == "default" || formulation == "H_CURL"){
      PtrParamNode potSpaceNode = infoNode->Get("magPotential");
      crSpaces[MAG_POTENTIAL] = FeSpace::CreateInstance(myParam_, potSpaceNode, FeSpace::HCURL, ptGrid_ );
      crSpaces[MAG_POTENTIAL]->Init(solStrat_);
      }else{
      EXCEPTION("The formulation " << formulation << "of magnetic edge PDE is not known!");
      }
      return crSpaces;
  }
} // end of namespace
