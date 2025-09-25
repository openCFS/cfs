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
    matModelCoefm_.clear();
    // init the nlFluxCoef
    nlFluxCoefm_.clear();
    nlScalCoefm_.clear();
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

    UInt is_pseudo_time_stepping = 1; // per default pseudo time stepping is used
    myParam_->GetValue("isPseudoTimeStepping", is_pseudo_time_stepping, ParamNode::PASS);
    if(nonLin_ && (modelName_ == "EBHysteresisModel")){
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
    if (dim_ == 2) {
      ret = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 2, Double>(), factor, coef, updatedGeo_);
    } else {
      ret = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3, Double>(), factor, coef, updatedGeo_);
    }
    
    return ret;
  }


  // *************************************************************
  // Definition of linearforms that involve coils
  // linear:    (p j,curlN)
  // nonlinear: - (p j,curlN); negative due to Newton scheme
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
          BaseMaterial * actSDMat = NULL;
          actSDMat    = materials_[actRegion];
          shared_ptr<ElemList> actSDList(new ElemList(ptGrid_));
          actSDList->SetRegion(actRegion);
          LinearForm* curInt = NULL;

          std::map<UInt, PtrCoefFct> jFct;
          CoefXprVecScalOp iVec = CoefXprVecScalOp(mp_, actPart.jUnitVec, actCoil.srcVal_, CoefXpr::OP_MULT);
          PtrCoefFct iFct = CoefFunction::Generate(mp_, part, iVec);

          // ============= IMPORTANT (START) =======================================
          // Here the value for the current density vector j is determined.
          // The linearform of the formulation is as follows
          //         (p j, curl(h')) ,
          // where p is a regularization parameter that occurs due to the
          // penalty method to end up in this formulation.
          // A common way to calculate this coefficient is
          //         p = mu_rmax*mu*10^(beta/2)*1/(beta^0.01),
          //         beta ... numerical parameter
          // where mu is the permeability in this region and n is a penalty parameter
          // that should be chosen rather small. However, since coils are usually 
          // in regions where mu ~ mu0 this factor is harcoded and only beta can be 
          // altered.  Anyway, this should work also for the nonlinear hysteretic case,
          // since the coil region is usually linear material.
          // (Honestly, I did not know how to implement it differently)
          // ============= IMPORTANT (END) ======================================= 
          Double beta = 1; // default value for beta
          myParam_->GetValue("penaltyFunctionParameter", beta, ParamNode::PASS); // get the beta from the .xml file and set the default value
          // Get the flag for pseudo time stepping
          // 1 ... pseudo time stepping ( simulate more time steps, but no time derivatives are involved )
          // 0 ... real time stepping (time derivatives are involved) 
          UInt is_pseudo_time_stepping = 1; // per default pseudo time stepping is used
          myParam_->GetValue("isPseudoTimeStepping", is_pseudo_time_stepping, ParamNode::PASS);
          // ===================================================================
          // PSEUDO TIME-STEPPING [START]
          // ===================================================================
          if (is_pseudo_time_stepping == 1){   
            //CoefXprVecScalOp jVec = CoefXprVecScalOp(mp_, iFct, boost::lexical_cast<std::string>(actPart.wireCrossSect*(1/((1/(std::pow(beta,0.01)))*1*1.256637061e-06*std::pow(10,beta/2)))), CoefXpr::OP_DIV);
            CoefXprVecScalOp jVec = CoefXprVecScalOp(mp_, iFct, boost::lexical_cast<std::string>(actPart.wireCrossSect), CoefXpr::OP_DIV);
            PtrCoefFct jVec_coefFct = CoefFunction::Generate(mp_, part, jVec);
            CoefXprVecScalOp p_times_jVec = CoefXprVecScalOp(mp_, jVec_coefFct, boost::lexical_cast<std::string>(1.256637061e-06*std::pow(10,beta/2)), CoefXpr::OP_MULT);
            jFct[0] = CoefFunction::Generate(mp_, part, p_times_jVec);
            coilCurrentDens_[actRegion] = jFct[0];
            curInt = GetCurrentDensityInt( 1.0, jFct[0] );
            curInt->SetName("(p j,curlN): CoilIntegrator");
            LinearFormContext * coilContext = new LinearFormContext( curInt );
            coilContext->SetEntities( actSDList );
            coilContext->SetFeFunction( feFunc );
            assemble_->AddLinearForm( coilContext );
          }
          else if (is_pseudo_time_stepping == 0){   
            if(nonLin_ && (modelName_ == "EBHysteresisModel")) { // GENERAL NONLINEAR/HYSTERESIS CASE
              // get resistivity in coil regions via conductivity
              PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
              PtrCoefFct sigma = NULL;
              PtrCoefFct rho = NULL;

              constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
              sigma = actSDMat->GetScalCoefFnc(MAG_CONDUCTIVITY_SCALAR, Global::REAL);
              rho = CoefFunction::Generate( mp_,  Global::REAL, CoefXprBinOp(mp_, constOne, sigma, CoefXpr::OP_DIV ) );

              CoefXprVecScalOp jVec = CoefXprVecScalOp(mp_, iFct, boost::lexical_cast<std::string>(actPart.wireCrossSect), CoefXpr::OP_DIV);
              PtrCoefFct jVec_coefFct = CoefFunction::Generate(mp_, part, jVec);
              CoefXprVecScalOp rho_times_jVec = CoefXprVecScalOp(mp_, jVec_coefFct, rho, CoefXpr::OP_MULT);
              jFct[0] = CoefFunction::Generate(mp_, part, rho_times_jVec);
              coilCurrentDens_[actRegion] = jFct[0];
              curInt = GetCurrentDensityInt( 1.0, jFct[0] );
              curInt->SetName("(rho j,curlN): CoilIntegrator");
              LinearFormContext * coilContext = new LinearFormContext( curInt );
              coilContext->SetEntities( actSDList );
              coilContext->SetFeFunction( feFunc );
              coilContext->SetTypeLinearForm(2); // declared as nonlinear
              assemble_->AddLinearForm( coilContext );
            } else { // TOTALLY LINEAR TRANSIENT PROBLEM (all regions are linear)
              // get resistivity in coil regions via conductivity
              PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
              PtrCoefFct sigma = NULL;
              PtrCoefFct rho = NULL;

              constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
              sigma = actSDMat->GetScalCoefFnc(MAG_CONDUCTIVITY_SCALAR, Global::REAL);
              rho = CoefFunction::Generate( mp_,  Global::REAL, CoefXprBinOp(mp_, constOne, sigma, CoefXpr::OP_DIV ) );

              CoefXprVecScalOp jVec = CoefXprVecScalOp(mp_, iFct, boost::lexical_cast<std::string>(actPart.wireCrossSect), CoefXpr::OP_DIV);
              PtrCoefFct jVec_coefFct = CoefFunction::Generate(mp_, part, jVec);
              CoefXprVecScalOp rho_times_jVec = CoefXprVecScalOp(mp_, jVec_coefFct, rho, CoefXpr::OP_MULT);
              jFct[0] = CoefFunction::Generate(mp_, part, rho_times_jVec);
              coilCurrentDens_[actRegion] = jFct[0];
              curInt = GetCurrentDensityInt( 1.0, jFct[0] );
              curInt->SetName("(rho j,curlN): CoilIntegrator (all regions are linear)");
              LinearFormContext * coilContext = new LinearFormContext( curInt );
              coilContext->SetEntities( actSDList );
              coilContext->SetFeFunction( feFunc );
              assemble_->AddLinearForm( coilContext );
            }
          }
          else {
            EXCEPTION("ERROR: is_pseudo_time_stepping is either 1 or 0")
          }
        } // loop: parts
      }
    } // loop: coils
  }


  // **********************************************************
  // Definition of all bilinearforms
  // linear:    (p curlN, curlN) + (mu N, N)
  // nonlinear: (p curlN, curlN) + (db/dh N, N)
  // **********************************************************
  void MagEdgeHPDE::DefineStandardIntegrators(){

    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;

    // global storage for permeability makes the postprocessing easier
    perm_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, false));
    resistivity_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, false));
    
    // get the beta from the .xml file and set the default value
    Double beta = 1;
    myParam_->GetValue("penaltyFunctionParameter", beta, ParamNode::PASS);

    // Get the flag for pseudo time stepping
    // 1 ... pseudo time stepping ( simulate more time steps, but no time derivatives are involved )
    // 0 ... real time stepping (time derivatives are involved) 
    UInt is_pseudo_time_stepping = 1; // per default pseudo time stepping is used
    myParam_->GetValue("isPseudoTimeStepping", is_pseudo_time_stepping, ParamNode::PASS);

    // Get dt
    Double dt = 1e-3;
    //dt = myParam_->Get("analysis")->Get("transient")->Get("deltaT")->MathParse<double>();
    //dt = myParam_->Get("delta_t")->MathParse<Double>();
    myParam_->GetValue("delta_t", dt, ParamNode::PASS);

    // get FEFunction and space
    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_FIELD_INTENSITY];
    shared_ptr<FeSpace> feSpace = feFunc->GetFeSpace();

    PtrCoefFct magFieldCoef = this->GetCoefFct(MAG_FIELD_INTENSITY);

    // currently we are only allowed to have one hysteresis region
    bool moreThan1HystRegion = false;
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
        matModelCoefm_[actRegion]->Init(magFieldCoef, modelName_, dim_); 
      }

      // ===================================================================
      // PSEUDO TIME-STEPPING [START]
      // ===================================================================
      if (is_pseudo_time_stepping == 1){
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
        PtrCoefFct muNL = NULL;
        PtrCoefFct mu = NULL;
        BaseBDBInt *massInt = NULL;
        if (nonLinTypes.Find(PERMEABILITY) != -1){ // NONLINEAR/HYSTERESIS CASE
          if (modelName_ == "JilesAthertonModel"){ // just because, nobody wonders that it does not work with this model
            EXCEPTION("Jiles-Atherton model not implemented for MagEdgeHPDE");
          }
          else if (modelName_ == "EBHysteresisModel"){ // (this is the model we want to use!)
            // get nonlinear/hysteretic material parameter
            std::map<std::string, double> ParameterMap;
            std::map<std::string, string> StringParameterMap; 
            ParameterMap.clear();
            StringParameterMap.clear();
            actSDMat->GetString(StringParameterMap["weights_file_path"], MAG_WEIGHTS_FILE_PATH_EB);
            /// here is sill set the strong param map for the EBhst
            if(actSDMat->GetAnhystMagModel() == "analytic_anhysteresis"){
              actSDMat->GetString(StringParameterMap["anhyst_type"], MAG_ANHYST_TYPE_EB);
              if(actSDMat->GetAnhystFormula() == "atan"){
                actSDMat->GetString(StringParameterMap["anhyst_formula"], MAG_ANHYST_FORMULA_EB);
                actSDMat->GetScalar(ParameterMap["Ps"], MAG_PS_EB, Global::REAL);
                actSDMat->GetScalar(ParameterMap["A"], MAG_A_EB, Global::REAL);
              }

              if(actSDMat->GetAnhystFormula() == "pacejka"){
                actSDMat->GetString(StringParameterMap["anhyst_formula"], MAG_ANHYST_FORMULA_EB);
                actSDMat->GetScalar(ParameterMap["m_sat"], MAG_MSAT_PACEJKA_EB, Global::REAL);
                actSDMat->GetScalar(ParameterMap["a"], MAG_A_PACEJKA_EB, Global::REAL);
                actSDMat->GetScalar(ParameterMap["b"], MAG_B_PACEJKA_EB, Global::REAL);
                actSDMat->GetScalar(ParameterMap["c"], MAG_C_PACEJKA_EB, Global::REAL);

              }
            }else if(actSDMat->GetAnhystMagModel() == "multiscale_anhysteresis"){
              actSDMat->GetString(StringParameterMap["anhyst_type"], MAG_ANHYST_TYPE_EB);
              actSDMat->GetScalar(ParameterMap["AS"], MAG_MSM_AS, Global::REAL);
              actSDMat->GetScalar(ParameterMap["K1"], MAG_MSM_K1, Global::REAL);
              actSDMat->GetScalar(ParameterMap["K2"], MAG_MSM_K2, Global::REAL);
              actSDMat->GetScalar(ParameterMap["lambda100"], MAG_MSM_LAMBDA100, Global::REAL);
              actSDMat->GetScalar(ParameterMap["lambda111"], MAG_MSM_LAMBDA111, Global::REAL);
              actSDMat->GetScalar(ParameterMap["Ps"], MAG_MSM_PS, Global::REAL);
            }
            actSDMat->GetScalar(ParameterMap["jacobian_method"], MAG_JACOBIAN_METHOD_EB, Global::REAL);
            actSDMat->GetScalar(ParameterMap["approx_type"], MAG_APPROX_TYPE, Global::REAL);
            ParameterMap["isMH"] = 0;
            actSDMat->GetString(StringParameterMap["pinning_forces_weights_file"], MAG_PINNING_FORCES_WEIGHTS_EB);

            
            matModelCoefm_[actRegion]->InitModel(ParameterMap,StringParameterMap, actSDList);

            if(actSDMat->GetAnhystMagModel() == "multiscale_anhysteresis"){
              // check if we have a dependency of the anhysteretic curve with mechanical stress
              PtrCoefFct stresscoef;
              //get coeff-Fnc to evaluate the temperature
              StdVector<std::string> dispDofNames = feFunc->GetResultInfo()->dofNames;
              shared_ptr<EntityList> ent = ptGrid_->GetEntityList( EntityList::ELEM_LIST, regionName );
              ReadMaterialDependency( "permeability", dispDofNames, ResultInfo::VECTOR, false,
                                      ent, stresscoef, updatedGeo_ );
              matModelCoefm_[actRegion]->RegisterStressDependence(stresscoef);
            }
            // evaluate dbdh = mu
            muNL = matModelCoefm_[actRegion];
            nlFluxCoefm_[actRegion].reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_, true)); 
            nlFluxCoefm_[actRegion]->AddRegion(actRegion, matModelCoefm_[actRegion]);
            nlScalCoefm_[actRegion].reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_, true)); 
            nlScalCoefm_[actRegion]->AddRegion(actRegion, matModelCoefm_[actRegion]);
            if (dim_ == 2) {
              massInt = new BDBInt<>(new IdentityOperator<FeHCurl,2,1,Double>(),muNL,1.0,updatedGeo_);
            } else {
              massInt = new BDBInt<>(new IdentityOperator<FeHCurl,3,1,Double>(),muNL,1.0,updatedGeo_);
            }
            
            massInt->SetName("(db/dh N,N): IdentityIntegrator");
          }
        }
        else{ // LINEAR CASE (difference is described above)
          mu = actSDMat->GetScalCoefFnc(MAG_PERMEABILITY_SCALAR, Global::REAL);
          perm_->AddRegion(actRegion, mu);
          if (dim_ == 2) {
            massInt = new BBInt<>(new IdentityOperator<FeHCurl,2,1,Double>(),mu,1.0,updatedGeo_);
          } else {
            massInt = new BBInt<>(new IdentityOperator<FeHCurl,3,1,Double>(),mu,1.0,updatedGeo_);
          }
          
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
        // STIFFNESS INTEGRATOR [START] (p curlN, curlN)
        // The coefficent p depends on the material case:
        // ~ linear:      p = (1/(std::pow(beta,0.01)))*1*mu_regularize*std::pow(10,beta/2)
        // ~ nonlinear:   p = (1/beta^0.01)*||mu_regularize||_inf*10^(beta/2)
        //                since mu in the nonlinear case is a tensor the inf-norm is a good choice for mu_regularize
        // ~ hysteresis: same as in the nonlinear case
        // ====================================================================
        // get material coefficient (in this case: artificial resistivity that regularizes the PDE and is defined as 1/sigma_artificial = mu*10^(n/2)
        // where n is an integer number)
        Double mu_regularize;
        materials_[actRegion]->GetScalar( mu_regularize, MAG_PERMEABILITY_SCALAR, Global::REAL );
        PtrCoefFct p;
        PtrCoefFct p_nl;
        PtrCoefFct max_svd; // max. svd value of Jacobian matrix db/dh
        BaseBDBInt* curlcurl = NULL;

        if (nonLinTypes.Find(PERMEABILITY) != -1){ // NONLINEAR/HYSTERESIS CASE
          if (modelName_ == "JilesAthertonModel"){ // just because, nobody wonders that it does not work with this model
            EXCEPTION("Jiles-Atherton model not implemented for MagEdgeHPDE");
          }
          else if (modelName_ == "EBHysteresisModel"){ // (this is the model we want to use!)
            moreThan1HystRegion = true;
            if (dim_ == 2) {
              p = CoefFunction::Generate(mp_, Global::REAL,lexical_cast<std::string>(mu_regularize*std::pow(10,beta/2)));
              curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,2, Double>(), nlScalCoefm_[actRegion],1.0, updatedGeo_);
            } else {
              curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), nlScalCoefm_[actRegion],1.0, updatedGeo_);
            }
            curlcurl->SetName("(p_nl curlN,CurlN): CurlCurlIntegrator");
          }
        } else{ // LINEAR CASE
          //p = CoefFunction::Generate(mp_, Global::REAL,lexical_cast<std::string>((1/(std::pow(beta,0.01)))*1*mu_regularize*std::pow(10,beta/2)));
          p = CoefFunction::Generate(mp_, Global::REAL,lexical_cast<std::string>(mu_regularize*std::pow(10,beta/2)));
          if (dim_ == 2) {
            curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,2, Double>(), p,1.0, updatedGeo_);
          } else {
            curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), p,1.0, updatedGeo_);
          }
          curlcurl->SetName("(p curlN,CurlN): CurlCurlIntegrator");
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
        PtrCoefFct constZero = CoefFunction::Generate( mp_, Global::REAL, "0.0");
        if (dim_ == 2) {
          AUX_curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,2, Double>(), constZero,0.0, updatedGeo_);
        } else {
          AUX_curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), constZero,0.0, updatedGeo_);
        }
        AUX_curlcurl->SetName("(0.0*curlN,CurlN): Auxiliary,CurlCurlIntegrator (just for postprocessing: get curlh from h)");
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
      }
      // ===================================================================
      // PSEUDO TIME-STEPPING [END]
      // ===================================================================

      // ===================================================================
      // REAL TIME-STEPPING [START]
      // ===================================================================
      if (is_pseudo_time_stepping == 0){
        // ====================================================================
        // MASS INTEGRATOR [START] d/dt(c N, N)
        // The coefficient c of this integrator depends on the material case:
        // ~ linear:     use just the value of the permeability of this region.
        // ~ nonlinear:  use the derivative scaled by the current delta_t
        //               (1/delta_t)*(db/dh) for the coefficient c.
        //               AND the matrix is in openCFS a STIFFNESS MATRIX then.
        // ====================================================================
        double dt = 1e-3;       
        myParam_->GetValue("delta_t", dt, ParamNode::PASS); 
        PtrCoefFct muNL = NULL;
        PtrCoefFct mu = NULL;
        BaseBDBInt *massInt = NULL;
        if(nonLin_ && (modelName_ == "EBHysteresisModel")) { // GENERAL NONLINEAR/HYSTERESIS CASE
          if (nonLinTypes.Find(PERMEABILITY) != -1){ // NONLINEAR/HYSTERESIS SUBREGION
            if (modelName_ == "EBHysteresisModel"){ // EB HYSTERESIS MODEL
              std::map<std::string, double> ParameterMap;
              std::map<std::string, string> StringParameterMap; 
              ParameterMap.clear();
              StringParameterMap.clear();
              actSDMat->GetString(StringParameterMap["weights_file_path"], MAG_WEIGHTS_FILE_PATH_EB);
              /// here is sill set the strong param map for the EBhst
              if(actSDMat->GetAnhystMagModel() == "analytic_anhysteresis"){
                actSDMat->GetString(StringParameterMap["anhyst_type"], MAG_ANHYST_TYPE_EB);
                if(actSDMat->GetAnhystFormula() == "atan"){
                  actSDMat->GetString(StringParameterMap["anhyst_formula"], MAG_ANHYST_FORMULA_EB);
                  actSDMat->GetScalar(ParameterMap["Ps"], MAG_PS_EB, Global::REAL);
                  actSDMat->GetScalar(ParameterMap["A"], MAG_A_EB, Global::REAL);
                }

                if(actSDMat->GetAnhystFormula() == "pacejka"){
                  actSDMat->GetString(StringParameterMap["anhyst_formula"], MAG_ANHYST_FORMULA_EB);
                  actSDMat->GetScalar(ParameterMap["m_sat"], MAG_MSAT_PACEJKA_EB, Global::REAL);
                  actSDMat->GetScalar(ParameterMap["a"], MAG_A_PACEJKA_EB, Global::REAL);
                  actSDMat->GetScalar(ParameterMap["b"], MAG_B_PACEJKA_EB, Global::REAL);
                  actSDMat->GetScalar(ParameterMap["c"], MAG_C_PACEJKA_EB, Global::REAL);

                }
              }else if(actSDMat->GetAnhystMagModel() == "multiscale_anhysteresis"){
                actSDMat->GetString(StringParameterMap["anhyst_type"], MAG_ANHYST_TYPE_EB);
                actSDMat->GetScalar(ParameterMap["AS"], MAG_MSM_AS, Global::REAL);
                actSDMat->GetScalar(ParameterMap["K1"], MAG_MSM_K1, Global::REAL);
                actSDMat->GetScalar(ParameterMap["K2"], MAG_MSM_K2, Global::REAL);
                actSDMat->GetScalar(ParameterMap["lambda100"], MAG_MSM_LAMBDA100, Global::REAL);
                actSDMat->GetScalar(ParameterMap["lambda111"], MAG_MSM_LAMBDA111, Global::REAL);
                actSDMat->GetScalar(ParameterMap["Ps"], MAG_MSM_PS, Global::REAL);
              }
              actSDMat->GetScalar(ParameterMap["jacobian_method"], MAG_JACOBIAN_METHOD_EB, Global::REAL);
              actSDMat->GetScalar(ParameterMap["approx_type"], MAG_APPROX_TYPE, Global::REAL);
              ParameterMap["isMH"] = 0;
              actSDMat->GetString(StringParameterMap["pinning_forces_weights_file"], MAG_PINNING_FORCES_WEIGHTS_EB);

              matModelCoefm_[actRegion]->InitModel(ParameterMap,StringParameterMap, actSDList);

              if(actSDMat->GetAnhystMagModel() == "multiscale_anhysteresis"){
                // check if we have a dependency of the anhysteretic curve with mechanical stress
                PtrCoefFct stresscoef;
                //get coeff-Fnc to evaluate the temperature
                StdVector<std::string> dispDofNames = feFunc->GetResultInfo()->dofNames;
                shared_ptr<EntityList> ent = ptGrid_->GetEntityList( EntityList::ELEM_LIST, regionName );
                ReadMaterialDependency( "permeability", dispDofNames, ResultInfo::VECTOR, false,
                                        ent, stresscoef, updatedGeo_ );
                matModelCoefm_[actRegion]->RegisterStressDependence(stresscoef);
              }
              // evaluate dbdh = mu
              muNL = matModelCoefm_[actRegion];
              nlFluxCoefm_[actRegion].reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_, true)); 
              nlFluxCoefm_[actRegion]->AddRegion(actRegion, matModelCoefm_[actRegion]);
              if (dim_ == 2) {
                massInt = new BDBInt<>(new IdentityOperator<FeHCurl,2,1,Double>(),muNL,(1.0/dt),updatedGeo_);
              } else {
                massInt = new BDBInt<>(new IdentityOperator<FeHCurl,3,1,Double>(),muNL,(1.0/dt),updatedGeo_);
              }
              massInt->SetName("((1/dt)*(db/dh) N,N): IdentityIntegrator, nonlinear time stepping");
              BiLinFormContext * massContext = new BiLinFormContext(massInt, STIFFNESS );
              // save defined integrator
              massContext->SetEntities( actSDList, actSDList );
              massContext->SetFeFunctions( feFunc, feFunc );
              assemble_->AddBiLinearForm( massContext );
              massInts_[actRegion] = massInt; // insert mass integrator to list of defined mass integrators
            } else {
              EXCEPTION("NO OTHER HYSTERESIS MODEL IMPLEMENTED!")
            }
          } else { // LINEAR CASE (SUBREGION)
            // get linear material coefficient and save it in a global variable for postprocessing
            mu = actSDMat->GetScalCoefFnc(MAG_PERMEABILITY_SCALAR, Global::REAL);
            perm_->AddRegion(actRegion, mu);
            // define linear integrator
            if (dim_ == 2) {
              massInt = new BBInt<>(new IdentityOperator<FeHCurl,2,1,Double>(),mu,(1.0/dt),updatedGeo_);
            } else {
              massInt = new BBInt<>(new IdentityOperator<FeHCurl,3,1,Double>(),mu,(1.0/dt),updatedGeo_);
            }
            massInt->SetName("d/dt(mu N,N): IdentityIntegrator, linear time stepping");
            BiLinFormContext * massContext = new BiLinFormContext(massInt, STIFFNESS );
            // save defined integrator
            massContext->SetEntities( actSDList, actSDList );
            massContext->SetFeFunctions( feFunc, feFunc );
            assemble_->AddBiLinearForm( massContext );
            massInts_[actRegion] = massInt; // insert mass integrator to list of defined mass integrators
          }
        } else { // GENERAL LINEAR CASE (All regions are linear)
            // get linear material coefficient and save it in a global variable for postprocessing
            mu = actSDMat->GetScalCoefFnc(MAG_PERMEABILITY_SCALAR, Global::REAL);
            perm_->AddRegion(actRegion, mu);
            // define linear integrator
            if (dim_ == 2) {
              massInt = new BBInt<>(new IdentityOperator<FeHCurl,2,1,Double>(),mu,1.0,updatedGeo_);
            } else {
              massInt = new BBInt<>(new IdentityOperator<FeHCurl,3,1,Double>(),mu,1.0,updatedGeo_);
            }
            massInt->SetName("d/dt(mu N,N): IdentityIntegrator, all regions, linear time stepping");
            BiLinFormContext * massContext = new BiLinFormContext(massInt, DAMPING );
            // save defined integrator
            massContext->SetEntities( actSDList, actSDList );
            massContext->SetFeFunctions( feFunc, feFunc );
            assemble_->AddBiLinearForm( massContext );
            massInts_[actRegion] = massInt; // insert mass integrator to list of defined mass integrators        
        }
        // ====================================================================
        // MASS INTEGRATOR [END] d/dt(c N, N)
        // ====================================================================

        // ====================================================================
        // STIFFNESS INTEGRATOR [START] (rho curlN, curlN)
        // The coefficent rho depends on the material case:
        // ~ linear: rho is the resistivity of the region
        // ~ nonlinear: rho is the resistivity of the region
        // ====================================================================
        PtrCoefFct sigma = NULL; // conductivity of the material
        PtrCoefFct rho = NULL; // resistivity of the material
        PtrCoefFct constOne = NULL; // coefFct with 1 in every element
        BaseBDBInt* curlcurl = NULL;

        // get resistivity in every region via conductivity
        constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
        sigma = actSDMat->GetScalCoefFnc(MAG_CONDUCTIVITY_SCALAR, Global::REAL);
        rho = CoefFunction::Generate( mp_,  Global::REAL, CoefXprBinOp(mp_, constOne, sigma, CoefXpr::OP_DIV ) );
        resistivity_->AddRegion(actRegion, rho);

        if (dim_ == 2) {
          curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,2, Double>(), rho,1.0, updatedGeo_);
        } else {
          curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), rho,1.0, updatedGeo_);
        }
        curlcurl->SetName("(rho curlN,CurlN): CurlCurlIntegrator");
        BiLinFormContext * stiffContext = new BiLinFormContext(curlcurl, STIFFNESS );
        stiffContext->SetEntities( actSDList, actSDList );
        stiffContext->SetFeFunctions( feFunc, feFunc );
        assemble_->AddBiLinearForm( stiffContext );
        // ====================================================================
        // STIFFNESS INTEGRATOR [END] (rho curlN, curlN)
        // ====================================================================

        // ====================================================================
        // AUXILIARY CURL-CURL INTEGRATOR [START]
        // just define, but never assembled into global matrices
        // ====================================================================
        BaseBDBInt* AUX_curlcurl = NULL;
        PtrCoefFct constZero = NULL; // coefFct with 0 in every element
        constZero = CoefFunction::Generate( mp_, Global::REAL, "1.0");
        if (dim_ == 2) {
          AUX_curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,2, Double>(), constZero,0.0, updatedGeo_);
        } else {
          AUX_curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), constZero,0.0, updatedGeo_);
        }
        AUX_curlcurl->SetName("(0*curlN,CurlN): Auxiliary,CurlCurlIntegrator (just for postprocessing: get curlh from h)");
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
  void MagEdgeHPDE::DefineRhsLoadIntegrators(){
    // ============= IMPORTANT (START) ================================
    // all linearforms defined here only occur in the nonlinear case.
    // In the linear case only (p j,curlN) is needed and that is
    // defined in DefineCoilIntegrators()! For this reason an "if" is
    // wrapped around all statements, such that this method is empty {}
    // in the linear case.
    //
    // ADDITION: The formulation now also supports permanent magnets,
    //           wich are of course also considered in the linear case
    //           It is a linear form on the RHS
    //                   (b_r,N)
    //           where b_r is a remanence flux density that describes 
    //           the flux density of a permanent magnet region.
    // ============= IMPORTANT (END) ==================================
    
    // get FEFunctions and space
    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_FIELD_INTENSITY];
    shared_ptr<FeSpace> feSpace_reduced = feFunc->GetFeSpace();

    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    StdVector<std::string> vecDofNames = feFunc->GetResultInfo()->dofNames;
    LinearForm * lin1_pts = NULL; // (p curlh,curlN), only in the nonlinear case necessary (PSEUDO TIME STEPPING (pts))
    LinearForm * lin2_pts = NULL; // (b(h),N), only in the nonlinear case necessary (PSEUDO TIME STEPPING (pts))
    LinearForm * lin1_rts = NULL; // (rho curlh,curlN), only in the nonlinear case necessary (REAL TIME STEPPING (rts))
    LinearForm * lin2_rts = NULL; // ((1/delta_t)*b(h)_n,N), only in the nonlinear case necessary (REAL TIME STEPPING (rts))
    LinearForm * lin3_rts = NULL; // ((1/delta_t)*b(h)_{n-1},N), only in the nonlinear case necessary (REAL TIME STEPPING (rts))
    LinearForm * lin_br = NULL; // (b_r,N), describes linear permanent magnets in the linear/nonlinear case
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;

    bool coefUpdateGeo = true;
    bool isHystereticMat = false;

    Double beta = 1;
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
        if(nonLin_ && (modelName_ == "EBHysteresisModel")) {
          PtrCoefFct fluxDensityNL = NULL;
          fluxDensityNL = matModelCoefm_[actRegion];
          // ===============================================================================================
          // lin1_pts: (p curlh,curlN) [START]
          // curlh:   is obtained by the mag. field intensity h from the last Newton iteration.
          // p: regularization parameter that depends on the current scalar permeability mu
          //    -    nonlinear subregion: uses nlScalCoefm_ to get the norm of dbdh
          //    -    linear subregion:    uses linear case to get p
          // ===============================================================================================
          // generate the coefFct that is the multiplication of the curlh and p (p*curlh)
          Double mu_regularize;
          materials_[actRegion]->GetScalar( mu_regularize, MAG_PERMEABILITY_SCALAR, Global::REAL );
          PtrCoefFct p;
          PtrCoefFct p_nl;
          PtrCoefFct p_times_curlh;
          if (nonLinTypes.Find(PERMEABILITY) != -1){ // NONLINEAR CASE, NONLINEAR SUBREGION
            CoefXprVecScalOp temp = CoefXprVecScalOp(mp_, GetCoefFct( MAG_FIELD_INTENSITY_CURL ), nlScalCoefm_[actRegion], CoefXpr::OP_MULT);
            PtrCoefFct p_times_curlh = CoefFunction::Generate(mp_, Global::REAL, temp);
            if (dim_ == 2) {
              lin1_pts = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 2,Double>(),-1.0, p_times_curlh, volRegions, coefUpdateGeo);
            } else {
              lin1_pts = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3,Double>(),-1.0, p_times_curlh, volRegions, coefUpdateGeo);
            }
            lin1_pts->SetName("(p_nl curlh,curlN): residual");
          } else{ // NONLINEAR CASE, LINEAR SUBREGION
            p = CoefFunction::Generate(mp_, Global::REAL,lexical_cast<std::string>(mu_regularize*std::pow(10,beta/2)));
            CoefXprVecScalOp temp = CoefXprVecScalOp(mp_, GetCoefFct( MAG_FIELD_INTENSITY_CURL ), p, CoefXpr::OP_MULT);
            PtrCoefFct p_times_curlh = CoefFunction::Generate(mp_, Global::REAL, temp);
            if (dim_ == 2) {
              lin1_pts = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 2,Double>(),-1.0, p_times_curlh, volRegions, coefUpdateGeo);
            } else {
              lin1_pts = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3,Double>(),-1.0, p_times_curlh, volRegions, coefUpdateGeo);
            }
            lin1_pts->SetName("(p curlh,curlN): residual");
          }
          LinearFormContext *ctx = new LinearFormContext(lin1_pts);
          ctx->SetEntities(actSDList);
          ctx->SetFeFunction(feFunc);
          assemble_->AddLinearForm(ctx);
          // ===============================================================================================
          // lin1_pts: (p curlh,curlN) [END]
          // ===============================================================================================

          // ===============================================================================================
          // lin2_pts: (b(h),N) [START]
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
              if (dim_ == 2) {
                lin2_pts = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,2>(),(-1.0), fluxDensityNL, coefUpdateGeo, false);
              } else {
                lin2_pts = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3>(),(-1.0), fluxDensityNL, coefUpdateGeo, false);
              }
            }
            lin2_pts->SetName("(b(h),N): residual");
            lin2_pts->SetSolDependent();
            LinearFormContext *ctx = new LinearFormContext( lin2_pts );
            ctx->SetEntities( actSDList );
            ctx->SetFeFunction(feFunc);
            assemble_->AddLinearForm(ctx);
          } else { // NONLINEAR CASE, BUT LINEAR SUBREGION
            if (dim_ == 2) {
              lin2_pts = new BUIntegrator<Double>(new IdentityOperator<FeHCurl, 2>(),(-1.0), GetCoefFct( MAG_FLUX_DENSITY ), volRegions, coefUpdateGeo);
            } else {
              lin2_pts = new BUIntegrator<Double>(new IdentityOperator<FeHCurl, 3>(),(-1.0), GetCoefFct( MAG_FLUX_DENSITY ), volRegions, coefUpdateGeo);
            }
            lin2_pts->SetName("(b(h)_linear,N): residual");
            LinearFormContext *ctx = new LinearFormContext(lin2_pts);
            ctx->SetEntities(actSDList);
            ctx->SetFeFunction(feFunc);
            assemble_->AddLinearForm(ctx);
          }
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
        if(nonLin_ && (modelName_ == "EBHysteresisModel")) {
          PtrCoefFct fluxDensityNL = NULL;
          fluxDensityNL = matModelCoefm_[actRegion];
          // ===============================================================================================
          // lin1_rts: (rho curlh,curlN) [START]
          // curlh:   is obtained by the mag. field intensity h from the last Newton iteration.
          // rho: resistivity of the region
          // E = rho curlh ... ELEC_FIELD_INTENSITY
          // ===============================================================================================
          // generate the coefFct that is the multiplication of the curlh and rho (rho*curlh)
          PtrCoefFct sigma = NULL; // conductivity of the material
          PtrCoefFct rho = NULL; // resistivity of the material
          PtrCoefFct constOne = NULL; // coefFct with 1 in every element

          // get resistivity in every region via conductivity
          constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
          sigma = actSDMat->GetScalCoefFnc(MAG_CONDUCTIVITY_SCALAR, Global::REAL);
          rho = CoefFunction::Generate( mp_,  Global::REAL, CoefXprBinOp(mp_, constOne, sigma, CoefXpr::OP_DIV ) );

          PtrCoefFct rho_times_curlh = CoefFunction::Generate(mp_, Global::REAL, CoefXprVecScalOp(mp_, GetCoefFct( MAG_FIELD_INTENSITY_CURL ), rho, CoefXpr::OP_MULT));
          if (dim_ == 2) {
            lin1_rts = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 2,Double>(),-1.0, GetCoefFct( ELEC_FIELD_INTENSITY ), volRegions, coefUpdateGeo);
          } else {
            lin1_rts = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3,Double>(),-1.0, GetCoefFct( ELEC_FIELD_INTENSITY ), volRegions, coefUpdateGeo);
          }
          lin1_rts->SetName("(rho curlh,curlN): residual");

          LinearFormContext *ctx = new LinearFormContext(lin1_rts);
          ctx->SetEntities(actSDList);
          ctx->SetFeFunction(feFunc);
          ctx->SetTypeLinearForm(2); // declared as nonlinear
          assemble_->AddLinearForm(ctx);
          // ===============================================================================================
          // lin1_pts: (p curlh,curlN) [END]
          // ===============================================================================================

          // ===============================================================================================
          // lin2_rts: ((1/delta_t)*b(h)_n,N),N) [START]
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
              if (dim_ == 2) {
                lin2_rts = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,2>(),(-1/dt), fluxDensityNL, coefUpdateGeo, false);
              } else {
                lin2_rts = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3>(),(-1/dt), fluxDensityNL, coefUpdateGeo, false);
              }
            }
            lin2_rts->SetName("((1/delta_t)*b(h)_nonlin_n,N): residual");
            lin2_rts->SetSolDependent();
            LinearFormContext *ctx = new LinearFormContext( lin2_rts );
            ctx->SetEntities( actSDList );
            ctx->SetFeFunction(feFunc);
            ctx->SetTypeLinearForm(2); // declared as nonlinear
            assemble_->AddLinearForm(ctx);
          } else { // NONLINEAR CASE, BUT LINEAR SUBREGION
            if (dim_ == 2) {
              lin2_rts = new BUIntegrator<Double>(new IdentityOperator<FeHCurl, 2>(),(-1/dt), GetCoefFct( MAG_FLUX_DENSITY ), volRegions, coefUpdateGeo);
            } else {
              lin2_rts = new BUIntegrator<Double>(new IdentityOperator<FeHCurl, 3>(),(-1/dt), GetCoefFct( MAG_FLUX_DENSITY ), volRegions, coefUpdateGeo);
            }
            lin2_rts->SetName("((1/delta_t)*b(h)_linear_n,N): residual");
            LinearFormContext *ctx = new LinearFormContext(lin2_rts);
            ctx->SetEntities(actSDList);
            ctx->SetFeFunction(feFunc);
            ctx->SetTypeLinearForm(2); // declared as nonlinear
            assemble_->AddLinearForm(ctx);
          }
          // ===============================================================================================
          // lin2_rts: (b(h),N) [END]
          // ===============================================================================================

          // ===============================================================================================
          // lin3_rts: ((1/delta_t)*b(h)_{n-1},N),N) [START]
          // the b(h) is obtained by the used material model where the input is the mag. field intensity h
          // from the last TIME step.
          // ===============================================================================================
          if (nonLinTypes.Find(PERMEABILITY) != -1 && modelName_ != "nonlinearCurve"){
            if (modelName_ == "JilesAthertonModel")
            {
              EXCEPTION("Jiles-Atherton model not implemented for MagneticScalarPotentialPDE");
            }
            else if (modelName_ == "EBHysteresisModel") // NONLINEAR CASE, AND NONLINEAR SUBREGION
            {
              isHystereticMat = true;
              if (dim_ == 2) {
                  lin3_rts = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,2>(),(1/dt), fluxDensityNL, coefUpdateGeo, false);
                } else {
                  lin3_rts = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3>(),(1/dt), fluxDensityNL, coefUpdateGeo, false);
                }
            }
            lin3_rts->SetName("((1/delta_t)*b(h)_nonlin_{n-1},N): residual");
            lin3_rts->SetSolDependent();
            LinearFormContext *ctx = new LinearFormContext( lin3_rts );
            ctx->SetEntities( actSDList );
            ctx->SetFeFunction(feFunc);
            ctx->SetTypeLinearForm(1); // declared as linear
            assemble_->AddLinearForm(ctx);
          } else { // NONLINEAR CASE, BUT LINEAR SUBREGION
            if (dim_ == 2) {
              lin3_rts = new BUIntegrator<Double>(new IdentityOperator<FeHCurl, 2>(),(1/dt), GetCoefFct( MAG_FLUX_DENSITY ), volRegions, coefUpdateGeo);
            } else {
              lin3_rts = new BUIntegrator<Double>(new IdentityOperator<FeHCurl, 3>(),(1/dt), GetCoefFct( MAG_FLUX_DENSITY ), volRegions, coefUpdateGeo);
            }
            lin3_rts->SetName("((1/delta_t)*b(h)_linear_{n-1},N): residual");
            LinearFormContext *ctx = new LinearFormContext(lin3_rts);
            ctx->SetEntities(actSDList);
            ctx->SetFeFunction(feFunc);
            ctx->SetTypeLinearForm(1); // declared as linear
            assemble_->AddLinearForm(ctx);
          }
          // ===============================================================================================
          // lin3_rts: (b(h),N) [END]
          // ===============================================================================================
        }
      } 
      // ===================================================================
      // REAL TIME-STEPPING [END]
      // =================================================================== 
    }// end loop over entities

    // ===============================================================================================
    // lin_br: (b_r,N) [START]
    // The b_r has to be defined in the input .xml file for the simulation of permanent magnets
    // ===============================================================================================
    ReadRhsExcitation( "fluxDensity", vecDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    PtrCoefFct br;
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
          ent[i]->GetType() == EntityList::SURF_ELEM_LIST ) {
        EXCEPTION("Prescribed magnetic flux density can only be specified in a volume!")
      }

      double remanence_scale = 0;
      myParam_->GetValue("remanenceScale", remanence_scale, ParamNode::PASS);
      PtrCoefFct remanence_scale_coef_fct = CoefFunction::Generate(mp_, Global::REAL,lexical_cast<std::string>(remanence_scale));
      if (remanence_scale != 0) {
        PtrCoefFct Br_norm = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_, coef[i], CoefXpr::OP_NORM));
        PtrCoefFct factor = CoefFunction::Generate(mp_, Global::REAL,CoefXprBinOp(mp_,Br_norm,remanence_scale_coef_fct,CoefXpr::OP_DIV));
        CoefXprVecScalOp temp = CoefXprVecScalOp(mp_, coef[i], factor, CoefXpr::OP_DIV);
        br = CoefFunction::Generate(mp_, Global::REAL, temp); 
      } else {
        br = coef[i];
      }
      Brmap_[ent[i]->GetRegion()] = br; // Here we store the flux density remanence field for every region to have it ready for postprocessing

      if (dim_ == 2) {
        lin_br = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,2,1,Double>(), -1.0, br, coefUpdateGeo);
      } else {
        lin_br = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3,1,Double>(), -1.0, br, coefUpdateGeo);
      }
      lin_br->SetName("fluxIntegrator: (b_r,N)");
      LinearFormContext *ctx = new LinearFormContext( lin_br );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(feFunc);
      /* assemble_->AddLinearForm(ctx); */
      feFunc->AddEntityList(ent[i]);
      bRHSRegions_[ent[i]->GetRegion()] = br;
    }
    // ===============================================================================================
    // lin_br: (b_r,N) [END]
    // ===============================================================================================
  }


  // **********************************************************
  // definition of the primary result 
  // magnetic field intensity h
  // **********************************************************
  void MagEdgeHPDE::DefinePrimaryResults() {

    StdVector<std::string> vecComponents;
    if (dim_ == 2) {
      vecComponents = "x", "y";
    } else {
      vecComponents = "x", "y", "z";
    }
    

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

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MAG_FIELD_INTENSITY] = "fluxNormal";
    idbcSolNameMap_[MAG_FIELD_INTENSITY] = "potential";

  }


  // **********************************************************
  // definition of postprocessed results 
  // - magnetic flux density:  b = mu h
  // - source current density: j = curlh
  // **********************************************************
  void MagEdgeHPDE::DefinePostProcResults() {

    StdVector<std::string> vecComponents;
    if (dim_ == 2) {
      vecComponents = "x", "y";
    } else {
      vecComponents = "x", "y", "z";
    }
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

    // =====================================================
    // ELEC_FIELD_INTENSITY (START)
    // =====================================================
    shared_ptr<ResultInfo> elecIntens(new ResultInfo);
    elecIntens->resultType = ELEC_FIELD_INTENSITY;
    elecIntens->SetVectorDOFs(dim_, isaxi_);
    elecIntens->dofNames = vecComponents;
    elecIntens->unit = "V/m";
    elecIntens->definedOn = ResultInfo::ELEMENT;
    elecIntens->entryType = ResultInfo::VECTOR;
    // The recipe about how to actually evaluate E, is defined in FinalizePostProcResults()
    shared_ptr<CoefFunctionMulti> eFunc(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
    DefineFieldResult(eFunc, elecIntens);
    // =====================================================
    // ELEC_FIELD_INTENSITY (END)
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
    shared_ptr<CoefFunctionMulti> eCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[ELEC_FIELD_INTENSITY]);
    shared_ptr<CoefFunctionMulti> iCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_TOTAL_CURRENT_DENSITY]);
    StdVector<RegionIdType>::iterator regIt = regions_.Begin();
    regIt = regions_.Begin();

    for (; regIt != regions_.End(); ++regIt){
      // =====================================================
      // MAG_FLUX_DENSTIY (START)
      // b = mu*h + b_r (linear case)
      // b = b(h) + b_r (nonlinear case, with hysteresis this is not true)
      // =====================================================
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[*regIt];  // Just to find out which linear/nonlinear type is defined in this region
      PtrCoefFct b;
      if( nonLinTypes.Find(PERMEABILITY) != -1 && modelName_ != "nonlinearCurve" ){ // hysteretic/nonlinear case
        if(Brmap_.find(*regIt) != Brmap_.end()){ // There is a remancence flux density in the region prescribed
          PtrCoefFct br = Brmap_[*regIt];
          PtrCoefFct b_temp = nlFluxCoefm_[*regIt];
          b = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, br, b_temp, CoefXpr::OP_ADD));
        } else { // There is NO remancence flux density in the region prescribed
          b = nlFluxCoefm_[*regIt];
        }
        bCoef->AddRegion(*regIt, b);
      }
      else{ // classical nonlinear case and linear case
        if(Brmap_.find(*regIt) != Brmap_.end()){ // There is a remancence flux density in the region prescribed
          PtrCoefFct br;
          br = Brmap_[*regIt];
          PtrCoefFct b_temp = CoefFunction::Generate( mp_, Global::REAL, CoefXprVecScalOp(mp_, GetCoefFct( MAG_FIELD_INTENSITY ), perm_, CoefXpr::OP_MULT));
          b = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, br, b_temp, CoefXpr::OP_ADD));
        } else { // There is NO remancence flux density in the region prescribed
          b = CoefFunction::Generate( mp_, Global::REAL, CoefXprVecScalOp(mp_, GetCoefFct( MAG_FIELD_INTENSITY ), perm_, CoefXpr::OP_MULT));
        }
        bCoef->AddRegion(*regIt, b);      
      }
      // =====================================================
      // MAG_FLUX_DENSTIY (END)
      // =====================================================

      // =====================================================
      // ELEC_FIELD_INTENSITY (START)
      // =====================================================
      PtrCoefFct e = CoefFunction::Generate( mp_, Global::REAL, CoefXprVecScalOp(mp_, GetCoefFct( MAG_FIELD_INTENSITY_CURL ), resistivity_, CoefXpr::OP_MULT));
      eCoef->AddRegion(*regIt, e);
      // =====================================================
      // ELEC_FIELD_INTENSITY (END)
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
