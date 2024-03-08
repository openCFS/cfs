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
          BaseMaterial * actMat = NULL;
          actMat    = materials_[actRegion];
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
          // Get the flag for pseudo time stepping
          // 1 ... pseudo time stepping ( simulate more time steps, but no time derivatives are involved )
          // 0 ... real time stepping (time derivatives are involved) 
          Double is_pseudo_time_stepping = 1; // per default pseudo time stepping is used
          myParam_->GetValue("isPseudoTimeStepping", is_pseudo_time_stepping, ParamNode::PASS);
          // ===================================================================
          // PSEUDO TIME-STEPPING [START]
          // ===================================================================
          if (is_pseudo_time_stepping == 1){   
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
          }

          if (is_pseudo_time_stepping == 0){   
            // get resistivity in coil regions via conductivity
            PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
            PtrCoefFct sigma = NULL;
            PtrCoefFct rho = NULL;

            constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
            sigma = actMat->GetScalCoefFnc(MAG_CONDUCTIVITY_SCALAR, Global::REAL);
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
            assemble_->AddLinearForm( coilContext );
          }
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
    resistivity_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, false));
    
    // get the penaltyParameter from the .xml file and set the default value
    Double penaltyParameter = 1e-6;
    myParam_->GetValue("penaltyFunctionParameter", penaltyParameter, ParamNode::PASS);
    std::cout << "penalty bilinear" << penaltyParameter << std::endl;
    std::cout << penaltyParameter << std::endl;

    // Get the flag for pseudo time stepping
    // 1 ... pseudo time stepping ( simulate more time steps, but no time derivatives are involved )
    // 0 ... real time stepping (time derivatives are involved) 
    Double is_pseudo_time_stepping = 1; // per default pseudo time stepping is used
    myParam_->GetValue("isPseudoTimeStepping", is_pseudo_time_stepping, ParamNode::PASS);

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
          // ====================================================================
          // LINEAR CASE (difference is described above)
          PtrCoefFct mu = NULL;
          BaseBDBInt *massInt = NULL;

          mu = actMat->GetScalCoefFnc(MAG_PERMEABILITY_SCALAR, Global::REAL);
          perm_->AddRegion(actRegion, mu);
          massInt = new BBInt<>(new IdentityOperator<FeHCurl,3,1,Double>(),mu,1.0,updatedGeo_);
          massInt->SetName("d/dt(mu N,N): IdentityIntegrator");
          BiLinFormContext * massContext = new BiLinFormContext(massInt, DAMPING );
          massContext->SetEntities( actSDList, actSDList );
          massContext->SetFeFunctions( feFunc, feFunc );
          assemble_->AddBiLinearForm( massContext );
          massInts_[actRegion] = massInt; // insert mass integrator to list of defined mass integrators
          // ====================================================================
          // MASS INTEGRATOR [END] d/dt(c N, N)
          // ====================================================================

          // ====================================================================
          // STIFFNESS INTEGRATOR [START] (rho curlN, curlN)
          // The coefficent rho depends on the material case:
          // ~ linear: rho is the resistivity of the region
          // ====================================================================
          PtrCoefFct sigma = NULL; // conductivity of the material
          PtrCoefFct rho = NULL; // resistivity of the material
          PtrCoefFct constOne = NULL; // coefFct with 1 in every element
          BaseBDBInt* curlcurl = NULL;

          // get resistivity in every region via conductivity
          constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
          sigma = actMat->GetScalCoefFnc(MAG_CONDUCTIVITY_SCALAR, Global::REAL);
          rho = CoefFunction::Generate( mp_,  Global::REAL, CoefXprBinOp(mp_, constOne, sigma, CoefXpr::OP_DIV ) );
          resistivity_->AddRegion(actRegion, rho);

          curlcurl = new BBInt<>(new  CurlOperator<FeHCurl,3, Double>(), rho,1.0, updatedGeo_);
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
          // ====================================================================
          BaseBDBInt* AUX_curlcurl = NULL;
          constOne = CoefFunction::Generate( mp_, Global::REAL, "1.0");
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
        }
        // ===================================================================
        // REAL TIME-STEPPING [END]
        // ===================================================================
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
    //
    // ADDITION: The formulation now also supports permanent magnets,
    //           wich are of course also considered in the linear case
    //           It is a linear form on the RHS
    //                   (b_r,N)
    //           where b_r is a remanence flux density that describes 
    //           the flux density of a perfect permanent magnet region.
    // ============= IMPORTANT (END) ==================================
    
    // get FEFunctions and space
    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_FIELD_INTENSITY];
    shared_ptr<FeSpace> feSpace_reduced = feFunc->GetFeSpace();

    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    StdVector<std::string> vecDofNames = feFunc->GetResultInfo()->dofNames;
    LinearForm * lin1 = NULL; // (rho_art curlh,curlN), only in the nonlinear case necessary
    LinearForm * lin2 = NULL; // (b(h),N), only in the nonlinear case necessary
    LinearForm * lin3 = NULL; // (b_r,N), describes permanent magnets
    RegionIdType actRegion;

    bool coefUpdateGeo = true;
    bool isHystereticMat = false;

    Double penaltyParameter = 1e-6;
    myParam_->GetValue("penaltyFunctionParameter", penaltyParameter, ParamNode::PASS);

    // iterate over the region (or materials)
    for (UInt iRegion = 0; iRegion < regions_.GetSize(); iRegion++) {
      // set current region and material
      actRegion = regions_[iRegion];
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];
      std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );

      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList(new ElemList(ptGrid_));
      actSDList->SetRegion(actRegion);

      if(nonLin_ && (modelName_ == "EBHysteresisModel")) {
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
      }
    }// end loop over entities

    // ===============================================================================================
    // lin3: (b_r,N) [START]
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
      lin3 = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3,1,Double>(), -1.0, br, coefUpdateGeo);

      lin3->SetName("fluxIntegrator: (b_r,N)");
      LinearFormContext *ctx = new LinearFormContext( lin3 );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(feFunc);
      assemble_->AddLinearForm(ctx);
      feFunc->AddEntityList(ent[i]);
      bRHSRegions_[ent[i]->GetRegion()] = br;
    }
    // ===============================================================================================
    // lin3: (b_r,N) [END]
    // ===============================================================================================
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
          PtrCoefFct b_temp = nlFluxCoef_;
          b = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, br, b_temp, CoefXpr::OP_ADD));
        } else { // There is NO remancence flux density in the region prescribed
          b = nlFluxCoef_;
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
