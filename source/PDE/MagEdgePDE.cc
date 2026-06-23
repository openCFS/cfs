#include <fstream>

#include "MagEdgePDE.hh"

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

// time stepping
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

#include "Driver/MultiHarmonicDriver.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
namespace CoupledField
{

  // declare class specific logging stream
  DEFINE_LOG(magEdgePde, "magEdgePde")

  // **************
  //  Constructor
  // **************
  MagEdgePDE::MagEdgePDE(Grid *aptgrid, PtrParamNode paramNode,
                         PtrParamNode infoNode,
                         shared_ptr<SimState> simState, Domain *domain)
      : MagBasePDE(aptgrid, paramNode, infoNode, simState, domain)
  {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_ = "magneticEdge";
    pdematerialclass_ = ELECTROMAGNETIC;
    formulation_ = MagBasePDE::EDGE;

    //! Always use updated Lagrangian formulation
    updatedGeo_ = true;

    // default is false
    useGradFields_ = paramNode->Get("useGradientFields")->As<bool>();
    onlyVacuum_ = paramNode->Get("onlyVacuum")->As<bool>();

    // check if we have a 3d setup
    if (!is3d_)
      EXCEPTION("MagEdgePDE is just implemented for 3D setups!");

    multiHarmCoef_.reset(new CoefFunctionHarmBalance<Complex>());

    // Be aware that the hysteresis via CoefFunctionMaterialModel does not use the isHysteresis_ flag!
    matModelCoefm_.clear(); //matModelCoef_.reset(new CoefFunctionMaterialModel<Complex>());

    // init the nonlinear_field_intensity_coef_
    //nonlinear_field_intensity_coef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_, true));
  }

  // *************
  //  Destructor
  // *************
  MagEdgePDE::~MagEdgePDE()
  {
  }

  // *****************************
  //  Definition of Integrators
  // *****************************
  void MagEdgePDE::DefineIntegrators()
  {
    this->DefineStandardIntegrators();
    // in MagBasePDE
    DefineCoilIntegrators(1.0);
  }

  void MagEdgePDE::DefineStandardIntegrators()
  {

    // FIXME This is a dirty fix, the problem is described in gitlab issue #5
    if (analysistype_ == MULTIHARMONIC)
    {
      reluc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, true));
    }

    RegionIdType actRegion;
    BaseMaterial *actMat = NULL;

    // initially, check for regularization factor
    Double regularizationFactor = 1e-6;
    myParam_->GetValue("penaltyFactor", regularizationFactor, ParamNode::PASS);

    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_POTENTIAL];
    shared_ptr<FeSpace> feSpace = feFunc->GetFeSpace();

    PtrCoefFct magFluxCoef = this->GetCoefFct(MAG_FLUX_DENSITY);
    // Create new harmonic balance coefficient function and register the regions and material
    UInt baseFreq = 0, N, M, nFFT;
    if (analysistype_ == MULTIHARMONIC)
    {
      baseFreq = dynamic_cast<MultiHarmonicDriver *>(domain_->GetSingleDriver())->baseFreq_;
      N = dynamic_cast<MultiHarmonicDriver *>(domain_->GetSingleDriver())->numHarmonics_N_;
      M = dynamic_cast<MultiHarmonicDriver *>(domain_->GetSingleDriver())->numHarmonics_M_;
      nFFT = dynamic_cast<MultiHarmonicDriver *>(domain_->GetSingleDriver())->numFFT_;
      multiHarmCoef_->Init(feFunc, feSpace, regions_, materials_, ptGrid_, magFluxCoef, N, M, baseFreq, nFFT, modelName_, magFluxCoef);
    }

    for (UInt iRegion = 0; iRegion < regions_.GetSize(); iRegion++)
    {
      actRegion = regions_[iRegion];
      actMat = materials_[actRegion];
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region", "name", regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      feSpace->SetRegionApproximation(actRegion, polyId, integId);

      if (useGradFields_)
      {
        feSpace->SetUseGradients(actRegion);
      }

          // Init material model for hysteretic transient analysis
      if (((analysistype_ == STATIC) || (analysistype_ == TRANSIENT)) && nonLin_ && (modelName_ != "nonlinearCurve"))
      {
        //matModelCoef_->Init(magFluxCoef, modelName_, dim_);
        matModelCoefm_[actRegion].reset(new CoefFunctionMaterialModel<Complex>());
        matModelCoefm_[actRegion]->Init(magFluxCoef, modelName_, dim_); 
      }

      // get possible nonlinearities defined in this region
      StdVector<NonLinType> matDepenTypes = regionMatDepTypes_[actRegion]; // material dependency
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];   // non-linearity

      // create new entity list
      shared_ptr<ElemList> actSDList(new ElemList(ptGrid_));
      actSDList->SetRegion(actRegion);

      // pass entitylist ot fespace / fefunction
      feFunc->AddEntityList(actSDList);

      std::cout << "nonLinTypes.Find(PERMEABILITY): " << nonLinTypes.Find(PERMEABILITY) << "\n";
      if (nonLinTypes.Find(PERMEABILITY) != -1) {
        if (modelName_ == "invEBHysteresisModel"){ // use inverse energy-based Hysteresis model
          // define needed variables
          PtrCoefFct nu_nonlinear_eb = NULL;
          PtrCoefFct field_intensity_eb = NULL;
          BaseBDBInt *stiffInt = NULL;

          // init. EB Material Model
          std::map<std::string, double> ParameterMap;
          std::map<std::string, std::string> StringParameterMap;
          if (actMat->GetAnhystMagModel() == "analytic_anhysteresis"){
            actMat->GetString(StringParameterMap["anhyst_type"], MAG_ANHYST_TYPE_INVEB);
            if (actMat->GetAnhystFormula() == "tan"){
              actMat->GetString(StringParameterMap["anhyst_formula"], MAG_ANHYST_FORMULA_INVEB);
              actMat->GetScalar(ParameterMap["Js"], MAG_JS_INVEB, Global::REAL);
              actMat->GetScalar(ParameterMap["A"], MAG_A_INVEB, Global::REAL);
            } else if (actMat->GetAnhystFormula() == "brauer") {
              actMat->GetString(StringParameterMap["anhyst_formula"], MAG_ANHYST_FORMULA_INVEB);
              actMat->GetScalar(ParameterMap["p_0"], MAG_P0_INVEB, Global::REAL);
              actMat->GetScalar(ParameterMap["p_1"], MAG_P1_INVEB, Global::REAL);
              actMat->GetScalar(ParameterMap["p_2"], MAG_P2_INVEB, Global::REAL);
            } else if (actMat->GetAnhystFormula() == "lookuptable") {
              actMat->GetString(StringParameterMap["anhyst_formula"], MAG_ANHYST_FORMULA_INVEB);
              actMat->GetString(StringParameterMap["lookup_table_file"], MAG_LOOKUP_TABLE_FILE_INVEB);
            }
          }
          actMat->GetString(StringParameterMap["weights_file_path"], MAG_WEIGHTS_FILE_PATH_EB);
          actMat->GetString(StringParameterMap["pinning_forces_weights_file"], MAG_PINNING_FORCES_WEIGHTS_INVEB);
          actMat->GetScalar(ParameterMap["approx_type"], MAG_APPROX_TYPE, Global::REAL);
          ParameterMap["isMH"] = 0;
          actMat->GetScalar(ParameterMap["jacobian_method"], MAG_JACOBIAN_METHOD_INVEB, Global::REAL);
          matModelCoefm_[actRegion]->InitModel(ParameterMap, StringParameterMap, actSDList);
          nu_nonlinear_eb = matModelCoefm_[actRegion]; //matModelCoef_;

          nlFieldIntensitym_[actRegion].reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_, true)); 
          nlFieldIntensitym_[actRegion]->AddRegion(actRegion, matModelCoefm_[actRegion]);
          //nonlinear_field_intensity_coef_->AddRegion(actRegion, matModelCoefm_[actRegion]);

          PtrCoefFct curCoef = NULL;
          curCoef = materials_[actRegion]->GetTensorCoefFnc(MAG_RELUCTIVITY_TENSOR, FULL, Global::REAL);

          // define integrators
          if (dim_ == 3){ // 3D case
            stiffInt = new BDBInt<>(new CurlOperator<FeHCurl, 3, Double>(), nu_nonlinear_eb, 1.0, updatedGeo_);
          } else {
            EXCEPTION("For Edge elements, only 3D is possible!")
          }
          stiffInt->SetName("(dH/dB curl A, curl A',nonlinear case, nonlinear subregion)");
          BiLinFormContext *stiffContext = new BiLinFormContext(stiffInt, STIFFNESS);
          stiffContext->SetEntities(actSDList, actSDList);
          stiffContext->SetFeFunctions(feFunc, feFunc);
          assemble_->AddBiLinearForm(stiffContext);
          bdbInts_.insert(std::pair<RegionIdType, BaseBDBInt *>(actRegion, stiffInt));
        } else if ((nonLinTypes.Find(HYSTERESIS) == -1) && ((nonLinTypes.GetSize() > 0) || (analysistype_ == MULTIHARMONIC))){
          //        std::cout <<  "StiffnessMatrix - Case1: NonLin or Multiharmonic" << std::endl;

          // ================================
          //  Nonlinear Stiffness Integrator
          // =================================

          // create stiffness integrator
          // BaseBOperator* bOp = new CurlOperator<FeHCurl,3, Double>();
          CoefFunctionOpt *cfo = NULL; // we might do optimization and then we have such a thing
          PtrCoefFct magFluxCoef = this->GetCoefFct(MAG_FLUX_DENSITY);
          PtrCoefFct nuNl = NULL;
          if (analysistype_ == MULTIHARMONIC)
          {
            // register element list of region
            bool nL = (nonLinTypes.GetSize() > 0) ? true : false;
            nuNl = multiHarmCoef_->GenerateMatCoefFnc(iRegion, "Reluctivity", nL, actSDList);
          }
          else
          {
            if (nonLinTypes.Find(PERMEABILITY) != -1 && matDepenTypes.Find(PERMEABILITY) != -1)
            {
              // case of temperature-dependent BH-curves (solution-dependent BH curve AND temperature-dependency!)
              // ========================================
              //  Temperature dependent BH curves - Nonlinear Stiffness Integrator
              // ========================================
              if (analysistype_ != TRANSIENT)
              {
                EXCEPTION("Temperature-dependent BH-curves only implemented for transient magnetic analysis");
              }
              // CoefFct for temperature dependent permeability
              PtrCoefFct tempcoef;

              // get coeff-Fnc to evaluate the temperature
              StdVector<std::string> dispDofNames = feFunc->GetResultInfo()->dofNames;
              shared_ptr<EntityList> ent = ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName);
              ReadMaterialDependency("permeability", dispDofNames, ResultInfo::SCALAR, false,
                                     ent, tempcoef, updatedGeo_);

              // inform the ElectroMagneticMaterial about the flux density and also the temperature-dependency
              nuNl = actMat->GetScalCoefFncNonLin(MAG_RELUCTIVITY_SCALAR, Global::REAL, magFluxCoef, tempcoef);
            }
            else
            {
              // ========================================
              //  Classic Nonlinear Stiffness Integrator
              // ========================================
              nuNl = actMat->GetScalCoefFncNonLin(MAG_RELUCTIVITY_SCALAR, Global::REAL, magFluxCoef);
            }
          }

          // replace in optimization case
          if (domain->HasDesign())
          {
            cfo = new CoefFunctionOpt(domain->GetDesign(), nuNl, MAG_RELUCTIVITY_SCALAR, this);
            nuNl.reset(cfo);
          }

          // compute permeability
          PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");

          if (analysistype_ == MULTIHARMONIC)
          {
            PtrCoefFct permeability = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, constOne, nuNl, CoefXpr::OP_DIV));
            matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, permeability);
          }
          else
          {
            PtrCoefFct permeability = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, constOne, nuNl, CoefXpr::OP_DIV));
            matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, permeability);
          }

          BaseBDBInt *stiff1 = NULL;
          if (analysistype_ == MULTIHARMONIC)
          {
            stiff1 = new BBInt<Complex>(new CurlOperator<FeHCurl, 3, Double>(), nuNl, (Complex)1.0, updatedGeo_);
          }
          else
          {
            stiff1 = new BBInt<>(new CurlOperator<FeHCurl, 3, Double>(), nuNl, 1.0, updatedGeo_);
          }

          stiff1->SetName("CurlCurlIntegrator-NL");

          BiLinFormContext *stiffContext = new BiLinFormContext(stiff1, STIFFNESS);
          stiffContext->SetEntities(actSDList, actSDList);
          stiffContext->SetFeFunctions(feFunc, feFunc);
          assemble_->AddBiLinearForm(stiffContext);
          // Important: Add bdb-integrator to global list, as we need them later
          // for calculation of postprocessing results
          bdbInts_.insert(std::pair<RegionIdType, BaseBDBInt *>(actRegion, stiff1));
          // add also material to global, distributed reluctivity coefficient function
          reluc_->AddRegion(actRegion, nuNl);

          // ================================================
          //  Nonlinear Stiffness Integrator (only Newton )
          // ================================================
          // Note: currently we set the nonlinear method hard-coded to NEWTON for
          // testing purpose
          if (nonLinMethod_ == NEWTON)
          {

            // BaseBOperator* bOp = new CurlOperator<FeHCurl,3, Double>();
            PtrCoefFct magFluxCoef = this->GetCoefFct(MAG_FLUX_DENSITY);
            PtrCoefFct nuDeriv = actMat->GetTensorCoefFncNonLin(MAG_RELUCTIVITY_DERIV, FULL, Global::REAL, magFluxCoef);

            // create stiffness integrator
            BiLinearForm *stiff2 = NULL;
            stiff2 = new BDBInt<>(new CurlOperator<FeHCurl, 3, Double>(), nuDeriv, 1.0, updatedGeo_);
            stiff2->SetName("CurlCurlIntegrator-NL-Newton");
            stiff2->SetNewtonBiLinearForm();

            BiLinFormContext *stiffContext2 = new BiLinFormContext(stiff2, STIFFNESS);
            stiffContext2->SetEntities(actSDList, actSDList);
            stiffContext2->SetFeFunctions(feFunc, feFunc);
            assemble_->AddBiLinearForm(stiffContext2);
          }
        }
      } else {
        //        std::cout << "StiffnessMatrix - Case2: Hyst or Linear" << std::endl;
        // ***************************************
        // HYSTERESIS + LINEAR PART
        // ***************************************
        PtrCoefFct curCoef;
        BaseBDBInt *curlcurl = NULL;
        if (nonLinTypes.Find(HYSTERESIS) != -1)
        {
          /* for both the delta material method as well as the std fixpoint method we have to know
           * which regions are affected by hysteresis
           */
          curCoef = GenerateHystereticCoefFunctions(actRegion);
          curlcurl = GeHystStiffInt(1.0, curCoef);

          //        std::cout << "Add full tensor" << std::endl;
          relucTensor_->AddRegion(actRegion, curCoef); // full tensor
          //        std::cout << "Add scalar - hyst - actRegionId = " << actRegion << std::endl;
          Global::ComplexPart coefType = (analysistype_ == MULTIHARMONIC) ? Global::COMPLEX : Global::REAL;
          PtrCoefFct curCoefDet = CoefFunction::Generate(mp_, coefType, CoefXprUnaryOp(mp_, curCoef, CoefXpr::OP_DET));
          reluc_->AddRegion(actRegion, curCoefDet);
        }
        else
        {
          // ===============================
          //  Standard Stiffness Integrator (Linear case)
          // ===============================
          PtrCoefFct curCoef = NULL;
          CoefFunctionOpt *cfo = NULL; // we might do optimization and then we have such a thing

          if (analysistype_ == MULTIHARMONIC)
          {
            // For MULTIHARMONIC mode, use the harmonic balance coefficient function
            // This matches the behavior in the else if branch for MULTIHARMONIC
            bool nL = (nonLinTypes.GetSize() > 0) ? true : false;
            curCoef = multiHarmCoef_->GenerateMatCoefFnc(iRegion, "Reluctivity", nL, actSDList);

            // compute permeability (COMPLEX for MULTIHARMONIC)
            PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");
            PtrCoefFct permeability = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, constOne, curCoef, CoefXpr::OP_DIV));
            matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, permeability);

            // Create COMPLEX integrator for MULTIHARMONIC
            curlcurl = new BBInt<Complex>(new CurlOperator<FeHCurl, 3, Double>(), curCoef, (Complex)1.0, updatedGeo_);
          }
          else
          {
            // Standard case: use material coefficient directly
            PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");
            PtrCoefFct mu0 = CoefFunction::Generate(mp_, Global::REAL, "795774.716369");
            PtrCoefFct permeability = NULL;
            if (onlyVacuum_)
            {
              curCoef = mu0;
              permeability = mu0;
            }
            else
            {
              curCoef = actMat->GetScalCoefFnc(MAG_RELUCTIVITY_SCALAR, Global::REAL);
              permeability = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, constOne, curCoef, CoefXpr::OP_DIV));
            }
            matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, permeability);

            if (domain->HasDesign())
            {
              cfo = new CoefFunctionOpt(domain->GetDesign(), curCoef, MAG_RELUCTIVITY_SCALAR, this);
              curCoef.reset(cfo);
            }

            curlcurl = new BBInt<>(new CurlOperator<FeHCurl, 3, Double>(), curCoef, 1.0, updatedGeo_);

            // when we have a CoefFunctionOpt, we tell it the proper form, which we only have now
            if (cfo)
              cfo->SetForm(curlcurl);
          }
          //        std::cout << "Add scalar - linear - actRegionId = " << actRegion << std::endl;
          // add also material to global, distributed reluctivity coefficient function
          reluc_->AddRegion(actRegion, curCoef);
        }

        curlcurl->SetName("CurlCurlIntegrator");

        BiLinFormContext *stiffContext = new BiLinFormContext(curlcurl, STIFFNESS);
        stiffContext->SetEntities(actSDList, actSDList);
        stiffContext->SetFeFunctions(feFunc, feFunc);
        assemble_->AddBiLinearForm(stiffContext);

        // Important: Add bdb-integrator to global list, as we need them later
        // for calculation of postprocessing results
        bdbInts_.insert(std::pair<RegionIdType, BaseBDBInt *>(actRegion, curlcurl));

      } // END OF NONLIN/LIN PART

      // ============================================================
      // Mass Matrix
      // ============================================================
      if ((matDepenTypes.Find(NLELEC_CONDUCTIVITY) != -1 && nonLinTypes.GetSize() == 0) ||
          (matDepenTypes.Find(NLELEC_CONDUCTIVITY) != -1 && analysistype_ == MULTIHARMONIC)) 
      {
        // =================================================
        // Pure temperature-dependent nonlinear Mass Matrix
        // =================================================
        // purely temperature dependent (no further non linearity)
        PtrCoefFct coef;
        shared_ptr<BaseFeFunction> myFct = feFunctions_[MAG_POTENTIAL];
        StdVector<std::string> dispDofNames = myFct->GetResultInfo()->dofNames;
        shared_ptr<EntityList> ent = ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName);

        // get coeff-Fnc to evaluate the temperature
        ReadMaterialDependency("electricConductivity", dispDofNames, ResultInfo::SCALAR, false,
                               ent, coef, updatedGeo_);
        // coef-Fnc for electric conductivity
        PtrCoefFct condNL = actMat->GetScalCoefFncNonLin(MAG_CONDUCTIVITY_SCALAR, Global::REAL, coef);
        if (analysistype_ == STATIC)
        {
          EXCEPTION("No temperature-dependent conductivity implemented for static magnetic analysis");
        }
        conduc_->AddRegion(actRegion, condNL);
        BaseBDBInt *massInt;
        BiLinFormContext *massContext;
        massInt = new BBIntMassEdge<>(new IdentityOperator<FeHCurl, 3, 1, Double>(),
                                      condNL, 1.0, updatedGeo_);
        massContext = new BiLinFormContext(massInt, DAMPING);
        massInt->SetName("MassIntegrator-NL");
        massContext->SetEntities(actSDList, actSDList);
        massContext->SetFeFunctions(feFunc, feFunc);
        assemble_->AddBiLinearForm(massContext);

        // insert mass integrator to list of defined mass integrators
        massInts_[actRegion] = massInt;
      } else if ( nonLinTypes.Find(PERMEABILITY) != -1  && modelName_ == "invEBHysteresisModel"  ) {
        PtrCoefFct regularization_parameter;
        BaseBDBInt *massInt;
        BiLinFormContext *massContext;

        regularization_parameter = CoefFunction::Generate(mp_, Global::REAL, lexical_cast<std::string>(1e-2 * 1000*4*M_PI*1e-7));
        massInt = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<FeHCurl, 3, Double>(), regularization_parameter, 1.0);
        massInt->SetName("MassIntegrator (regularization_parameter)");
        massContext = new BiLinFormContext(massInt, STIFFNESS);
        massContext->SetEntities(actSDList, actSDList);
        massContext->SetFeFunctions(feFunc, feFunc);
        assemble_->AddBiLinearForm(massContext);
        // insert mass integrator to the list of defined mass integrators
        massInts_[actRegion] = massInt;
      } else {
        //        std::cout <<  "MassMatrix - Case2: Standard linear Mass Matrix" << std::endl;
        // =================================================
        // Standard linear Mass Matrix
        // =================================================
        PtrCoefFct conductivityCoeff;
        Double conductivity = 0.0;
        if (materials_[actRegion]->GetSymmetryType(MAG_CONDUCTIVITY_TENSOR) == BaseMaterial::GENERAL)
        {
          conductivityCoeff = materials_[actRegion]->GetTensorCoefFnc(MAG_CONDUCTIVITY_TENSOR, FULL, Global::REAL);
          Matrix<Double> conduc;
          materials_[actRegion]->GetTensor(conduc, MAG_CONDUCTIVITY_TENSOR, Global::REAL);
          conductivity = (conduc[0][0] + conduc[1][1] + conduc[2][2]) / 3;
        }
        else
        {
          materials_[actRegion]->GetScalar(conductivity, MAG_CONDUCTIVITY_SCALAR, Global::REAL);
          conductivityCoeff = CoefFunction::Generate(mp_, Global::REAL, lexical_cast<std::string>(conductivity));
        }
        // regularize
        bool scaleByEdgeSize = false;
        if (conductivity < 1e-10 || analysistype_ == STATIC)
        {
          Matrix<Double> reluc;
          // get tensor of permeability and determine max. value
          materials_[actRegion]->GetTensor(reluc, MAG_RELUCTIVITY_TENSOR, Global::REAL);
          conductivityCoeff = CoefFunction::Generate(mp_, Global::REAL, lexical_cast<std::string>(regularizationFactor * reluc[0][0]));
          scaleByEdgeSize = true;
          regularizedRegions_.insert(actRegion);
          WARN("Regularization activated for region '" << regionName << "' (conductivity=" << conductivity << ") - using artificial conductivity " << regularizationFactor * reluc[0][0]);
        }

        conduc_->AddRegion(actRegion, conductivityCoeff);
        BaseBDBInt *massInt;
        BiLinFormContext *massContext;
        if (analysistype_ == STATIC)
        {
          // we have to guarantee, that we add some mass to curl-curl integrator.
          // Additionally, the integrator gets scaled by the edge size for a uniform
          // conditioning
          massInt = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<FeHCurl, 3, Double>(),
                                        conductivityCoeff, 1.0);
          massInt->SetName("MassIntegrator");
          massContext = new BiLinFormContext(massInt, STIFFNESS);
        }
        else
        {
          // here we add the "normal" mass integrator, which gets not scaled by the
          // edge size

          if (scaleByEdgeSize)
          { // this is for non-conducting regions
            massInt = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<FeHCurl, 3, Double>(),
                                          conductivityCoeff, 1.0, updatedGeo_);
            massContext = new BiLinFormContext(massInt, STIFFNESS);
          }
          else
          {
            if (conductivityCoeff->GetDimType() == CoefFunction::TENSOR)
            {
              massInt = new BDBInt<>(new IdentityOperator<FeHCurl, 3, 1, Double>(), conductivityCoeff, 1.0, updatedGeo_);
            }
            else
            {
              massInt = new BBIntMassEdge<>(new IdentityOperator<FeHCurl, 3, 1, Double>(), conductivityCoeff, 1.0, updatedGeo_);
            }
            massContext = new BiLinFormContext(massInt, DAMPING);
          }
          massInt->SetName("MassIntegrator");
        }
        massContext->SetEntities(actSDList, actSDList);
        massContext->SetFeFunctions(feFunc, feFunc);
        assemble_->AddBiLinearForm(massContext);
        // insert mass integrator to list of defined mass integrators
        massInts_[actRegion] = massInt;
      } // End of nonlin/lin mass matrix part

      // ====================================================================
      // check for velocity
      // ====================================================================
      std::string velocityId = curRegNode->Get("velocityId")->As<std::string>();
      if (velocityId != "")
      {
        bool coefUpdateGeo;
        ReadRegionVelocityField(velocityId, actSDList, actRegion, coefUpdateGeo);

        // coef-Fnc for electric conductivity
        Matrix<Double> reluc;
        Double conductivity = 0.0;

        // get conductivity
        materials_[actRegion]->GetScalar(conductivity, MAG_CONDUCTIVITY_SCALAR, Global::REAL);
        assert(conductivity != 0.0);
        // PtrCoefFct coeff = CoefFunction::Generate(mp_, Global::REAL, lexical_cast<std::string>(conductivity));
        PtrCoefFct coeff = materials_[actRegion]->GetScalCoefFnc(MAG_CONDUCTIVITY_SCALAR, Global::REAL);

        // Create the integrators
        BaseBDBInt *velocityStiff = NULL;

        // ConvectiveOperator doesn't work with FeHCurl, works at the moment just with FeH1, I am working on it
        if (isComplex_)
        {
          EXCEPTION("MagEdgePDE: VelocityStiffInt not implemented for complex case!");
        }
        else
        {
          if (dim_ == 2)
            velocityStiff = new ABInt<>(new IdentityOperator<FeHCurl, 2, 1>(), new CurlOperatorMag<FeHCurl, 2, Double>(), coeff, 1.0, coefUpdateGeo);
          else
            velocityStiff = new ABInt<>(new IdentityOperator<FeHCurl, 3, 1>(), new CurlOperatorMag<FeHCurl, 3, Double>(), coeff, 1.0, coefUpdateGeo);
        }
        assert(velocityStiff != NULL);

        velocityStiff->SetBCoefFunctionOpB(VelocityCoef_);
        velocityStiff->SetName("VelocityStiffInt");
        velocityInts_[actRegion] = velocityStiff;

        BiLinFormContext *VelocityContextStiff = new BiLinFormContext(velocityStiff, STIFFNESS);
        VelocityContextStiff->SetEntities(actSDList, actSDList);
        VelocityContextStiff->SetFeFunctions(feFunctions_[MAG_POTENTIAL], feFunc);
        assemble_->AddBiLinearForm(VelocityContextStiff);
      } // end velocityId
    } // end for regions
  }

  LinearForm *MagEdgePDE::GetCurrentDensityInt(Double factor, PtrCoefFct coef, std::string coilId)
  {
    LinearForm *ret = NULL;
    // -------------------
    //  EDGE FORMULATION
    // -------------------
    // ===  3D CASE ===
    if (dim_ == 3)
    {
      if (isComplex_)
        ret = new BUIntegrator<Complex>(new IdentityOperator<FeHCurl, 3, 1, Complex>(), factor, coef, updatedGeo_);
      else
        ret = new BUIntegrator<Double>(new IdentityOperator<FeHCurl, 3, 1, Double>(), factor, coef, updatedGeo_);
    }
    else
    {
      // ===  2D / AXI CASE ===
      EXCEPTION("MagEdgePDE is just implemented for 3D setups!");
    }
    return ret;
  }

  LinearForm *MagEdgePDE::GetRHSMagnetizationInt(Double factor, PtrCoefFct rhsMag, bool fullEvaluation)
  {
    LinearForm *lin = NULL;
    bool coefUpdateGeo = true;

    // -------------------
    //  EDGE FORMULATION
    // -------------------
    // ===  3D CASE ===
    if (dim_ == 3)
    {
      if (isComplex_)
      {
        lin = new BUIntegrator<Complex>(new CurlOperator<FeHCurl, 3, Complex>(),
                                        Complex(factor), rhsMag, coefUpdateGeo, fullEvaluation);
      }
      else
      {
        lin = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3, Double>(),
                                       factor, rhsMag, coefUpdateGeo, fullEvaluation);
      }
    }
    else
    {
      // ===  2D / AXI CASE ===
      EXCEPTION("MagEdgePDE is just implemented for 3D setups!");
    }
    return lin;
  }

  BaseBDBInt *MagEdgePDE::GeHystStiffInt(Double factor, PtrCoefFct tensorReluctivity)
  {
    BaseBDBInt *stiffInt = NULL;

    // -------------------
    //  EDGE FORMULATION
    // -------------------
    // ===  3D CASE ===
    if (dim_ == 3)
    {
      stiffInt = new BDBInt<>(new CurlOperator<FeHCurl, 3, Double>(), tensorReluctivity, factor, updatedGeo_);
    }
    else
    {
      // ===  2D / AXI CASE ===
      EXCEPTION("MagEdgePDE is just implemented for 3D setups!");
    }

    return stiffInt;
  }

  void MagEdgePDE::DefineNcIntegrators()
  {
    StdVector<NcInterfaceInfo>::iterator ncIt = ncInterfaces_.Begin(), endIt = ncInterfaces_.End();
    for (; ncIt != endIt; ++ncIt)
    {
      switch (ncIt->type)
      {
      case NC_MORTAR:
        EXCEPTION("No Mortar nonconforming interface for magnetic PDE with edge elements.\n"
                  "Try using H1 nodal elements in MagneticPDE")
        break;
      case NC_NITSCHE:
        if (dim_ == 2)
          EXCEPTION("MagEdgePDE only works for 3D geometry!")
        else if (analysistype_ == MULTIHARMONIC)
        {
          DefineNitscheCoupling<3, 1>(MAG_POTENTIAL, *ncIt, reluc_);
        }
        else
        {
          DefineNitscheCoupling<3, 1>(MAG_POTENTIAL, *ncIt);
        }
        break;
      default:
        EXCEPTION("Unknown type of ncInterface");
        break;
      }
    }
  }

  void MagEdgePDE::DefineRhsLoadIntegrators()
  {
    LOG_TRACE(magEdgePde) << "Defining rhs load integrators for magEdgePDE";

    // Get FESpace and FeFunction of mechanical displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[MAG_POTENTIAL];

    StdVector<shared_ptr<EntityList>> ent;
    StdVector<PtrCoefFct> coef;
    LinearForm *lin = NULL;
    StdVector<std::string> vecDofNames = myFct->GetResultInfo()->dofNames;
    bool coefUpdateGeo = true;

    double magStrictFactor = 1.0;
    //    if ( isMagnetoStrictCoupled_ == true ){
    //      magStrictFactor = -1.0;
    //    }

    // ===============================================================
    // ENERGY-BASED HYSTERETIC A-FORMULATION 3D (start)
    // ===============================================================
    LinearForm *linearform_h_nl_curln = NULL;
    LinearForm *linearform_h_lin_curln = NULL;
    LinearForm *linearform_A_A = NULL;
    RegionIdType actRegion;
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

      PtrCoefFct field_intensity_nl = NULL;
      field_intensity_nl = matModelCoefm_[actRegion]; //matModelCoef_;
      // ===============================================================================================
      // NONLINEAR CASE AND NONLINEAR REGION: \int H(B) \curlN' (start)
      // ===============================================================================================
      if (nonLinTypes.Find(PERMEABILITY) != -1) {
        if (modelName_ == "invEBHysteresisModel") { // use inverse energy-based Hysteresis model
          if (dim_ == 3){
            linearform_h_nl_curln = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3, Double>(), (-1.0), field_intensity_nl, updatedGeo_);
            linearform_A_A = new BUIntegrator<Double>(new IdentityOperator<FeHCurl, 3, 1, Double>(), (0.0), GetCoefFct(MAG_POTENTIAL), updatedGeo_);
          } else {
            EXCEPTION("For Edge elements, only 3D is possible!")
          }

          linearform_h_nl_curln->SetName("(H(B),curl N'): nonlinear problem; nonlinear subregion RHS");
          linearform_h_nl_curln->SetSolDependent();
          LinearFormContext *ctx = new LinearFormContext(linearform_h_nl_curln);
          ctx->SetEntities(actSDList);
          ctx->SetFeFunction(myFct);
          assemble_->AddLinearForm(ctx);
          linearform_A_A->SetName("(A, N'): nonlinear problem; linear subregion RHS");
          linearform_A_A->SetSolDependent();
          ctx = new LinearFormContext(linearform_A_A);
          ctx->SetEntities(actSDList);
          ctx->SetFeFunction(myFct);
          assemble_->AddLinearForm(ctx);
        }
      } else {
        if (modelName_ == "invEBHysteresisModel") { // NONLINEAR CASE BUT LINEAR REGION: \int H curlA'
          if (dim_ == 3) {
            linearform_h_lin_curln = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3, Double>(), (-1.0), GetCoefFct(MAG_FIELD_INTENSITY), updatedGeo_);
            linearform_A_A = new BUIntegrator<Double>(new IdentityOperator<FeHCurl, 3, 1, Double>(), (-1.0), GetCoefFct(MAG_POTENTIAL), updatedGeo_);
          } else {
            EXCEPTION("For Edge elements, only 3D is possible!")
          }
          linearform_h_lin_curln->SetName("(H,curl N'): nonlinear problem; linear subregion RHS");
          linearform_h_lin_curln->SetSolDependent();
          LinearFormContext *ctx = new LinearFormContext(linearform_h_lin_curln);
          ctx->SetEntities(actSDList);
          ctx->SetFeFunction(myFct);
          assemble_->AddLinearForm(ctx);
          linearform_A_A->SetName("(A, N'): nonlinear problem; linear subregion RHS");
          linearform_A_A->SetSolDependent();
          ctx = new LinearFormContext(linearform_A_A);
          ctx->SetEntities(actSDList);
          ctx->SetFeFunction(myFct);
          assemble_->AddLinearForm(ctx);
        }
      }
      // ===============================================================================================
      // NONLINEAR CASE AND NONLINEAR REGION: \int H(B) \curlN' (end)
      // ===============================================================================================
    }
    // ===============================================================
    // ENERGY-BASED HYSTERETIC A-FORMULATION 3D (end)
    // ===============================================================

    // =================================
    //  Magnetization
    // =================================
    // magnetization should contain ALL regions; hysteresisCoefs only those where hysteresis is defined
    // and mRHSRegions_ only those where constant magnetization is prescribed
    std::map<RegionIdType, PtrCoefFct> regionCoefs = magnetization_->GetRegionCoefs();
    std::map<RegionIdType, PtrCoefFct> regionHystCoefs = hysteresisCoefs_->GetRegionCoefs();
    std::map<RegionIdType, shared_ptr<CoefFunction>>::iterator it;
    for (it = regionCoefs.begin(); it != regionCoefs.end(); it++)
    {

      // get regionIdType
      RegionIdType curReg = it->first;
      PtrCoefFct curHystCoef = regionHystCoefs[curReg];
      PtrCoefFct curFixedMagCoef = mRHSRegions_[curReg];

      // get SDList
      shared_ptr<ElemList> actSDList(new ElemList(ptGrid_));
      actSDList->SetRegion(curReg);

      shared_ptr<CoefFunction> rhsMag;
      // set fullevaluation to trigger evaluation at each integration point
      // the nonlinear parameter "evaluation depth" determines if each
      // integration point gets mapped to midpoint (> fullevaluation = false)
      // or if hyst operator really is evaluated at the actual int. point
      bool fullevaluation = true;
      bool forceSolDependent = false;

      if (curHystCoef != NULL)
      {
        std::cout << "Hysteresis region found" << std::endl;
        // NEW: we do not pass the hysteresis coefficient function
        // directly but instead a special class that returns the
        // correctly weighted term
        // even though, we have a similar function for output,
        // we need a separate coefFunction here as rhsMag might
        // be evaluated at another timestep/interation step as outputMag
        // shared_ptr<CoefFunction> rhsMag = it->second->GenerateRHSCoefFnc("MagMagnetization");
        rhsMag = curHystCoef->GenerateRHSCoefFnc("MagMagnetization");
        forceSolDependent = true;
        mRHSRegions_[curReg] = rhsMag;
      }
      else if (curFixedMagCoef != NULL)
      {
        // std::cout << "Use constant magnetization " << std::endl;
        rhsMag = mRHSRegions_[curReg];
      }
      else
      {
        // std::cout << "Neither hysteresis nor constant magnetization prescribed " << std::endl;
        continue;
      }

      lin = GetRHSMagnetizationInt(magStrictFactor, rhsMag, fullevaluation);

      lin->SetName("rhs_magnetization");
      if (forceSolDependent)
      {
        lin->SetSolDependent();
      }
      LinearFormContext *ctx = new LinearFormContext(lin);
      ctx->SetEntities(actSDList);
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      // Add entity list will add nothing, if entities were already assigned
      myFct->AddEntityList(actSDList);
    }
    
    // ==================
    //  FLUX DENSITY
    // ==================
    LOG_DBG(magEdgePde) << "Reading prescribed flux density";

    ReadRhsExcitation("fluxDensity", vecDofNames, ResultInfo::VECTOR, isComplex_,
                      ent, coef, coefUpdateGeo);
    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
          ent[i]->GetType() == EntityList::SURF_ELEM_LIST)
      {
        EXCEPTION("Prescribed magnetic flux density can only be specified in a volume!")
      }

      bool forceSolDependency = false;
      shared_ptr<CoefFunction> rhsCoefFunctionToBeUsed;

      // check for hysteresis on current entity
      std::map<RegionIdType, PtrCoefFct> regionCoefs = hysteresisCoefs_->GetRegionCoefs();

      if (regionCoefs.count(ent[i]->GetRegion()))
      {
        //        std::cout << "RHS-Hysteresis" << std::endl;
        // HYST CASE > compare with MagneticPDE.cc
        // GenerateRHSCoefFnc("MagPolarization",coef[i]) returns coef[i] - MagPolarization = B - J
        // together with the reluctivityTensor we get M on rhs
        rhsCoefFunctionToBeUsed = regionCoefs[ent[i]->GetRegion()]->GenerateRHSCoefFnc("MagPolarization", coef[i]);
        forceSolDependency = true;

        if (isComplex_)
        {
          lin = new BDUIntegrator<CurlOperator<FeHCurl, 3, Double>, Complex>(Complex(magStrictFactor),
                                                                             rhsCoefFunctionToBeUsed, relucTensor_, coefUpdateGeo);
        }
        else
        {
          lin = new BDUIntegrator<CurlOperator<FeHCurl, 3, Double>, Double>(magStrictFactor, rhsCoefFunctionToBeUsed, relucTensor_, coefUpdateGeo);
        }
      }
      else
      {
        // classical implementation
        Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
        PtrCoefFct factor = CoefFunction::Generate(mp_, part,
                                                   CoefXprBinOp(mp_, reluc_, coef[i], CoefXpr::OP_MULT));

        if (isComplex_)
        {
          lin = new BUIntegrator<Complex>(new CurlOperator<FeHCurl, 3, Complex>(),
                                          Complex(magStrictFactor), factor, coefUpdateGeo);
        }
        else
        {
          lin = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3, Double>(),
                                         magStrictFactor, factor, coefUpdateGeo);
        }
      }

      lin->SetName("FluxIntegrator");
      if (forceSolDependency)
      {
        lin->SetSolDependent();
      }

      LinearFormContext *ctx = new LinearFormContext(lin);
      ctx->SetEntities(ent[i]);
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);

      bRHSRegions_[ent[i]->GetRegion()] = coef[i];
    } // for

    // ==================
    //  FIELD INTENSITY
    // ==================
    LOG_DBG(magEdgePde) << "Reading prescribed field intensity";

    ReadRhsExcitation("fieldIntensity", vecDofNames, ResultInfo::VECTOR, isComplex_,
                      ent, coef, coefUpdateGeo);
    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST ||
          ent[i]->GetType() == EntityList::SURF_ELEM_LIST)
      {
        EXCEPTION("Prescribed magnetic field intensity can only be specified in a volume!")
      }

      if (isComplex_)
      {
        lin = new BUIntegrator<Complex>(new CurlOperator<FeHCurl, 3, Complex>(),
                                        Complex(1.0), coef[i], coefUpdateGeo);
      }
      else
      {
        lin = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3, Double>(),
                                       1.0, coef[i], coefUpdateGeo);
      }
      lin->SetName("FieldIntensityIntegrator");
      LinearFormContext *ctx = new LinearFormContext(lin);
      ctx->SetEntities(ent[i]);
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);

    } // for

    //    std::cout << "DONE Define RHS Integrators" << std::endl;
  }

  // ======================================================
  // TIME-STEPPING SECTION > MagBasePDE
  // ======================================================

  void MagEdgePDE::DefinePrimaryResults()
  {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    // === MAGNETIC VECTOR POTENTIAL ===
    shared_ptr<ResultInfo> potInfo(new ResultInfo);
    potInfo->resultType = MAG_POTENTIAL;
    potInfo->dofNames = vecComponents;
    potInfo->unit = "Vs/m";
    potInfo->definedOn = ResultInfo::ELEMENT;
    potInfo->entryType = ResultInfo::VECTOR;

    feFunctions_[MAG_POTENTIAL]->SetResultInfo(potInfo);
    DefineFieldResult(feFunctions_[MAG_POTENTIAL], potInfo);

    // === COIL CURRENT ===
    if (hasVoltCoils_)
    {
      shared_ptr<ResultInfo> currentInfo(new ResultInfo);
      currentInfo->resultType = COIL_CURRENT;
      currentInfo->dofNames = "";
      currentInfo->unit = "A";
      currentInfo->definedOn = ResultInfo::COIL;
      currentInfo->entryType = ResultInfo::SCALAR;

      feFunctions_[COIL_CURRENT]->SetResultInfo(currentInfo);
      DefineFieldResult(feFunctions_[COIL_CURRENT], currentInfo);
    }

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MAG_POTENTIAL] = "fluxParallel";
    idbcSolNameMap_[MAG_POTENTIAL] = "potential";

    // === PERMEABILITY ===
    shared_ptr<ResultInfo> permeability(new ResultInfo);
    permeability->resultType = MAG_ELEM_PERMEABILITY;
    permeability->dofNames = "";
    permeability->unit = "Vs/Am";
    permeability->definedOn = ResultInfo::ELEMENT;
    permeability->entryType = ResultInfo::SCALAR;

    if (analysistype_ == MULTIHARMONIC)
    {
      shared_ptr<CoefFunctionMulti> permFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, true));
      matCoefs_[MAG_ELEM_PERMEABILITY] = permFct;
      DefineFieldResult(permFct, permeability);
    }
    else
    {
      // In multiharmonic analysis we have complex permeability
      shared_ptr<CoefFunctionMulti> permFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, false));
      matCoefs_[MAG_ELEM_PERMEABILITY] = permFct;
      DefineFieldResult(permFct, permeability);
    }

    // creates the velocity
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
        if (myParam_->Get("formulation")->As<std::string>() == "A-V")
        {
          WARN("No implementation for axi-symmetric model in the MagEdge Case");
        };
      }
      else
      {
        vecDofNames = "x", "y";
      }
    }

    //// === VELOCITY ===
    shared_ptr<ResultInfo> velocity(new ResultInfo);
    velocity->resultType = MEAN_FLUIDMECH_VELOCITY;
    velocity->dofNames = vecDofNames;
    velocity->unit = "m/s";

    velocity->definedOn = ResultInfo::NODE;
    velocity->entryType = ResultInfo::VECTOR;

    VelocityCoef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
    DefineFieldResult(VelocityCoef_, velocity);

    results_.Push_back(velocity);
    availResults_.insert(velocity);
  }

  void MagEdgePDE::DefinePostProcResults()
  {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];

    // === TIME DERIVATIVES OF PRIMARY RESULTS ===
    if (analysistype_ == TRANSIENT || analysistype_ == HARMONIC || analysistype_ == MULTIHARMONIC)
    {
      // === MAGNETIC VECTOR POTENTIAL - 1ST DERIVATIVE ===
      shared_ptr<ResultInfo> aDot(new ResultInfo);
      aDot->resultType = MAG_POTENTIAL_DERIV1;
      aDot->dofNames = vecComponents;
      aDot->unit = "V/m";
      aDot->definedOn = ResultInfo::ELEMENT;
      aDot->entryType = ResultInfo::VECTOR;
      availResults_.insert(aDot);
      DefineTimeDerivResult(MAG_POTENTIAL_DERIV1, 1, MAG_POTENTIAL);

      // === COIL CURRENT OF COILS EXCITED BY VOLTAGE - 1ST DERIVATIVE ===
      if (hasVoltCoils_)
      {
        shared_ptr<ResultInfo> iDot(new ResultInfo);
        iDot->resultType = COIL_CURRENT_DERIV1;
        iDot->dofNames = "";
        iDot->unit = "A/s";
        iDot->definedOn = ResultInfo::COIL;
        iDot->entryType = ResultInfo::SCALAR;
        availResults_.insert(iDot);
        DefineTimeDerivResult(COIL_CURRENT_DERIV1, 1, COIL_CURRENT);
      }
    }

    if (!fluxDensityDefined_)
    {
      DefineMagFluxDensity();
    }

    PtrCoefFct bFunc = this->GetCoefFct(MAG_FLUX_DENSITY);

    // === MAGNETIC NORMAL FLUX DENSITY ===
    shared_ptr<ResultInfo> normFlux(new ResultInfo);
    normFlux->resultType = MAG_NORMAL_FLUX_DENSITY;
    normFlux->dofNames = "";
    normFlux->unit = "Vs/m^2";
    normFlux->entryType = ResultInfo::SCALAR;
    normFlux->definedOn = ResultInfo::ELEMENT;
    shared_ptr<CoefFunctionSurf> sNormFDens;
    sNormFDens.reset(new CoefFunctionSurf(true, 1.0, normFlux));
    DefineFieldResult(sNormFDens, normFlux);
    surfCoefFcts_[sNormFDens] = bFunc;

    // === MAGNETIC_FLUX ===
    shared_ptr<ResultInfo> flux(new ResultInfo);
    flux->resultType = MAG_FLUX;
    flux->dofNames = "";
    flux->unit = "Vs";
    flux->entryType = ResultInfo::SCALAR;
    flux->definedOn = ResultInfo::SURF_REGION;
    shared_ptr<ResultFunctor> fluxFct;
    if (isComplex_)
    {
      fluxFct.reset(new ResultFunctorIntegrate<Complex>(sNormFDens,
                                                        feFct, flux));
    }
    else
    {
      fluxFct.reset(new ResultFunctorIntegrate<Double>(sNormFDens,
                                                       feFct, flux));
    }
    resultFunctors_[MAG_FLUX] = fluxFct;
    availResults_.insert(flux);

    // === MAGNETIC RHS ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = MAG_RHS_LOAD;
    rhs->dofNames = vecComponents;
    rhs->unit = "-";
    rhs->entryType = ResultInfo::VECTOR;
    rhs->definedOn = ResultInfo::ELEMENT;
    rhsFeFunctions_[MAG_POTENTIAL]->SetResultInfo(rhs);
    DefineFieldResult(rhsFeFunctions_[MAG_POTENTIAL], rhs);

    // === RESULTS RELATED TO TIME DERIVATIVES ===
    shared_ptr<CoefFunctionFormBased> jFunc;
    shared_ptr<CoefFunction> jPowerDensFunc;
    if (analysistype_ != STATIC)
    {
      // === EDDY CURRENT DENSITY ===
      shared_ptr<BaseFeFunction> aDotFct =
          timeDerivFeFunctions_[MAG_POTENTIAL_DERIV1];
      shared_ptr<ResultInfo> eddy(new ResultInfo);
      eddy->resultType = MAG_EDDY_CURRENT_DENSITY;
      eddy->dofNames = vecComponents;
      eddy->unit = "A/m^2";
      eddy->definedOn = ResultInfo::ELEMENT;
      eddy->entryType = ResultInfo::VECTOR;

      if (isComplex_)
      {
        jFunc.reset(new CoefFunctionFlux<Complex>(aDotFct, eddy, -1.0));
      }
      else
      {
        jFunc.reset(new CoefFunctionFlux<Double>(aDotFct, eddy, -1.0));
      }
      DefineFieldResult(jFunc, eddy);

      // === EDDY CURRENT (SURFACE RESULT) ===
      shared_ptr<ResultInfo> ec(new ResultInfo());
      ec->resultType = MAG_EDDY_CURRENT;
      ec->dofNames = "";
      ec->unit = "A";
      ec->definedOn = ResultInfo::SURF_REGION;
      ec->entryType = ResultInfo::SCALAR;
      availResults_.insert(ec);

      // first, create normal mapping
      shared_ptr<CoefFunctionSurf> ncd(new CoefFunctionSurf(true, 1.0, ec));
      surfCoefFcts_[ncd] = jFunc;

      // then, integrate values
      shared_ptr<ResultFunctor> eddyCurrentFunc;
      if (isComplex_)
      {
        eddyCurrentFunc.reset(new ResultFunctorIntegrate<Complex>(ncd,
                                                                  feFct, ec));
      }
      else
      {
        eddyCurrentFunc.reset(new ResultFunctorIntegrate<Double>(ncd,
                                                                 feFct, ec));
      }
      resultFunctors_[MAG_EDDY_CURRENT] = eddyCurrentFunc;

      // === EDDY POWER DENSITY ===
      shared_ptr<ResultInfo> epd(new ResultInfo());
      epd->resultType = MAG_EDDY_POWER_DENSITY;
      epd->dofNames = "";
      epd->unit = "W/m^3";
      epd->definedOn = ResultInfo::ELEMENT;
      epd->entryType = ResultInfo::SCALAR;
      shared_ptr<CoefFunctionFormBased> epdFunctor;
      if (isComplex_)
      {
        epdFunctor.reset(new CoefFunctionBdBKernel<Complex>(aDotFct, 1.0));
      }
      else
      {
        epdFunctor.reset(new CoefFunctionBdBKernel<Double>(aDotFct, 1.0));
      }
      DefineFieldResult(epdFunctor, epd);
      massFormCoefs_.insert(epdFunctor);

      // === EDDY POWER ===
      shared_ptr<ResultInfo> ep(new ResultInfo());
      ep->resultType = MAG_EDDY_POWER;
      ep->dofNames = "";
      ep->unit = "W";
      ep->definedOn = ResultInfo::REGION;
      ep->entryType = ResultInfo::SCALAR;
      availResults_.insert(ep);
      shared_ptr<ResultFunctor> epFunctor;
      if (isComplex_)
      {
        epFunctor.reset(new EnergyResultFunctor<Complex>(aDotFct, ep, 1.0));
      }
      else
      {
        epFunctor.reset(new EnergyResultFunctor<Double>(aDotFct, ep, 1.0));
      }
      resultFunctors_[MAG_EDDY_POWER] = epFunctor;
      massFormFunctors_.insert(epFunctor);

      // === ELECTRIC FIELD INTENSITY ===
      shared_ptr<ResultInfo> elecIntens(new ResultInfo);
      elecIntens->resultType = ELEC_FIELD_INTENSITY;
      elecIntens->SetVectorDOFs(dim_, isaxi_);
      elecIntens->dofNames = vecComponents;
      elecIntens->unit = "V/m";
      elecIntens->definedOn = ResultInfo::ELEMENT;
      elecIntens->entryType = ResultInfo::VECTOR;

      // assemble coefficient function E = - dA/dt
      PtrCoefFct eIFuncTmp = CoefFunction::Generate(mp_, part,
                                                    CoefXprBinOp(mp_, aDotFct, aDotFct, CoefXpr::OP_MULT));
      PtrCoefFct eIFunc = CoefFunction::Generate(mp_, part,
                                                 CoefXprBinOp(mp_, "-1.0", aDotFct, CoefXpr::OP_MULT));
      DefineFieldResult(eIFunc, elecIntens);

      // === COIL LINKED FLUX - 1ST DERIVATIVE, CORRESPONDS TO INDUCED VOLTAGE ===
      shared_ptr<ResultInfo> psiDotRes(new ResultInfo());
      psiDotRes->resultType = COIL_INDUCED_VOLTAGE;
      psiDotRes->dofNames = "";
      psiDotRes->unit = "V";
      psiDotRes->definedOn = ResultInfo::COIL;
      psiDotRes->entryType = ResultInfo::SCALAR;

      availResults_.insert(psiDotRes);
      shared_ptr<ResultFunctor> psiDotFunc;
      shared_ptr<CoefFunctionMulti> psiDotDens(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1,
                                                                     isComplex_));
      if (isComplex_)
      {
        psiDotFunc.reset(new ResultFunctorIntegrate<Complex>(psiDotDens, aDotFct, psiDotRes));
      }
      else
      {
        psiDotFunc.reset(new ResultFunctorIntegrate<Double>(psiDotDens, aDotFct, psiDotRes));
      }
      resultFunctors_[COIL_INDUCED_VOLTAGE] = psiDotFunc;
      // it is an integrated result but we need to save the coef function
      // somewhere for the finalization
      fieldCoefs_[COIL_INDUCED_VOLTAGE] = psiDotDens;
    }

    // === COIL CURRENT DENSITY ===
    shared_ptr<ResultInfo> ccd(new ResultInfo);
    ccd->resultType = MAG_COIL_CURRENT_DENSITY;
    ccd->dofNames = vecComponents;
    ccd->unit = "A/m^2";
    ccd->definedOn = ResultInfo::ELEMENT;
    ccd->entryType = ResultInfo::VECTOR;
    availResults_.insert(ccd);
    shared_ptr<CoefFunctionMulti> ccdCoef(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1,
                                                                isComplex_));
    DefineFieldResult(ccdCoef, ccd);

    // === TOTAL CURRENT DENSITY ===
    shared_ptr<ResultInfo> tcd(new ResultInfo);
    tcd->resultType = MAG_TOTAL_CURRENT_DENSITY;
    tcd->dofNames = vecComponents;
    tcd->unit = "A/m^2";
    tcd->definedOn = ResultInfo::ELEMENT;
    tcd->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionMulti> tcdCoef(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1,
                                                                isComplex_));
    DefineFieldResult(tcdCoef, tcd);

    GenerateLorentzForceResults(vecComponents, tcdCoef, bFunc, part, feFct);

    // === RELUCTIVITY  ===
    shared_ptr<ResultInfo> reluc(new ResultInfo);
    reluc->resultType = MAG_ELEM_RELUCTIVITY;
    reluc->dofNames = "";
    reluc->unit = "Am/Vs";
    reluc->definedOn = ResultInfo::ELEMENT;
    reluc->entryType = ResultInfo::SCALAR;
    DefineFieldResult(reluc_, reluc);

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
    if (modelName_ == "invEBHysteresisModel")
    {
      // The recipe how to actually evaluate H, is defined in FinalizePostProcResults()
      shared_ptr<CoefFunctionMulti> hFunc(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
      DefineFieldResult(hFunc, magIntens);
    }
    else
    {
      DefineFieldResult(fieldIntensity_, magIntens);
      availResults_.insert(magIntens);
    }
    // =====================================================
    // MAG_FIELD_INTENSITY (END)
    // =====================================================

    shared_ptr<ResultInfo> magJ(new ResultInfo);
    magJ->resultType = MAG_POLARIZATION;
    magJ->SetVectorDOFs(dim_, isaxi_);
    magJ->unit = "Vs/m^2";
    magJ->definedOn = ResultInfo::ELEMENT;
    magJ->entryType = ResultInfo::VECTOR;

    DefineFieldResult(polarization_, magJ);
    availResults_.insert(magJ);

    shared_ptr<ResultInfo> magM(new ResultInfo);
    magM->resultType = MAG_MAGNETIZATION;
    magM->SetVectorDOFs(dim_, isaxi_);
    magM->unit = "A/m";
    magM->definedOn = ResultInfo::ELEMENT;
    magM->entryType = ResultInfo::VECTOR;

    DefineFieldResult(magnetization_, magM);
    availResults_.insert(magM);

    GenerateMaxwellForce(vecComponents, bFunc, feFct);

    GenerateVWPForce(vecComponents, bFunc, feFct);

    // === MAGNETIC ENERGY ===
    shared_ptr<ResultInfo> energy(new ResultInfo);
    energy->resultType = MAG_ENERGY;
    energy->dofNames = "";
    energy->unit = "Ws";
    energy->definedOn = ResultInfo::REGION;
    energy->entryType = ResultInfo::SCALAR;
    availResults_.insert(energy);
    shared_ptr<ResultFunctor> energyFunc;
    if (isComplex_)
    {
      energyFunc.reset(new EnergyResultFunctor<Complex>(feFct, energy, 0.5));
    }
    else
    {
      energyFunc.reset(new EnergyResultFunctor<Double>(feFct, energy, 0.5));
    }
    resultFunctors_[MAG_ENERGY] = energyFunc;
    stiffFormFunctors_.insert(energyFunc);

    // === COIL LINKED FLUX  ===
    shared_ptr<ResultInfo> psiRes(new ResultInfo());
    psiRes->resultType = COIL_LINKED_FLUX;
    psiRes->dofNames = "";
    psiRes->unit = "Vs/m^2";
    psiRes->definedOn = ResultInfo::COIL;
    psiRes->entryType = ResultInfo::SCALAR;

    availResults_.insert(psiRes);
    shared_ptr<ResultFunctor> psiFunc;
    shared_ptr<CoefFunctionMulti> psiDens(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1,
                                                                isComplex_));
    if (isComplex_)
    {
      psiFunc.reset(new ResultFunctorIntegrate<Complex>(psiDens, feFct, psiRes));
    }
    else
    {
      psiFunc.reset(new ResultFunctorIntegrate<Double>(psiDens, feFct, psiRes));
    }
    resultFunctors_[COIL_LINKED_FLUX] = psiFunc;
    // it is an integrated result but we need to save the coef function
    // somewhere for the finalization
    fieldCoefs_[COIL_LINKED_FLUX] = psiDens;

    // === CORE LOSS DENSITY ===
    shared_ptr<ResultInfo> cldRes(new ResultInfo());
    cldRes->resultType = MAG_CORE_LOSS_DENSITY;
    cldRes->dofNames = "";
    cldRes->unit = "W/kg";
    cldRes->definedOn = ResultInfo::ELEMENT;
    cldRes->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> coreLossDensCoef(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
    DefineFieldResult(coreLossDensCoef, cldRes);

    // === CORE LOSS ===
    shared_ptr<ResultInfo> clRes(new ResultInfo());
    clRes->resultType = MAG_CORE_LOSS;
    clRes->dofNames = "";
    clRes->unit = "W";
    clRes->definedOn = ResultInfo::REGION;
    clRes->entryType = ResultInfo::SCALAR;
    // DefineFieldResult( coreLossCoef, clRes );
    availResults_.insert(clRes);
    shared_ptr<ResultFunctor> coreLossFunc;
    shared_ptr<CoefFunctionMulti> coreLossCoef(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
    if (isComplex_)
    {
      coreLossFunc.reset(new ResultFunctorIntegrate<Complex>(coreLossCoef, feFct, clRes));
    }
    else
    {
      coreLossFunc.reset(new ResultFunctorIntegrate<Double>(coreLossCoef, feFct, clRes));
    }
    resultFunctors_[MAG_CORE_LOSS] = coreLossFunc;
    // it is an integrated result but we need to save the coef function
    // somewhere for the finalization
    fieldCoefs_[MAG_CORE_LOSS] = coreLossCoef;

    // === JOULE LOSS Power DENSITY INTEGRATED OVER PERIOD	 ===
    shared_ptr<ResultInfo> jld(new ResultInfo);
    jld->resultType = MAG_JOULE_LOSS_POWER_DENSITY;
    jld->dofNames = "";
    jld->unit = "W/m^3";
    jld->definedOn = ResultInfo::ELEMENT;
    jld->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> jldCoef(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
    DefineFieldResult(jldCoef, jld);

    // optimization results are provided in DesignSpace::ExtractResults()
    DefineFieldResult(PSEUDO_DENSITY, ResultInfo::SCALAR, ResultInfo::ELEMENT, "", true);
    DefineFieldResult(PHYSICAL_PSEUDO_DENSITY, ResultInfo::SCALAR, ResultInfo::ELEMENT, "", true);
    DefineFieldResult(RHS_PSEUDO_DENSITY, ResultInfo::SCALAR, ResultInfo::ELEMENT, "", true);
    DefineFieldResult(PHYSICAL_RHS_PSEUDO_DENSITY, ResultInfo::SCALAR, ResultInfo::ELEMENT, "", true);

    if (analysistype_ != STATIC)
    {
      shared_ptr<ResultInfo> jldN(new ResultInfo);
      jldN->resultType = MAG_JOULE_LOSS_POWER_DENSITY_ON_NODES;
      jldN->dofNames = "";
      jldN->unit = "W/m^3";
      jldN->definedOn = ResultInfo::NODE;
      jldN->entryType = ResultInfo::SCALAR;
      shared_ptr<CoefFunctionMulti> jldNCoef(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
      DefineFieldResult(jldNCoef, jldN);

      // === JOULE LOSS POWER ===
      shared_ptr<ResultInfo> jldRes(new ResultInfo());
      jldRes->resultType = MAG_JOULE_LOSS_POWER;
      jldRes->dofNames = "";
      jldRes->unit = "W";
      jldRes->definedOn = ResultInfo::REGION;
      jldRes->entryType = ResultInfo::SCALAR;
      availResults_.insert(jldRes);
      shared_ptr<ResultFunctor> jldFunc;
      shared_ptr<CoefFunctionMulti> coreLossCoef(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
      if (isComplex_)
      {
        jldFunc.reset(new ResultFunctorIntegrate<Complex>(jldCoef, feFct, jldRes));
      }
      else
      {
        jldFunc.reset(new ResultFunctorIntegrate<Double>(jldCoef, feFct, jldRes));
      }
      resultFunctors_[MAG_JOULE_LOSS_POWER] = jldFunc;
      // it is an integrated result but we need to save the coef function
      // somewhere for the finalization
      fieldCoefs_[MAG_CORE_LOSS] = coreLossCoef;
    }
  }

  void MagEdgePDE::FinalizePostProcResults()
  {

    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();

    // Initialized results involving coils

    // === EDDY CURRENT DENSITY ===
    shared_ptr<CoefFunctionFormBased> jEddyCoef =
        dynamic_pointer_cast<CoefFunctionFormBased>(fieldCoefs_[MAG_EDDY_CURRENT_DENSITY]);
    std::map<RegionIdType, BaseBDBInt *>::iterator massIt = massInts_.begin();
    for (; massIt != massInts_.end(); ++massIt)
    {
      RegionIdType region = massIt->first;
      BaseBDBInt *massInt = massIt->second;

      // only assign region to jEddy, if it not a "regularized" region, i.e.
      // only regions with "physical" conductivity get assigned
      if (regularizedRegions_.find(region) == regularizedRegions_.end())
      {
        jEddyCoef->AddIntegrator(massInt, region);
      }
    }

    // === COIL CURRENT DENSITY ===
    shared_ptr<CoefFunctionMulti> ccdCoef =
        dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_COIL_CURRENT_DENSITY]);
    // loop over all coil coefficients and add contribution to coef
    std::map<RegionIdType, PtrCoefFct>::iterator coilIt = coilCurrentDens_.begin();
    for (; coilIt != coilCurrentDens_.end(); ++coilIt)
    {
      ccdCoef->AddRegion(coilIt->first, coilIt->second);
    }

    // === TOTAL CURRENT DENSITY ===
    PtrCoefFct jEddy = GetCoefFct(MAG_EDDY_CURRENT_DENSITY);
    shared_ptr<CoefFunctionMulti> tcdCoef =
        dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_TOTAL_CURRENT_DENSITY]);
    // loop over all regions and assemble total current density:
    //  - if region is coil -> take coil current
    //  - if region is no coil and analyis is transient/harmonic -> eddy

    StdVector<RegionIdType>::iterator regIt = regions_.Begin();
    for (; regIt != regions_.End(); ++regIt)
    {
      RegionIdType actRegion = *regIt;
      if (coilCurrentDens_.find(actRegion) != coilCurrentDens_.end())
      {
        // region is a coil
        tcdCoef->AddRegion(actRegion, coilCurrentDens_[actRegion]);
      }
      else
      {
        // region is no coil
        if (analysistype_ == TRANSIENT || analysistype_ == HARMONIC || analysistype_ == MULTIHARMONIC)
        {
          tcdCoef->AddRegion(actRegion, jEddy);
        }
      }
    }

    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;

    // === H field ===
    if (!isHysteresis_)
    {
      // magnetization already defined in InitMagnetization (see MagBasePDE.cc)
      // however, this function does not consider bRHSRegions_
      // > if B is set on rhs, redefine H by setting new flag
      bool allowReplacement = true;
      // assemble coefficient function field intensity = reluctivity * (flux density - remanence)
      // the remanence is the RHS load flux density
      shared_ptr<CoefFunctionMulti> hCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_FIELD_INTENSITY]);
      regIt = regions_.Begin();
      for (; regIt != regions_.End(); ++regIt)
      {
        if (bRHSRegions_.find(*regIt) != bRHSRegions_.end())
        {

          if (mRHSRegions_.find(*regIt) != mRHSRegions_.end())
          {
            WARN("Cannot prescribe RHS Flux and fixed magnetization at the same time. RHS Flux will overwrite fixed magnetization");
          }

          PtrCoefFct bPM = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, this->GetCoefFct(MAG_FLUX_DENSITY), bRHSRegions_[*regIt], CoefXpr::OP_SUB));
          PtrCoefFct hPM = CoefFunction::Generate(mp_, part, CoefXprVecScalOp(mp_, bPM, reluc_, CoefXpr::OP_MULT));
          hCoef->AddRegion(*regIt, hPM, allowReplacement);
        }
      }
    }

    // === INDUCED VOLTAGE ===
    // This code assembles a CoefFunctionMulti containing all coil regions, which
    // is integrated by a ResultFunctorIntegrate. Ref: M. Kaltenbacher, Numer. Sim.
    // of. Mech. Sens. and Act., 2nd edition, p. 212, Eq. (7.19)
    // Every coil part has an own CoefFunction for the current density, which is
    // why the coil regions are searched to find the corresponding part.
    // It is assumed implicitly that a part cannot contain a region contained by
    // another part. If that happens, the CoefFunctionMulti throws.
    if (analysistype_ == TRANSIENT || analysistype_ == HARMONIC || analysistype_ == MULTIHARMONIC)
    {
      PtrCoefFct temp = GetCoefFct(COIL_INDUCED_VOLTAGE);
      shared_ptr<CoefFunctionMulti> psiDotDens =
          dynamic_pointer_cast<CoefFunctionMulti>(temp);
      CoilRegionMap::iterator coilRegsIt = coilRegions_.begin();
      for (; coilRegsIt != coilRegions_.end(); ++coilRegsIt)
      {
        std::map<RegionIdType, shared_ptr<Coil::Part>>::iterator partIt =
            coilRegsIt->second->parts_.find(coilRegsIt->first);
        CoefXprVecScalOp eJscaledOp = CoefXprVecScalOp(mp_, partIt->second->jUnitVec,
                                                       boost::lexical_cast<std::string>(partIt->second->wireCrossSect), CoefXpr::OP_DIV);
        PtrCoefFct eJscaled = CoefFunction::Generate(mp_, part, eJscaledOp);
        CoefXprBinOp integrandOp = CoefXprBinOp(mp_, eJscaled,
                                                GetCoefFct(MAG_POTENTIAL_DERIV1), CoefXpr::OP_MULT);
        PtrCoefFct integrand = CoefFunction::Generate(mp_, part, integrandOp);
        psiDotDens->AddRegion(coilRegsIt->first, integrand);
      }
    }

    // === COIL LINKED FLUX ===
    // same as for induced voltage, but with the vector potential instead
    // of the first time derivative of the vector potential
    PtrCoefFct temp = GetCoefFct(COIL_LINKED_FLUX);
    shared_ptr<CoefFunctionMulti> psiDotDens =
        dynamic_pointer_cast<CoefFunctionMulti>(temp);
    CoilRegionMap::iterator coilRegsIt = coilRegions_.begin();
    for (; coilRegsIt != coilRegions_.end(); ++coilRegsIt)
    {
      std::map<RegionIdType, shared_ptr<Coil::Part>>::iterator partIt =
          coilRegsIt->second->parts_.find(coilRegsIt->first);
      CoefXprVecScalOp eJscaledOp = CoefXprVecScalOp(mp_, partIt->second->jUnitVec,
                                                     boost::lexical_cast<std::string>(partIt->second->wireCrossSect), CoefXpr::OP_DIV);
      PtrCoefFct eJscaled = CoefFunction::Generate(mp_, part, eJscaledOp);
      CoefXprBinOp integrandOp = CoefXprBinOp(mp_, eJscaled,
                                              GetCoefFct(MAG_POTENTIAL), CoefXpr::OP_MULT);
      PtrCoefFct integrand = CoefFunction::Generate(mp_, part, integrandOp);
      psiDotDens->AddRegion(coilRegsIt->first, integrand);
    }

    // === CORE LOSS DENSITY ===
    // This is a "density" per kg, not per m^3. That's why we cannot integrate it over
    // the volume to get the core loss (see below).

    // TODO Core loss computation is disabled temporarily for harmonic analysis because otherwise we
    //  get an error in AddRegion, due to the multiplication of a real- and a complex-valued
    //  coefficient function. CHECK THIS !!!
    if ((analysistype_ != HARMONIC) && (analysistype_ != MULTIHARMONIC))
    {
      shared_ptr<CoefFunctionMulti> cldCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_CORE_LOSS_DENSITY]);
      BaseMaterial *actMat = NULL;
      regIt = regions_.Begin();
      for (; regIt != regions_.End(); ++regIt)
      {
        RegionIdType actRegion = *regIt;
        actMat = materials_[actRegion];
        PtrCoefFct coreLossFct = actMat->GetScalCoefFncNonLin(MAG_CORE_LOSS_PER_MASS, Global::REAL, GetCoefFct(MAG_FLUX_DENSITY));
        cldCoef->AddRegion(actRegion, coreLossFct);
      }

      // === CORE LOSS ===
      // The core loss per kg has to be multiplied by the density to get it per m^3.
      // Then we can integrate over the volume in order to get the core loss.
      shared_ptr<CoefFunctionMulti> clCoef =
          dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_CORE_LOSS]);
      actMat = NULL;
      regIt = regions_.Begin();
      for (; regIt != regions_.End(); ++regIt)
      {
        RegionIdType actRegion = *regIt;
        actMat = materials_[actRegion];
        // core loss density in W/kg
        PtrCoefFct coreLossFct = actMat->GetScalCoefFncNonLin(MAG_CORE_LOSS_PER_MASS, Global::REAL, GetCoefFct(MAG_FLUX_DENSITY));
        // core loss density in W/m^3, which deserves to be called density
        PtrCoefFct densFct = actMat->GetScalCoefFnc(DENSITY, Global::REAL);
        PtrCoefFct coreLossDensCoef = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, coreLossFct, densFct, CoefXpr::OP_MULT));
        clCoef->AddRegion(actRegion, coreLossDensCoef);
      }
    }

    // === EDDY CURRENT (JOULE) LOSS DENSITY INTEGRATED===
    /*  The electric power averaved over
     *  one period T of the time history
     *    P_mean = 1/T \int_0^T E(t)*J(t) dt
     *  with the electric field E(t) and the total current density J(t).
     *  The general expression is:
     *  omega*Im(A'*J)
     *  where J is the total current density,
     *  and A' is the conjugate complex of the magnetic vector potential.
     */
    if ((analysistype_ == HARMONIC) || (analysistype_ == MULTIHARMONIC))
    {
      shared_ptr<CoefFunctionMulti> eddyLossCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_JOULE_LOSS_POWER_DENSITY]);
      shared_ptr<CoefFunctionMulti> eddyLossCoefN = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_JOULE_LOSS_POWER_DENSITY_ON_NODES]);
      regIt = regions_.Begin();
      // for the sake of simplicity we should real with the total current density
      for (; regIt != regions_.End(); ++regIt)
      {
        RegionIdType actRegion = *regIt;
        // OP_MULT_CONJ computes A in B'
        PtrCoefFct conjAinJ = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, GetCoefFct(MAG_TOTAL_CURRENT_DENSITY), GetCoefFct(MAG_POTENTIAL), CoefXpr::OP_MULT_CONJ));
        PtrCoefFct coefIm = CoefFunction::Generate(mp_, part, CoefXprUnaryOp(mp_, conjAinJ, CoefXpr::OP_IM));
        // frequency factor
        PtrCoefFct coefFreqFactor = NULL;
        if (analysistype_ == MULTIHARMONIC)
        {
          // consider that the quantities (A,J) in multiharmonic are 1/2 of the harmonic ones,
          // the total dissipated power is obtained by summing all the harmonics
          // the negative factor is to get out positive values ???
          coefFreqFactor = CoefFunction::Generate(mp_, part, "-2*pi*f");
        }
        else
        {
          // for harmonic
          coefFreqFactor = CoefFunction::Generate(mp_, part, "-pi*f");
        }
        eddyLossCoef->AddRegion(actRegion, CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, coefFreqFactor, coefIm, CoefXpr::OP_MULT)));
        eddyLossCoefN->AddRegion(actRegion, CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, coefFreqFactor, coefIm, CoefXpr::OP_MULT)));
      }
    }

    if (analysistype_ == TRANSIENT)
    {
      shared_ptr<CoefFunctionMulti> eddyLossCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_JOULE_LOSS_POWER_DENSITY]);
      shared_ptr<CoefFunctionMulti> eddyLossCoefN = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_JOULE_LOSS_POWER_DENSITY_ON_NODES]);
      regIt = regions_.Begin();
      for (; regIt != regions_.End(); ++regIt)
      {
        RegionIdType actRegion = *regIt;
        // J_total * E = J_total * dA/dt
        PtrCoefFct partTmp = CoefFunction::Generate(mp_, part,
                                                    CoefXprBinOp(mp_, GetCoefFct(ELEC_FIELD_INTENSITY), GetCoefFct(MAG_TOTAL_CURRENT_DENSITY), CoefXpr::OP_MULT));
        // PtrCoefFct partT = CoefFunction::Generate( mp_, part, CoefXprBinOp( mp_, "t", partTmp, CoefXpr::OP_MULT ) );

        eddyLossCoef->AddRegion(actRegion, partTmp);
        eddyLossCoefN->AddRegion(actRegion, partTmp);
      }
    }

    // =====================================================
    // MAG_FIELD_INTENSITY (START)
    // H = nu*B (linear case)
    // H = H(B) (nonlinear/hysteresis case)
    // =====================================================
    if (modelName_ == "invEBHysteresisModel")
    {
      shared_ptr<CoefFunctionMulti> hCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_FIELD_INTENSITY]);
      StdVector<RegionIdType>::iterator regIt = regions_.Begin(); // check which kind of nonlinearity is specified in the different regions
      regIt = regions_.Begin();
      for (; regIt != regions_.End(); ++regIt)
      {
        StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[*regIt]; // Just to find out which linear/nonlinear type is defined in this region
        PtrCoefFct h;                                                   // to store field intensity h
        if (nonLinTypes.Find(PERMEABILITY) != -1)
        { // hysteretic/nonlinear case
          h = nlFieldIntensitym_[*regIt]; //nonlinear_field_intensity_coef_;
          hCoef->AddRegion(*regIt, h);
        }
        else
        {
          h = CoefFunction::Generate(mp_, Global::REAL, CoefXprVecScalOp(mp_, GetCoefFct(MAG_FLUX_DENSITY), reluc_, CoefXpr::OP_MULT));
          hCoef->AddRegion(*regIt, h);
        }
      }
    }
    // =====================================================
    // MAG_FIELD_INTENSITY (END)
    // =====================================================
  }

  std::map<SolutionType, shared_ptr<FeSpace>>
  MagEdgePDE::CreateFeSpaces(const std::string &formulation,
                             PtrParamNode infoNode)
  {
    // ok default case so we create grid based approximation H1 elements
    // and standard Gauss integration
    std::map<SolutionType, shared_ptr<FeSpace>> crSpaces;
    if (formulation == "default" || formulation == "H_CURL")
    {
      PtrParamNode potSpaceNode = infoNode->Get("magPotential");
      crSpaces[MAG_POTENTIAL] =
          FeSpace::CreateInstance(myParam_, potSpaceNode, FeSpace::HCURL, ptGrid_);
      crSpaces[MAG_POTENTIAL]->Init(solStrat_);
    }
    else
    {
      EXCEPTION("The formulation " << formulation
                                   << "of magnetic edge PDE is not known!");
    }

    // in addition query, if special treatment of anisotropic elements
    // is activated
    if (myParam_->Has("thinElements"))
    {
      Double aspectRatio = 0.0;
      aspectRatio = myParam_->Get("thinElements")->Get("maxAspectRatio")->As<Double>();
      dynamic_pointer_cast<FeSpaceHCurlHi>(crSpaces[MAG_POTENTIAL])
          ->TreatThinElements(aspectRatio);
    }

    // ===================================================================
    // Create FeSpaceConst for coil currents (of coils excited by voltage)
    // ===================================================================
    if (hasVoltCoils_)
    {
      PtrParamNode voltSpaceNode = infoNode->Get("coilCurrent");
      crSpaces[COIL_CURRENT] =
          FeSpace::CreateInstance(myParam_, voltSpaceNode, FeSpace::CONSTANT, ptGrid_);
      crSpaces[COIL_CURRENT]->Init(solStrat_);
    }

    return crSpaces;
  }

} // end of namespace
