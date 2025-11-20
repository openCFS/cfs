#include <fstream>

#include "MagEdgeAdjTOPDE.hh"

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
#include "CoupledPDE/IterCoupledPDE.hh"

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
  DEFINE_LOG(magEdgeAdjTOPde, "magEdgeAdjTOPde")

  // **************
  //  Constructor
  // **************
  MagEdgeAdjTOPDE::MagEdgeAdjTOPDE(Grid *aptgrid, PtrParamNode paramNode,
                         PtrParamNode infoNode,
                         shared_ptr<SimState> simState, Domain *domain)
      : MagBasePDE(aptgrid, paramNode, infoNode, simState, domain)
  {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_ = "magneticEdgeAdjTO";
    pdematerialclass_ = ELECTROMAGNETIC;
    formulation_ = MagBasePDE::EDGE;

    //! Always use updated Lagrangian formulation
    updatedGeo_ = true; // true;

    // default is false
    useGradFields_ = paramNode->Get("useGradientFields")->As<bool>();
    onlyVacuum_ = paramNode->Get("onlyVacuum")->As<bool>();

    // check if we have a 3d setup
    //    bool is3d = domain_->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>() == "3d";
    if (!is3d_)
      EXCEPTION("MagEdgeAdjTOPDE is just implemented for 3D setups!");

    multiHarmCoef_.reset(new CoefFunctionHarmBalance<Complex>());

    // initialize material coef functions covering all regions :
    reluc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));
    perme_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));

    // Be aware that the hysteresis via CoefFunctionMaterialModel does not use the isHysteresis_ flag!
    matModelCoefm_.clear(); //matModelCoef_.reset(new CoefFunctionMaterialModel<Complex>());

    // init the nonlinear_field_intensity_coef_
    //nonlinear_field_intensity_coef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_, true));
  }

  // *************
  //  Destructor
  // *************
  MagEdgeAdjTOPDE::~MagEdgeAdjTOPDE()
  {
  }

  // *****************************
  //  Definition of Integrators
  // *****************************
  void MagEdgeAdjTOPDE::DefineIntegrators()
  {
    this->DefineStandardIntegrators();
    // in MagBasePDE
    DefineCoilIntegrators(1.0);
  }

  void MagEdgeAdjTOPDE::DefineStandardIntegrators()
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

    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_POTENTIAL_ADJ];
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

      // std::cout << "nonLinTypes.Find(PERMEABILITY): " << nonLinTypes.Find(PERMEABILITY) << "\n";
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
          PtrCoefFct curCoefDet = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_, curCoef, CoefXpr::OP_DET));
          reluc_->AddRegion(actRegion, curCoefDet);
        }
        else
        {

          // std::cout << " - Linear material case - " << std::endl;
          // ===============================
          //  Standard Stiffness Integrator
          // ===============================
          // curCoef = actMat->GetScalCoefFnc(MAG_RELUCTIVITY_SCALAR, Global::REAL);
          // actMat->GetTensorCoefFnc(MAG_RELUCTIVITY,FULL,Global::REAL );

          // ===============================
          //  Standard Stiffness Integrator
          // ===============================
          PtrCoefFct curCoef = NULL;

          CoefFunctionOpt *cfo = NULL; // we might do optimization and then we have such a thing

          // compute permeability
          PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");
          PtrCoefFct mu0 = CoefFunction::Generate(mp_, Global::REAL, "795774.716369");
          PtrCoefFct permeability = NULL;

          // TOAdjoint : read the necessary material parameter from h5 file !
          // std::cout << " - [TOAdjoint] : read the necessary material parameter from h5 file !" << std::endl;

          shared_ptr<ResultInfo> resultInfo = GetResultInfo(MAG_ELEM_PERMEABILITY); 
          shared_ptr<EntityList> entity = ptGrid_->GetEntityList( EntityList::ELEM_LIST, regionName );
        
          ReadMaterialDependency( "matPermeability", resultInfo->dofNames, resultInfo->entryType, false,
                                  entity, permeability, updatedGeo_ );
        
          perme_->AddRegion(actRegion, permeability);
          matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, permeability);
          
          PtrCoefFct reluctivity = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, constOne, 
                                                         permeability, CoefXpr::OP_DIV));
          reluc_->AddRegion(actRegion, reluctivity);  
          matCoefs_[MAG_ELEM_RELUCTIVITY]->AddRegion(actRegion, reluctivity);

          curlcurl = new BBInt<>(new CurlOperator<FeHCurl, 3, Double>(), reluctivity, 1.0, updatedGeo_);

          // when we have a CoefFunctionOpt, we tell it the proper form, which we only have now
          if (cfo)
            cfo->SetForm(curlcurl);

          //get the curl of the electric field of forward simulation
          std::string pdeName = "magneticEdge";
          BfieldForward_[actRegion] = iterCplPde_->GetCouplingCoefFct(MAG_FLUX_DENSITY, actSDList, pdeName, updatedGeo_);
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
        shared_ptr<BaseFeFunction> myFct = feFunctions_[MAG_POTENTIAL_ADJ];
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
          EXCEPTION("MagEdgeAdjTOPDE: VelocityStiffInt not implemented for complex case!");
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
        VelocityContextStiff->SetFeFunctions(feFunctions_[MAG_POTENTIAL_ADJ], feFunc);
        assemble_->AddBiLinearForm(VelocityContextStiff);
      } // end velocityId
    } // end for regions
  }

  LinearForm *MagEdgeAdjTOPDE::GetCurrentDensityInt(Double factor, PtrCoefFct coef, std::string coilId)
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
      EXCEPTION("MagEdgeAdjTOPDE is just implemented for 3D setups!");
    }
    return ret;
  }

  LinearForm *MagEdgeAdjTOPDE::GetRHSMagnetizationInt(Double factor, PtrCoefFct rhsMag, bool fullEvaluation)
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
      EXCEPTION("MagEdgeAdjTOPDE is just implemented for 3D setups!");
    }
    return lin;
  }

  BaseBDBInt *MagEdgeAdjTOPDE::GeHystStiffInt(Double factor, PtrCoefFct tensorReluctivity)
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
      EXCEPTION("MagEdgeAdjTOPDE is just implemented for 3D setups!");
    }

    return stiffInt;
  }

  void MagEdgeAdjTOPDE::DefineNcIntegrators()
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
          EXCEPTION("MagEdgeAdjTOPDE only works for 3D geometry!")
        else if (analysistype_ == MULTIHARMONIC)
        {
          DefineNitscheCoupling<3, 1>(MAG_POTENTIAL_ADJ, *ncIt, reluc_);
        }
        else
        {
          DefineNitscheCoupling<3, 1>(MAG_POTENTIAL_ADJ, *ncIt);
        }
        break;
      default:
        EXCEPTION("Unknown type of ncInterface");
        break;
      }
    }
  }

  void MagEdgeAdjTOPDE::DefineRhsLoadIntegrators()
  {
    LOG_TRACE(magEdgeAdjTOPde) << "Defining rhs load integrators for magEdgePDE";

    // Get FESpace and FeFunction of mechanical displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[MAG_POTENTIAL_ADJ];

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
            linearform_A_A = new BUIntegrator<Double>(new IdentityOperator<FeHCurl, 3, 1, Double>(), (0.0), GetCoefFct(MAG_POTENTIAL_ADJ), updatedGeo_);
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
            linearform_A_A = new BUIntegrator<Double>(new IdentityOperator<FeHCurl, 3, 1, Double>(), (-1.0), GetCoefFct(MAG_POTENTIAL_ADJ), updatedGeo_);
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
        // std::cout << "Hysteresis region found" << std::endl;
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
    //
    //    // =================================
    //    //  Magnetization -> from hysteresis (VOLUME) OLD
    //    // =================================
    //    //check for hysteresis
    //    if ( isHysteresis_ ){
    //      LOG_DBG(magEdgeAdjTOPde) << "Putting magnetization to rhs";
    //
    //      std::map<RegionIdType,PtrCoefFct > regionCoefs = magnetization_->GetRegionCoefs();
    //      std::map<RegionIdType,PtrCoefFct > regionHystCoefs = hysteresisCoefs_->GetRegionCoefs();
    //      std::map<RegionIdType, shared_ptr<CoefFunction> > ::iterator it;
    //      for( it = regionCoefs.begin(); it != regionCoefs.end(); it++) {
    //
    //        // get regionIdType
    //        RegionIdType curReg = it->first;
    //        PtrCoefFct curHystCoef = regionHystCoefs[curReg];
    //
    //        if(curHystCoef == NULL){
    //        //if(it->second == NULL){
    //          continue;
    //        }
    //        // get SDList
    //        shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
    //        actSDList->SetRegion( curReg );
    //
    //        // set fullevaluation to trigger evaluation at each integration point
    //        // the nonlinear parameter "evaluation depth" determines if each
    //        // integration point gets mapped to midpoint (> fullevaluation = false)
    //        // or if hyst operator really is evaluated at the actual int. point
    //        bool fullevaluation = true;
    //
    //        // NEW: we do not pass the hysteresis coefficient function
    //        // directly but instead a special class that returns the
    //        // correctly weighted term
    //        // even though, we have a similar function for output,
    //        // we need a separate coefFunction here as rhsMag might
    //        // be evaluated at another timestep/interation step as outputMag
    //        shared_ptr<CoefFunction> rhsMag = curHystCoef->GenerateRHSCoefFnc("MagMagnetization");
    //
    //        lin = GetRHSMagnetizationInt( magStrictFactor, rhsMag, fullevaluation );
    //
    //        lin->SetName("rhs_magnetization");
    //        lin->SetSolDependent();
    //        LinearFormContext *ctx = new LinearFormContext( lin );
    //        ctx->SetEntities( actSDList );
    //        ctx->SetFeFunction(myFct);
    //        assemble_->AddLinearForm(ctx);
    //        // Add entity list will add nothing, if entities were already assigned
    //        myFct->AddEntityList(actSDList);
    //      }
    //    }
    //    std::cout << "DONE Putting magnetization to rhs" << std::endl;
    // ==================
    //  FLUX DENSITY
    // ==================
    LOG_DBG(magEdgeAdjTOPde) << "Reading prescribed flux density";

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
    LOG_DBG(magEdgeAdjTOPde) << "Reading prescribed field intensity";

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

    // ==================
    //  MAGNETIC FLUX INTEGRATOR FOR ADJOINT RHS 
    // ==================

    LOG_DBG(magEdgeAdjTOPde) << "Reading precalculated magnetic flux";

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
        // Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
        // PtrCoefFct factor = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, reluc_, coef[i], CoefXpr::OP_MULT));
        PtrCoefFct factor = coef[i];

        if (isComplex_)
        {
          lin = new BUIntegrator<Complex>(new CurlOperator<FeHCurl, 3, Complex>(),
                                          Complex(-1.0), factor, coefUpdateGeo);
        }
        else
        {
          lin = new BUIntegrator<Double>(new CurlOperator<FeHCurl, 3, Double>(),
                                         -2.0, factor, coefUpdateGeo);
        }
      }

      lin->SetName("RHSFluxIntegrator");
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
    }
  }

  // ======================================================
  // TIME-STEPPING SECTION > MagBasePDE
  // ======================================================

  void MagEdgeAdjTOPDE::DefinePrimaryResults()
  {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    // === MAGNETIC VECTOR POTENTIAL ===
    shared_ptr<ResultInfo> potInfo(new ResultInfo);
    potInfo->resultType = MAG_POTENTIAL_ADJ;
    potInfo->dofNames = vecComponents;
    potInfo->unit = "Vs/m";
    potInfo->definedOn = ResultInfo::ELEMENT;
    potInfo->entryType = ResultInfo::VECTOR;

    feFunctions_[MAG_POTENTIAL_ADJ]->SetResultInfo(potInfo);
    DefineFieldResult(feFunctions_[MAG_POTENTIAL_ADJ], potInfo);

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MAG_POTENTIAL_ADJ] = "fluxParallel";
    idbcSolNameMap_[MAG_POTENTIAL_ADJ] = "potential";

    // === PERMEABILITY ===
    shared_ptr<ResultInfo> permeability(new ResultInfo);
    permeability->resultType = MAG_ELEM_PERMEABILITY;
    permeability->dofNames = "";
    permeability->unit = "Vs/Am";
    permeability->definedOn = ResultInfo::ELEMENT;
    permeability->entryType = ResultInfo::SCALAR;

    // In multiharmonic analysis we have complex permeability
    shared_ptr<CoefFunctionMulti> permFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, false));
    matCoefs_[MAG_ELEM_PERMEABILITY] = permFct;
    DefineFieldResult(permFct, permeability);

    // === RELUCTIVITY ===
    shared_ptr<ResultInfo> reluctivity(new ResultInfo);
    reluctivity->resultType = MAG_ELEM_RELUCTIVITY;
    reluctivity->dofNames = "";
    reluctivity->unit = "Vs/Am";
    reluctivity->definedOn = ResultInfo::ELEMENT;
    reluctivity->entryType = ResultInfo::SCALAR;

    // In multiharmonic analysis we have complex permeability
    shared_ptr<CoefFunctionMulti> relucFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, false));
    matCoefs_[MAG_ELEM_RELUCTIVITY] = relucFct;
    DefineFieldResult(relucFct, reluctivity);

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
  }

  void MagEdgeAdjTOPDE::DefinePostProcResults()
  {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL_ADJ];

    if (!fluxDensityDefined_)
    {
      DefineMagFluxDensity();
    }

    PtrCoefFct bFunc = this->GetCoefFct(MAG_FLUX_DENSITY_ADJ);

    // === MAGNETIC RHS ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = MAG_RHS_LOAD_ADJ;
    rhs->dofNames = vecComponents;
    rhs->unit = "-";
    rhs->entryType = ResultInfo::VECTOR;
    rhs->definedOn = ResultInfo::ELEMENT;
    rhsFeFunctions_[MAG_POTENTIAL_ADJ]->SetResultInfo(rhs);
    DefineFieldResult(rhsFeFunctions_[MAG_POTENTIAL_ADJ], rhs);

    
  }

  void MagEdgeAdjTOPDE::FinalizePostProcResults()
  {

    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();

    StdVector<RegionIdType>::iterator regIt = regions_.Begin();

    shared_ptr<CoefFunctionMulti> magFluxDensCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_FLUX_DENSITY_ADJ]);

    // ======================== Computation of gardients for SIMP =================
    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL_ADJ];

    // === Gradient via adjoint (element-by-element for topology optimization) ===
    // === Parameter 1 is permeability
    shared_ptr<ResultInfo> gradAdjParam1;
    gradAdjParam1.reset(new ResultInfo);
    gradAdjParam1->resultType = GRAD_PARAM1;
    gradAdjParam1->dofNames = "";
    gradAdjParam1->unit = "";
    gradAdjParam1->entryType = ResultInfo::SCALAR;
    gradAdjParam1->definedOn = ResultInfo::ELEMENT;
    shared_ptr<CoefFunctionMulti> optMult1(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
    DefineFieldResult(optMult1, gradAdjParam1);

    shared_ptr<CoefFunctionMulti> scalMult1 = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[GRAD_PARAM1]);
    shared_ptr<CoefFunctionMulti> bField = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_FLUX_DENSITY_ADJ]);

    regIt = regions_.Begin();
      for (; regIt != regions_.End(); ++regIt) {
        RegionIdType actRegion = *regIt;
        PtrCoefFct bForward = NULL;
        PtrCoefFct mult1 = NULL;
        if(BfieldForward_.find(*regIt) != BfieldForward_.end()){
          // Magnetic flux density = MAG_FLUX_DENSITY from forward simulation
          bForward = BfieldForward_[*regIt];
          // get coef of MAG_FLUX_DENSITY_ADJ: is adjoint solution!
          PtrCoefFct bFieldAdj = this->GetCoefFct(MAG_FLUX_DENSITY_ADJ);
          PtrCoefFct mu = perme_->GetRegionCoef(actRegion);
          PtrCoefFct nu = reluc_->GetRegionCoef(actRegion);
          PtrCoefFct nu_squared = CoefFunction::Generate( mp_, Global::REAL,
                                        CoefXprBinOp( mp_, nu, nu, CoefXpr::OP_MULT ) );
          
          PtrCoefFct bProduct = CoefFunction::Generate( mp_, Global::REAL,
                                        CoefXprBinOp( mp_, bForward, bFieldAdj, CoefXpr::OP_MULT) );
          mult1 = CoefFunction::Generate( mp_, Global::REAL,
                             CoefXprBinOp(mp_, nu_squared, bProduct, CoefXpr::OP_MULT) );
          scalMult1->AddRegion(actRegion, mult1);

          //curl of B-field from forward simulation
          // bField->AddRegion(*regIt, bForward);
        }
      } // loop over regions
  }

  std::map<SolutionType, shared_ptr<FeSpace>>
  MagEdgeAdjTOPDE::CreateFeSpaces(const std::string &formulation,
                             PtrParamNode infoNode)
  {
    // ok default case so we create grid based approximation H1 elements
    // and standard Gauss integration
    std::map<SolutionType, shared_ptr<FeSpace>> crSpaces;
    if (formulation == "default" || formulation == "H_CURL")
    {
      PtrParamNode potSpaceNode = infoNode->Get("magPotentialAdj");
      crSpaces[MAG_POTENTIAL_ADJ] =
          FeSpace::CreateInstance(myParam_, potSpaceNode, FeSpace::HCURL, ptGrid_);
      crSpaces[MAG_POTENTIAL_ADJ]->Init(solStrat_);
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
      dynamic_pointer_cast<FeSpaceHCurlHi>(crSpaces[MAG_POTENTIAL_ADJ])
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
