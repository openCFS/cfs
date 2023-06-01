#include <fstream>

#include "FullWaveMaxwellAVPDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/Coil.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"

#include "Driver/Assemble.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "FeBasis/HCurl/FeSpaceHCurlHi.hh"
#include "FeBasis/HCurl/HCurlElems.hh"
#include "FeBasis/H1/FeSpaceH1Hi.hh"
#include "FeBasis/H1/H1Elems.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

#include "Domain/CoefFunction/CoefFunctionExpression.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefXpr.hh"

// forms
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/LinForms/KXInt.hh"
#include "Forms/Operators/CurlOperator.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/DivOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"

// time stepping
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
namespace CoupledField
{

  // declare class specific logging stream
  DEFINE_LOG(FullWaveMaxwellAVPDE, "FullWaveMaxwellAVPDE")

  // **************
  //  Constructor
  // **************
  FullWaveMaxwellAVPDE::FullWaveMaxwellAVPDE(Grid *aptgrid, PtrParamNode paramNode,
                                         PtrParamNode infoNode,
                                         shared_ptr<SimState> simState, Domain *domain)
      : SinglePDE(aptgrid, paramNode, infoNode, simState, domain)
  {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_ = "fullwave-AV";
    // We reuse the Darwin material class (same as classical magnetic material plus permittivity)
    pdematerialclass_ = ELECTROMAGNETIC_DARWIN;

    //! Always use updated Lagrangian formulation
    updatedGeo_ = true; // true;

    formulationType_ = paramNode->Get("formulation")->As<std::string>();

    // check if we have a 3d setup
    bool is3d = domain_->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>() == "3d";
    if (!is3d)
      EXCEPTION("FullWaveMaxwellAVPDE is just implemented for 3D setups!");

    // initialize material coef functions covering all regions
    reluc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));
    conduc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
    eps_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
  }

  // *************
  //  Destructor
  // *************
  FullWaveMaxwellAVPDE::~FullWaveMaxwellAVPDE()
  {
  }

  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void FullWaveMaxwellAVPDE::InitNonLin()
  {

    SinglePDE::InitNonLin();
  }

  // *****************************
  //  Definition of Integrators
  // *****************************
  void FullWaveMaxwellAVPDE::DefineIntegrators()
  {

    RegionIdType actRegion;
    BaseMaterial *actMat = NULL;

    // initially, check for regularization factor
    Double regularizationFactor = 1e-6;
    myParam_->GetValue("penaltyFactor", regularizationFactor, ParamNode::PASS);

    shared_ptr<BaseFeFunction> magVecPotFeFunc = feFunctions_[MAG_POTENTIAL];
    shared_ptr<BaseFeFunction> elecScalPotFeFunc = feFunctions_[ELEC_POTENTIAL];
    shared_ptr<FeSpace> magVecPotFeSpace = magVecPotFeFunc->GetFeSpace();
    shared_ptr<FeSpace> elecScalPotFeSpace = elecScalPotFeFunc->GetFeSpace();

    PtrCoefFct magFluxCoef = this->GetCoefFct(MAG_FLUX_DENSITY);

    for (UInt iRegion = 0; iRegion < regions_.GetSize(); iRegion++)
    {
      actRegion = regions_[iRegion];
      actMat = materials_[actRegion];
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region", "name", regionName.c_str());

      // Get flag if we need to consider an electric scalar potential in this region
      bool isConducRegion = curRegNode->Get("isConducRegion")->As<bool>();

      // Get polynomial and integration order for magnetic vector potential
      std::string magVecPolyId = curRegNode->Get("magVecPolyId")->As<std::string>();
      std::string magVecIntegId = curRegNode->Get("magVecIntegId")->As<std::string>();
      magVecPotFeSpace->SetRegionApproximation(actRegion, magVecPolyId, magVecIntegId);

      // Get polynomial and integration order for electric scalar potential
      std::string elecScalPolyId = curRegNode->Get("elecScalPolyId")->As<std::string>();
      std::string elecScalIntegId = curRegNode->Get("elecScalIntegId")->As<std::string>();
      elecScalPotFeSpace->SetRegionApproximation(actRegion, elecScalPolyId, elecScalIntegId);

      // Get possible nonlinearities defined in this region
      StdVector<NonLinType> matDepenTypes = regionMatDepTypes_[actRegion]; // material dependency
      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];   // non-linearity

      // Create new entity list, based on the region
      shared_ptr<ElemList> actSDList(new ElemList(ptGrid_));
      actSDList->SetRegion(actRegion);

      // Pass entitylist to fespace / fefunction for magnetic vector and scalar potential
      magVecPotFeFunc->AddEntityList(actSDList);
      elecScalPotFeFunc->AddEntityList(actSDList);

      if (matDepenTypes.Find(NLELEC_CONDUCTIVITY) != -1)
      {
        EXCEPTION("FullWaveMaxwellAVPDE does not support nonlinear"
                  "(temperatur dependent) electric conductivity yet!\n"
                  "But the implementation is not big deal...I promise.\n"
                  "Mostly c&p form MagEdgePDE.");
      }

      if (nonLinTypes.GetSize() > 0)
      {
        // =================================================================================
        //  NONLINEAR SECTION
        // =================================================================================
        EXCEPTION("FullWaveMaxwellAVPDE does not support nonlinear reluctivity yet!";)
      }
      else
      {
        /* ==============================================
         * Handling of material parameters
           ============================================== */
        // Magnetic Reluctivity
        PtrCoefFct nuNl = NULL;
        nuNl = actMat->GetScalCoefFnc(MAG_RELUCTIVITY_SCALAR, Global::REAL);
        // Add material to global, distributed reluctivity coefficient function
        reluc_->AddRegion(actRegion, nuNl);

        // Magnetic Permeability
        PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");
        PtrCoefFct permeability = CoefFunction::Generate(mp_, Global::REAL, CoefXprBinOp(mp_, constOne, nuNl, CoefXpr::OP_DIV));
        matCoefs_[MAG_ELEM_PERMEABILITY]->AddRegion(actRegion, permeability);

        // Electric Conductivity
        Double conductivity;
        materials_[actRegion]->GetScalar(conductivity, MAG_CONDUCTIVITY_SCALAR, Global::REAL);
        PtrCoefFct conducCoef = materials_[actRegion]->GetScalCoefFnc(MAG_CONDUCTIVITY_SCALAR, Global::REAL);

        // Electric Permittivity
        PtrCoefFct eps = actMat->GetScalCoefFnc(MAG_PERMITTIVITY_SCALAR, Global::REAL);
        // Add material to global, distributed coefficient function
        eps_->AddRegion(actRegion, eps);


        // =================================================================================
        // =================================================================================
        //  STIFFNESS SECTION
        // =================================================================================
        // =================================================================================

        /* ==============================================
         * K_A_A^(nu):
         * nu curl(A) \cdot curl(A)
           ============================================== */
        BaseBDBInt *K_A_A_nu = NULL;
        K_A_A_nu = new BBInt<>(new CurlOperator<FeHCurl, 3, Double>(), nuNl, 1.0, updatedGeo_);
        K_A_A_nu->SetName("K_A_A_nu");

        BiLinFormContext *K_A_A_nuContext = new BiLinFormContext(K_A_A_nu, STIFFNESS);
        K_A_A_nuContext->SetEntities(actSDList, actSDList);
        K_A_A_nuContext->SetFeFunctions(magVecPotFeFunc, magVecPotFeFunc);
        assemble_->AddBiLinearForm(K_A_A_nuContext);
        // Add bdb-integrator to global list, needed for flux density evaluation
        bdbInts_.insert(std::pair<RegionIdType, BaseBDBInt *>(actRegion, K_A_A_nu));

        if (isConducRegion)
        {
          /* ==============================================
           * K_A_phi^(sigma):
           * \sigma A \cdot grad(V)
         ============================================== */
          BiLinearForm *K_A_phi_sigma = NULL;
          K_A_phi_sigma = new ABInt<>(new IdentityOperator<FeHCurl, 3, 1, Double>(), new GradientOperator<FeH1, 3, 1, Double>(),
                                      conducCoef, 1.0, updatedGeo_);
          K_A_phi_sigma->SetName("K_A_phi_sigma");

          BiLinFormContext *K_A_phi_sigmaContext = new BiLinFormContext(K_A_phi_sigma, STIFFNESS);
          K_A_phi_sigmaContext->SetEntities(actSDList, actSDList);
          K_A_phi_sigmaContext->SetFeFunctions(magVecPotFeFunc, elecScalPotFeFunc);
          assemble_->AddBiLinearForm(K_A_phi_sigmaContext);
        }
        
        /* ==============================================
          * K_phi_A^(epsilon):
          * \epsilon grad(V) \cdot  A
        ============================================== */
        BiLinearForm *K_phi_A_epsilon = NULL;
        K_phi_A_epsilon = new ABInt<>(new GradientOperator<FeH1, 3, 1, Double>(), new IdentityOperator<FeHCurl, 3, 1, Double>(),
                                    eps, 1.0, updatedGeo_);
        K_phi_A_epsilon->SetName("K_phi_A_epsilon");

        BiLinFormContext *K_phi_A_epsilonContext = new BiLinFormContext(K_phi_A_epsilon, STIFFNESS);
        K_phi_A_epsilonContext->SetEntities(actSDList, actSDList);
        K_phi_A_epsilonContext->SetFeFunctions(elecScalPotFeFunc, magVecPotFeFunc);
        assemble_->AddBiLinearForm(K_phi_A_epsilonContext);
        
        // =================================================================================
        // =================================================================================
        //  DAMPING SECTION
        // =================================================================================
        // =================================================================================

        /* ==============================================
          * D_A_A^(sigma):
          * \sigma (A) \cdot (A)
          ============================================== */
        if (isConducRegion)
        {
          BaseBDBInt *D_A_A_sigma;
          //        	D_A_A_sigma = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<FeHCurl,3,Double>(), conducCoef,1.0, updatedGeo_);
          D_A_A_sigma = new BBIntMassEdge<>(new IdentityOperator<FeHCurl, 3, 1, Double>(), conducCoef, 1.0, updatedGeo_);
          D_A_A_sigma->SetName("K_A_A_sigma");

          BiLinFormContext *D_A_A_sigmaContext;
          D_A_A_sigmaContext = new BiLinFormContext(D_A_A_sigma, DAMPING);
          D_A_A_sigmaContext->SetEntities(actSDList, actSDList);
          D_A_A_sigmaContext->SetFeFunctions(magVecPotFeFunc, magVecPotFeFunc);
          assemble_->AddBiLinearForm(D_A_A_sigmaContext);
        }

        /* ==============================================
          * D_A_phi^(epsilon):
          * \epsilon (A) \cdot grad(phi)
          ============================================== */
        BaseBDBInt *D_A_phi_epsilon;
        //        	D_A_phi_epsilon = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<FeHCurl,3,Double>(), conducCoef,1.0, updatedGeo_);
        D_A_phi_epsilon = new ABInt<>(new IdentityOperator<FeHCurl, 3, 1, Double>(), new GradientOperator<FeH1, 3, 1, Double>(),
                                      eps, 1.0, updatedGeo_);
        D_A_phi_epsilon->SetName("D_A_phi_sigma");

        BiLinFormContext *D_A_phi_epsilonContext;
        D_A_phi_epsilonContext = new BiLinFormContext(D_A_phi_epsilon, DAMPING);
        D_A_phi_epsilonContext->SetEntities(actSDList, actSDList);
        D_A_phi_epsilonContext->SetFeFunctions(magVecPotFeFunc, elecScalPotFeFunc);
        assemble_->AddBiLinearForm(D_A_phi_epsilonContext);
      

        // =================================================================================
        // =================================================================================
        //  MASS SECTION
        // =================================================================================
        // =================================================================================
        BaseBDBInt *M_A_A_epsilon;
        //        	M_A_A_epsilon = new BBIntMassEdge<>(new ScaledByEdgeIdentityOperator<FeHCurl,3,Double>(), eps,1.0, updatedGeo_);
        M_A_A_epsilon = new BBIntMassEdge<>(new IdentityOperator<FeHCurl, 3, 1, Double>(), eps, 1.0, updatedGeo_);
        M_A_A_epsilon->SetName("K_A_A_sigma");

        BiLinFormContext *M_A_A_epsilonContext;
        M_A_A_epsilonContext = new BiLinFormContext(M_A_A_epsilon, MASS);
        M_A_A_epsilonContext->SetEntities(actSDList, actSDList);
        M_A_A_epsilonContext->SetFeFunctions(magVecPotFeFunc, magVecPotFeFunc);
        assemble_->AddBiLinearForm(M_A_A_epsilonContext);
      } // END OF NONLIN/LIN PART
    }   // end for regions
  }     // end DefineIntegrators

  void FullWaveMaxwellAVPDE::DefineNcIntegrators()
  {
    StdVector<NcInterfaceInfo>::iterator ncIt = ncInterfaces_.Begin(), endIt = ncInterfaces_.End();
    for (; ncIt != endIt; ++ncIt)
    {
      switch (ncIt->type)
      {
      case NC_MORTAR:
        EXCEPTION("FullWaveMaxwellAVPDE: Mortar NC interfaces not tested!");
      case NC_NITSCHE:
      {
        /*
         * that's kind of a dirty hack because for Nitsche NC, we need to access the
         * electric conductivity as MAG_CONDUCTIVITY_SCALAR. But this should only be done in
         * the FullWaveMaxwellAVPDE
         */
        shared_ptr<CoefFunctionMulti> identifier = NULL;
        identifier.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, true));
        if (dim_ == 2)
          EXCEPTION("FullWaveMaxwellAVPDE possible only in 3D!")
        else
          // DefineNitscheCoupling<3,1>(ELEC_POTENTIAL, *ncIt, identifier );
          DefineNitscheCoupling<3, 1>(MAG_POTENTIAL, *ncIt);
        break;
      }
      default:
        EXCEPTION("Unknown type of ncInterface");
        break;
      }
    }
  }

  void FullWaveMaxwellAVPDE::DefineSurfaceIntegrators()
  {
  }

  void FullWaveMaxwellAVPDE::DefineRhsLoadIntegrators()
  {
    //LOG_TRACE(FullWaveMaxwellAVPDE) << "Defining rhs load integrators for FullWaveMaxwellAVPDE";

    shared_ptr<BaseFeFunction> magVecPotFeFunc = feFunctions_[MAG_POTENTIAL];
    shared_ptr<BaseFeFunction> elecScalPotFeFunc = feFunctions_[ELEC_POTENTIAL];

    StdVector<shared_ptr<EntityList>> ent;
    StdVector<PtrCoefFct> coef;

    StdVector<std::string> vecDofNames = magVecPotFeFunc->GetResultInfo()->dofNames;
    StdVector<std::string> scalDofNames = elecScalPotFeFunc->GetResultInfo()->dofNames;

    bool coefUpdateGeo = true;
    // ==================
    //  FLUX DENSITY
    // ==================
    //LOG_DBG(FullWaveMaxwellAVPDE) << "Reading prescribed flux density";

    ReadRhsExcitation("fluxDensity", vecDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo);
    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
      EXCEPTION("Currently no rhs for fluxDensity possible in FullWaveMaxwellAVPDE");
    }

    // ==================
    //  FIELD INTENSITY
    // ==================
    //LOG_DBG(FullWaveMaxwellAVPDE) << "Reading prescribed field intensity";

    ReadRhsExcitation("fieldIntensity", vecDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo);
    for (UInt i = 0; i < ent.GetSize(); ++i)
    {
      EXCEPTION("Currently no rhs for fieldIntensity possible in FullWaveMaxwellAVPDE");
    }
  }

  void FullWaveMaxwellAVPDE::DefineSolveStep()
  {
    solveStep_ = new StdSolveStep(*this);
  }


  // ======================================================
  // TIME-STEPPING SECTION
  // ======================================================

  void FullWaveMaxwellAVPDE::InitTimeStepping()
  {
    // Use complete implicit scheme
    Double gamma = 1.0;
    GLMScheme *scheme = new Trapezoidal(gamma);
    TimeSchemeGLM::NonLinType nlType = (nonLin_) ? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType));
    feFunctions_[MAG_POTENTIAL]->SetTimeScheme(myScheme);

    // Important: Create a new time scheme just for the elec potential unknowns, as otherwise the
    // size of the vectors does not match!
    GLMScheme *scheme2 = new Trapezoidal(gamma);
    shared_ptr<BaseTimeScheme> myScheme2(new TimeSchemeGLM(scheme2, 0, nlType));
    feFunctions_[ELEC_POTENTIAL]->SetTimeScheme(myScheme2);
  }

  // ******************************************************
  //   Query parameter object for information about coils
  // ******************************************************
  void FullWaveMaxwellAVPDE::ReadCoils()
  {
    PtrParamNode coilNode = myParam_->Get("coilList", ParamNode::PASS);
    PtrParamNode coilInfoNode = myInfo_->Get("coilList", ParamNode::PASS);
    if (coilNode)
    {
      EXCEPTION("Currently no coils are supported for FullWaveMaxwellAVPDE");
    }
  }

  void FullWaveMaxwellAVPDE::DefinePrimaryResults()
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

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MAG_POTENTIAL] = "fluxParallel";
    idbcSolNameMap_[MAG_POTENTIAL] = "potential";

    // === ELECTRIC SCALAR POTENTIAL ===
    shared_ptr<BaseFeFunction> scalFct = feFunctions_[ELEC_POTENTIAL];
    shared_ptr<ResultInfo> res2(new ResultInfo);
    res2->resultType = ELEC_POTENTIAL;
    res2->dofNames = "";
    res2->unit = "V";
    res2->definedOn = ResultInfo::NODE;
    res2->entryType = ResultInfo::SCALAR;
    results_.Push_back(res2);
    availResults_.insert(res2);
    scalFct->SetResultInfo(res2);
    DefineFieldResult(scalFct, res2);

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    // hdbcSolNameMap_[ELEC_POTENTIAL] = "fluxParallel";
    idbcSolNameMap_[ELEC_POTENTIAL] = "elecPotential";

    // === PERMEABILITY ===
    shared_ptr<ResultInfo> permeability(new ResultInfo);
    permeability->resultType = MAG_ELEM_PERMEABILITY;
    permeability->dofNames = "";
    permeability->unit = "Vs/Am";
    permeability->definedOn = ResultInfo::ELEMENT;
    permeability->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> permFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, false));
    matCoefs_[MAG_ELEM_PERMEABILITY] = permFct;
    DefineFieldResult(permFct, permeability);
  }

  void FullWaveMaxwellAVPDE::DefinePostProcResults()
  {
    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    shared_ptr<BaseFeFunction> magVecPotFeFct = feFunctions_[MAG_POTENTIAL];
    shared_ptr<BaseFeFunction> elecScalPotFeFct = feFunctions_[ELEC_POTENTIAL];

    // === TIME DERIVATIVES OF PRIMARY RESULTS ===
    if (analysistype_ != STATIC)
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

      // === GRADIENT ELEC SCALAR POTENTIAL ===
      shared_ptr<ResultInfo> gradV(new ResultInfo);
      gradV->resultType = GRAD_ELEC_POTENTIAL;
      gradV->dofNames = vecComponents;
      gradV->unit = "V/m";
      gradV->definedOn = ResultInfo::ELEMENT;
      gradV->entryType = ResultInfo::VECTOR;
      availResults_.insert(gradV);
      shared_ptr<CoefFunctionFormBased> gradVFunc;
      if (isComplex_)
      {
        gradVFunc.reset(new CoefFunctionBOp<Complex>(elecScalPotFeFct, gradV));
      }
      else
      {
        gradVFunc.reset(new CoefFunctionBOp<Double>(elecScalPotFeFct, gradV));
      }
      DefineFieldResult(gradVFunc, gradV);
      stiffFormCoefsAux1_.insert(gradVFunc);

      // === ELECTRIC FIELD INTENSITY TRANSVERSAL ===
      shared_ptr<ResultInfo> elecIntensT(new ResultInfo);
      elecIntensT->resultType = ELEC_FIELD_INTENSITY_TRANSVERSAL;
      elecIntensT->SetVectorDOFs(dim_, isaxi_);
      elecIntensT->dofNames = vecComponents;
      elecIntensT->unit = "V/m";
      elecIntensT->definedOn = ResultInfo::ELEMENT;
      elecIntensT->entryType = ResultInfo::VECTOR;
      shared_ptr<CoefFunctionMulti> elecIntensTFunc(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
      DefineFieldResult(elecIntensTFunc, elecIntensT);

      // === ELECTRIC FIELD INTENSITY LONGITUDINAL ===
      shared_ptr<ResultInfo> elecIntensL(new ResultInfo);
      elecIntensL->resultType = ELEC_FIELD_INTENSITY_LONGITUDINAL;
      elecIntensL->SetVectorDOFs(dim_, isaxi_);
      elecIntensL->dofNames = vecComponents;
      elecIntensL->unit = "V/m";
      elecIntensL->definedOn = ResultInfo::ELEMENT;
      elecIntensL->entryType = ResultInfo::VECTOR;
      shared_ptr<CoefFunctionMulti> elecIntensLFunc(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
      DefineFieldResult(elecIntensLFunc, elecIntensL);

      // === ELECTRIC FIELD INTENSITY ===
      shared_ptr<ResultInfo> elecIntens(new ResultInfo);
      elecIntens->resultType = ELEC_FIELD_INTENSITY;
      elecIntens->SetVectorDOFs(dim_, isaxi_);
      elecIntens->dofNames = vecComponents;
      elecIntens->unit = "V/m";
      elecIntens->definedOn = ResultInfo::ELEMENT;
      elecIntens->entryType = ResultInfo::VECTOR;
      shared_ptr<CoefFunctionMulti> elecIntensFunc(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
      DefineFieldResult(elecIntensFunc, elecIntens);

      // === DISPLACEMENT CURRENT DENSITY ===
      shared_ptr<ResultInfo> displcurrIntens(new ResultInfo);
      displcurrIntens->resultType = DISPLACEMENT_CURRENT_FIELD_INTENSITY;
      displcurrIntens->SetVectorDOFs(dim_, isaxi_);
      displcurrIntens->dofNames = vecComponents;
      displcurrIntens->unit = "";
      displcurrIntens->definedOn = ResultInfo::ELEMENT;
      displcurrIntens->entryType = ResultInfo::VECTOR;
      shared_ptr<CoefFunctionMulti> displcurrIntensFunc(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
      DefineFieldResult(displcurrIntensFunc, displcurrIntens);

      // === DISPLACEMENT CURRENT ===
      shared_ptr<ResultInfo> displcurr(new ResultInfo);
      displcurr->resultType = DISPLACEMENT_CURRENT_SURF;
      displcurr->SetVectorDOFs(dim_, isaxi_);
      displcurr->dofNames = "";
      displcurr->unit = "";
      displcurr->definedOn = ResultInfo::SURF_REGION;
      displcurr->entryType = ResultInfo::SCALAR;
      availResults_.insert(displcurr);
      // first, create normal mapping
      shared_ptr<CoefFunctionSurf> displcurrSurf(new CoefFunctionSurf(true, 1.0, displcurr));
      surfCoefFcts_[displcurrSurf] = displcurrIntensFunc;
      // then, integrate values
      shared_ptr<ResultFunctor> displcurrSurfFunc;
      if (isComplex_)
      {
        displcurrSurfFunc.reset(new ResultFunctorIntegrate<Complex>(displcurrSurf, magVecPotFeFct, displcurr));
      }
      else
      {
        displcurrSurfFunc.reset(new ResultFunctorIntegrate<Double>(displcurrSurf, magVecPotFeFct, displcurr));
      }
      resultFunctors_[DISPLACEMENT_CURRENT_SURF] = displcurrSurfFunc;
    }

    // === EDDY CURRENT DENSITY ===
    shared_ptr<ResultInfo> eddyJ(new ResultInfo);
    eddyJ->resultType = MAG_EDDY_CURRENT_DENSITY;
    eddyJ->dofNames = vecComponents;
    eddyJ->unit = "A/m^2";
    eddyJ->definedOn = ResultInfo::ELEMENT;
    eddyJ->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionMulti> eddyJFunc(
        new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_));
    DefineFieldResult(eddyJFunc, eddyJ);

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
    surfCoefFcts_[ncd] = eddyJFunc;
    // then, integrate values
    shared_ptr<ResultFunctor> eddyCurrentFunc;
    if (isComplex_)
    {
      eddyCurrentFunc.reset(new ResultFunctorIntegrate<Complex>(ncd, magVecPotFeFct, ec));
    }
    else
    {
      eddyCurrentFunc.reset(new ResultFunctorIntegrate<Double>(ncd, magVecPotFeFct, ec));
    }
    resultFunctors_[MAG_EDDY_CURRENT] = eddyCurrentFunc;

    // === MAGNETIC FLUX DENSITY ===
    shared_ptr<ResultInfo> fluxDens(new ResultInfo);
    fluxDens->resultType = MAG_FLUX_DENSITY;
    fluxDens->dofNames = vecComponents;
    fluxDens->unit = "Vs/m^2";
    fluxDens->definedOn = ResultInfo::ELEMENT;
    fluxDens->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionFormBased> bFunc;
    if (isComplex_)
    {
      bFunc.reset(new CoefFunctionBOp<Complex>(magVecPotFeFct, fluxDens));
    }
    else
    {
      bFunc.reset(new CoefFunctionBOp<Double>(magVecPotFeFct, fluxDens));
    }
    DefineFieldResult(bFunc, fluxDens);
    stiffFormCoefs_.insert(bFunc);

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
      fluxFct.reset(new ResultFunctorIntegrate<Complex>(sNormFDens, magVecPotFeFct, flux));
    }
    else
    {
      fluxFct.reset(new ResultFunctorIntegrate<Double>(sNormFDens, magVecPotFeFct, flux));
    }
    resultFunctors_[MAG_FLUX] = fluxFct;
    availResults_.insert(flux);
  }

  void FullWaveMaxwellAVPDE::FinalizePostProcResults()
  {
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;

    shared_ptr<BaseFeFunction> magVecPotFeFct = feFunctions_[MAG_POTENTIAL];
    shared_ptr<BaseFeFunction> elecScalPotFeFct = feFunctions_[ELEC_POTENTIAL];

    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();

    // === ELECTRIC FIELD INTENSITY ===
    // Assemble coefficient function for
    // E_T = - \nabla V
    // E_L = -\frac{\partial A}{\partial t}
    // E = E_T + E_L

    // === EDDY CURRENT DENSITY ===
    // Assemble coefficient function for
    // J = \gamma (E_T + E_L)

    shared_ptr<CoefFunctionMulti> elecIntensTCoef =
        dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[ELEC_FIELD_INTENSITY_TRANSVERSAL]);
    shared_ptr<CoefFunctionMulti> elecIntensLCoef =
        dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[ELEC_FIELD_INTENSITY_LONGITUDINAL]);
    shared_ptr<CoefFunctionMulti> elecIntensCoef =
        dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[ELEC_FIELD_INTENSITY]);
    shared_ptr<CoefFunctionMulti> jDensLCoef =
        dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_EDDY_CURRENT_DENSITY]);
    shared_ptr<CoefFunctionMulti> dispCurrDensCoef =
        dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[DISPLACEMENT_CURRENT_FIELD_INTENSITY]);

    StdVector<RegionIdType>::iterator regIt = regions_.Begin();
    PtrCoefFct constOne = CoefFunction::Generate(mp_, Global::REAL, "-1.0");
    regIt = regions_.Begin();
    for (; regIt != regions_.End(); ++regIt)
    {
      std::string regionName = ptGrid_->GetRegion().ToString(*regIt);
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region", "name", regionName.c_str());

      // Material parameters
      Double conductivity, permittivity;
      materials_[*regIt]->GetScalar(conductivity, MAG_CONDUCTIVITY_SCALAR, Global::REAL);
      PtrCoefFct conducCoef = materials_[*regIt]->GetScalCoefFnc(MAG_CONDUCTIVITY_SCALAR, Global::REAL);
      materials_[*regIt]->GetScalar(permittivity, MAG_PERMITTIVITY_SCALAR, Global::REAL);
      PtrCoefFct permittivityCoef = materials_[*regIt]->GetScalCoefFnc(MAG_PERMITTIVITY_SCALAR, Global::REAL);

      // Electric Field
      PtrCoefFct longitudinal = CoefFunction::Generate(mp_,
                                                       part, CoefXprVecScalOp(mp_, GetCoefFct(GRAD_ELEC_POTENTIAL), constOne, CoefXpr::OP_MULT));
      PtrCoefFct transversal = CoefFunction::Generate(mp_,
                                                      part, CoefXprVecScalOp(mp_, GetCoefFct(MAG_POTENTIAL_DERIV1), constOne, CoefXpr::OP_MULT));
      PtrCoefFct sum = CoefFunction::Generate(mp_,
                                              part, CoefXprBinOp(mp_, transversal, longitudinal, CoefXpr::OP_ADD));

      elecIntensTCoef->AddRegion(*regIt, transversal);
      elecIntensLCoef->AddRegion(*regIt, longitudinal);
      elecIntensCoef->AddRegion(*regIt, sum);

      // Current Density
      PtrCoefFct curr_tmp = CoefFunction::Generate(mp_,
                                                   part, CoefXprBinOp(mp_, longitudinal, transversal, CoefXpr::OP_ADD));
      PtrCoefFct curr = CoefFunction::Generate(mp_,
                                               part, CoefXprVecScalOp(mp_, curr_tmp, conducCoef, CoefXpr::OP_MULT));

      jDensLCoef->AddRegion(*regIt, curr);

      // Displacement Current Density dD/dt = epsilon * dE/dt
      // and in frequency-domain hat(D) = j * 2*pi*f * hat(E)
      // but remember: due to the Darwin approximation, dE/dt only consists of the time-derivative of the scalar potential
      PtrCoefFct tmp = CoefFunction::Generate(mp_,
                                              part, CoefXprVecScalOp(mp_, longitudinal, permittivityCoef, CoefXpr::OP_MULT));
      PtrCoefFct coefFreqFactor = NULL;
      coefFreqFactor = CoefFunction::Generate(mp_, Global::COMPLEX, "-2*pi*f");
      PtrCoefFct dispcurr = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, coefFreqFactor, tmp, CoefXpr::OP_MULT));

      dispCurrDensCoef->AddRegion(*regIt, dispcurr);
    }
  }

  std::map<SolutionType, shared_ptr<FeSpace>>
  FullWaveMaxwellAVPDE::CreateFeSpaces(const std::string &formulation,
                                     PtrParamNode infoNode)
  {
    std::map<SolutionType, shared_ptr<FeSpace>> crSpaces;
    PtrParamNode magVecPpotSpaceNode = infoNode->Get("magPotential");
    crSpaces[MAG_POTENTIAL] = FeSpace::CreateInstance(myParam_, magVecPpotSpaceNode, FeSpace::HCURL, ptGrid_);
    crSpaces[MAG_POTENTIAL]->Init(solStrat_);

    PtrParamNode elecScalPotSpaceNode = infoNode->Get("elecPotential");
    crSpaces[ELEC_POTENTIAL] = FeSpace::CreateInstance(myParam_, elecScalPotSpaceNode, FeSpace::H1, ptGrid_);
    crSpaces[ELEC_POTENTIAL]->Init(solStrat_);

    return crSpaces;
  }

} // end of namespace
