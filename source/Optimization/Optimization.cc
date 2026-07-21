#include <assert.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ResultHandler.hh"
#include "Domain/Domain.hh"
#include "Driver/Assemble.hh"
#include "Driver/SolveSteps/BaseSolveStep.hh"
#include "Driver/BaseDriver.hh"
#include "Driver/FormsContexts.hh"
#include "Driver/HarmonicDriver.hh"
#include "Driver/EigenFrequencyDriver.hh"
#include "Forms/LinForms/LinearForm.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignStructure.hh"
#include "Optimization/Design/FeatureMappingDesign.hh"
#include "Optimization/Design/SpaghettiDesign.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Function.hh"
#include "Optimization/LevelSet.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "Optimization/Optimizer/EvaluateOnly.hh"
#include "Optimization/Optimizer/GradientCheck.hh"
#include "Optimization/Optimizer/OptimalityCondition.hh"
#include "Optimization/Optimizer/ShapeOptimizer.hh"
#include "Optimization/Optimizer/MMA.hh"
#include "Optimization/Optimizer/DumasMMA.hh"
#include "Optimization/ParamMat.hh"
#include "Optimization/HeatSIMP.hh"
#include "Optimization/PiezoSIMP.hh"
#include "Optimization/MagSIMP.hh"
#include "Optimization/AcouSIMP.hh"
#include "Optimization/PiezoParamMat.hh"
#include "Optimization/FeatureMappingParamMat.hh"
#include "Optimization/SIMP.hh"
#include "Optimization/ShapeGrad.hh"
#include "Optimization/ShapeOpt.hh"
#include "Optimization/SplineBoxOpt.hh"
#include "Optimization/Transform.hh"
#include "Optimization/Tune.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/BasePDE.hh"
#include "Utils/ToolsFull.hh"
#include "Utils/PythonKernel.hh"
#include <filesystem>
#include "def_use_ipopt.hh"
#include "def_use_scpip.hh"
#include "def_use_snopt.hh"
#include "def_use_embedded_python.hh"
#include "def_use_sgp.hh"
#include "def_use_dumas.hh"

#ifdef USE_SGP
// check if Intel MKL is available
// need it for SGP
  #include <def_use_blas.hh>
  #include "Optimization/Optimizer/SGPHolder.hh"
#endif

// IPOPT, SCPIP and SnOpt are not necessarily linked
#ifdef USE_IPOPT
  #include "Optimization/Optimizer/IPOPTHolder.hh"
  #include "Optimization/Optimizer/FeasPP.hh"
#endif
#ifdef USE_SCPIP
  #include "Optimization/Optimizer/SCPIP.hh"
#endif
#ifdef USE_SNOPT
  #include "Optimization/Optimizer/SnOpt.hh"
#endif
#ifdef USE_EMBEDDED_PYTHON
  #include "Optimization/Optimizer/PythonOptimizer.hh"
#endif

using namespace CoupledField;
using namespace std;
namespace fs = std::filesystem;

DEFINE_LOG(opt, "opt")

// instantiation of the static elements
Enum<Optimization::Optimizer>        Optimization::optimizer;
Enum<App::Type>                      Optimization::application;
Enum<Optimization::CommitMode>       Optimization::commitMode;

Context*                             Optimization::context;
ContextManager                       Optimization::manager;

Optimization::Optimization()
{
  this->lastStoredResult_ = -1;
  this->design = NULL;
  this->baseOptimizer_ = NULL;
  this->currentIteration = 0; // a 1 or 0 can make a lot of difference! 0 is initial design!
  this->writeCounter_ = 0;
  this->problemSolvedCounter = 0;
  this->problemWithinIteration = 0;
  this->grid = domain->GetGrid();
  assert(domain->GetInfoRoot()->Get(ParamNode::SUMMARY)->Has("timer"));
  this->cfs_timer_ = domain->GetInfoRoot()->Get(ParamNode::SUMMARY)->Get("timer")->AsTimer();

  Optimization::manager.Init(); // there is also an init in DesignSpace

  optInfoNode = domain->GetInfoRoot()->Get("optimization");   // store our info results here
  PtrParamNode header = optInfoNode->Get(ParamNode::HEADER);
  optParamNode = domain->GetParamRoot()->Get("optimization"); // read our parameters from the xml file

  // in transient optimization one can specify the initial value as a solution to a static problem and a weight for it (just in tracking)
  firstStepStatic = optParamNode->Has("firstStepStatic");
  if(firstStepStatic)
    otherStepWeight = 1.0 - optParamNode->Get("firstStepStatic/weight")->As<Double>();
  else
    otherStepWeight = 1.0;

  // the tool to solve the optimization problem
  optimizer_ = optimizer.Parse(optParamNode->Get("optimizer/type")->As<string>());
  maxIterations = optParamNode->Get("optimizer/maxIterations")->As<Integer>();

  // might read a multiObjective problem
  objectives.Read(optParamNode->Get("costFunction"));
  PtrParamNode ino = optInfoNode->Get(ParamNode::HEADER)->Get("objective");
  // postoine objectives.ToInfo(ino) to PostInit()

  isMultiObjective_ = Objective::type.Parse(optParamNode->Get("costFunction/type")->As<std::string>()) == Objective::MULTI_OBJECTIVE;
  if(isMultiObjective_)
  {
    multiObjectiveType_ = Objective::multiObjType.Parse(optParamNode->Get("costFunction/multiObjective/type")->As<std::string>());
    multiObjectiveBeta_ = optParamNode->Get("costFunction/multiObjective/beta")->As<double>();
    ino->Get("multiObjective/type")->SetValue(Objective::multiObjType.ToString(multiObjectiveType_));
    if(multiObjectiveType_ != Objective::WEIGHTED_SUM)
      ino->Get("multiObjective/beta")->SetValue(multiObjectiveBeta_);
  }

  // multiple excitations are are toggled via attribute. Only if enabled we read the optional element
  // actually part of costFunction - but we store in Optimization itself!
  // theoretically we might have multiple multipleExcitation in the xml file for multi sequence cases.
  // however this is not implemented yet
  bool dme = optParamNode->Get("costFunction/multiple_excitation")->As<bool>();
  this->me = new MultipleExcitation(dme, dme ? optParamNode->Get("costFunction/multipleExcitation", ParamNode::PASS) : PtrParamNode());
  if(dme)
    me->ToInfo(header->Get("multipleExcitations"));

  if(manager.any().bloch && !dme)
    header->SetWarning("Bloch mode analysis but not multiple excitation activated");

  // slope constraints to be processed in SIMP -> Constraints::PostProc
  ParamNodeList list = optParamNode->GetList("constraint");
  constraints.Read(list);
  PtrParamNode in = header->Get("constraints");

  // the commit stuff
  string cm = optParamNode->Has("commit") ? optParamNode->Get("commit/mode")->As<string>() : "forward";
  this->commitMode_ = commitMode.Parse(cm);
  this->commitStride = optParamNode->Has("commit") ? optParamNode->Get("commit/stride")->As<Integer>() : 1;
  optInfoNode->Get("commit/mode")->SetValue(cm);
  optInfoNode->Get("commit/stride")->SetValue(commitStride);

  // write the HALTOPT filename, helps to memorize how to write it :)
  optInfoNode->Get("haltopt_file")->SetValue(fs::current_path().string() + "/HALTOPT");

  // remove a stop file, if found
  if(fs::exists("HALTOPT"))
  {
    bool good = fs::remove("HALTOPT");
    if(!good) throw new Exception("Could not remove file 'HALTOPT' after detection");
  }
}

Optimization::~Optimization()
{
  delete design; design = NULL;
  delete baseOptimizer_; baseOptimizer_ = NULL;
  delete me; me = NULL;
}

void Optimization::PostInit()
{
  // CalcGraynessScaling() cannot be called within constructor
  objectives.ToInfo(optInfoNode->Get(ParamNode::HEADER)->Get("objective"));


  // during Optimization construction there were no pdes (at least for multi sequence), now in PostInit() we need to read them
  for(unsigned int i = 1; i < manager.context.GetSize(); i++) // 0 is set below
    manager.SwitchContext(i); // driver and pdes are created once and then stored

  manager.SwitchContext(0); // go back to first which is what we expect.

  assert(context->pde != NULL);
}

void Optimization::PostInitSecond()
{
  log.AddToHeader("iter");

  if(manager.any().harmonic)
    log.AddToHeader("freq");

  // this is Daniel specific
  if(isMultiObjective_)
    log.AddToHeader("multiObjectiveValue");

  LOG_DBG(opt) << "PIS: me=" << me->IsEnabled() << " #me=" << me->excitations.GetSize() << " #obj=" << objectives.data.GetSize();

  for(unsigned int i = 0; i < objectives.data.GetSize(); i++)
  {
    const Objective* o = dynamic_cast<Objective*>(objectives.data[i]);
    log.AddToHeader(o->GetType() == Function::SLACK_FNCT ? Function::slackFnct.ToString(o->GetSlackFnct()) : o->ToString());
    if(o->GetType() == Function::BANDGAP)
    {
      log.AddToHeader("max_ef_" + to_string(o->bandgap.lower_ev) + "_wv");
      log.AddToHeader("min_ef_" + to_string(o->bandgap.upper_ev) + "_wv");
    }
    switch(o->GetTerm())
    {
    case Objective::PENALTY:
      log.AddToHeader(o->ToString() + "_penalty");
      break;
    case Objective::LN1P:
      log.AddToHeader(o->ToString() + "_ln1p");
      break;
    case Objective::POWER:
      log.AddToHeader(o->ToString() + "_power");
      break;
    case Objective::LINEAR:
      break; // no need to add
    }
  }

  if(me->IsEnabled() && me->excitations.GetSize() > 1)
    for(Excitation& ex : me->excitations)
      log.AddToHeader("objective_" + ex.GetFullLabel());

  log.AddToHeader("duration");

  if(design->HasAlphaVariable())
    log.AddToHeader("alpha");

  if(design->HasSlackVariable() && !objectives.Has(Function::SLACK))
    log.AddToHeader("slack");

  log.Init(this, optParamNode->Get("log")->As<string>(), optParamNode->Get("logging", ParamNode::PASS)); // is fail save

  // add the bandgap stuff in front of the constraints
  for(unsigned int i = 0; i < log.bloch_info.GetSize(); i++) // might be emtpy!
    log.AddToHeader(std::get<0>(log.bloch_info[i]));

  // constraints.ToInfo() is called in PostInitSecond()
  for(unsigned int i = 0; i < constraints.all.GetSize(); i++)
  {
    Condition* g = constraints.all[i];
    if(!g->IsLocalCondition())  {
      if(g->GetType() != Function::EIGENFREQUENCY || log.plot_ev)
        log.AddToHeader(g->ToString());
      if(g->GetType() == Function::EIGENFREQUENCY && g->GetExcitation()->DoBloch() && !g->DoFullBloch())
        log.AddToHeader("ef_" + std::to_string(g->GetEigenValueID()) + "_wv");
    }
    else {
      log.AddToHeader((g->GetBound() != Condition::LOWER_BOUND ? "max_abs_" : "min_abs_") + g->ToString());
      if(progOpts->DoDetailedInfo())
        log.AddToHeader("mean_abs_" + g->ToString());
      log.AddToHeader("infeas_count_" + g->ToString());
    }
    LOG_DBG2(opt) << "PIS: i=" << i << " g=" << g->ToString() << " gme=" << g->ToString() << " e=" << g->GetExcitation()->GetFullLabel() << " ei=" << g->GetExcitation()->index;
  }

  PtrParamNode opt = optParamNode->Get("optimizer");

  switch(optimizer_)
  {
    case IPOPT_SOLVER:
         #ifdef USE_IPOPT
           baseOptimizer_ = new IPOPTHolder(this, opt);
         #else
           throw Exception("openCFS was compiled w/o IPOPT!");
         #endif
         break;

    case SCPIP_SOLVER:
         #ifdef USE_SCPIP
           baseOptimizer_ = new SCPIP(this, opt);
         #else
           throw Exception("openCFS was compiled w/o SCPIP");
         #endif
         break;

    case FEAS_PP_SOLVER:
         #ifndef USE_IPOPT
           throw Exception("openCFS needs to be compiled with IPOPT to use feasPP");
         #else
           baseOptimizer_ = new FeasPP(this, opt);
         #endif
         break;

    case MMA_SOLVER:
         baseOptimizer_ = new MMA(this, opt); // our self written textbook MMA
         break;

    case DUMAS_MMA:
    case DUMAS_GCMMA:
         #ifdef USE_DUMAS
           baseOptimizer_ = new DumasMMA(this, opt, optimizer_); // C++ variant of Niels Aages' PETSc MMA or a GCMMA variant
         #else
           throw Exception("openCFS was compiled w/o Dumas (MMA/GCMMA)");
         #endif
         break;

    case SNOPT_SOLVER:
         #ifdef USE_SNOPT
           baseOptimizer_ = new SnOpt(this, opt);
         #else
           throw Exception("openCFS was compiled w/o SnOpt");
         #endif
         break;

    case PYTHON_SOLVER:
         #ifdef USE_EMBEDDED_PYTHON
           baseOptimizer_ = new PythonOptimizer(this, opt);
         #else
           throw Exception("openCFS was compiled w/o USE_EMBEDDED_PYTHON");
         #endif
         break;

    case SGP_SOLVER:
         #ifdef USE_SGP
           baseOptimizer_ = new SGPHolder(this, opt);
         #else
           throw Exception("openCFS was compiled w/o external SGP lib!");
         #endif

         break;

    case OCM_SOLVER:
         baseOptimizer_ = new OptimalityCondition(this, opt);
         break;

    case SHAPE_SOLVER:
         baseOptimizer_ = new ShapeOptimizer(this, opt);
         break;

    case EVALUATE_INITIAL_DESIGN:
         baseOptimizer_ = new EvaluateOnly(this, opt);
         break;

    case GRADIENT_CHECK:
         baseOptimizer_ = new GradientCheck(this, opt);
         break;
  }

  baseOptimizer_->PostInit();

  constraints.ToInfo(optInfoNode->Get(ParamNode::HEADER)->Get("constraints"));

  unsigned int n = design->GetNumberOfVariables();
  if(log.design)
  {
    for(unsigned int i = 0; i < n; ++i)
      log.AddToHeader("design");
  }
  if(log.designGradient)
  {
    for(unsigned int i = 0; i < objectives.data.GetSize(); ++i)
      for(unsigned int j = 0; j < n; ++j)
        log.AddToHeader("designGradient_" + objectives.data[i]->GetName());
  }
  if (this->log.designConstraintGradients)
  {
    for(unsigned int i = 0; i < constraints.all.GetSize(); i++)
    {
      Condition* g = constraints.all[i];
      if(!g->IsLocalCondition())
        for (unsigned int j = 0; j < n; ++j)
          log.AddToHeader("constraintGradient" + g->ToString());
    }
  }
  design->SetOptimizer(baseOptimizer_);
  // add plot logging of the optimizer
  baseOptimizer_->LogFileHeader(log);

  python->CallHook(PythonKernel::OPT_POST_INIT);
}


void Optimization::SetEnums()
{
  Function::type.SetName("Function::Type");
  Function::type.Add(Function::MULTI_OBJECTIVE, "multiObjective");
  Function::type.Add(Function::COMPLIANCE, "compliance");
  Function::type.Add(Function::OUTPUT, "output");
  Function::type.Add(Function::SQUARED_OUTPUT, "squaredOutput");
  Function::type.Add(Function::DYNAMIC_OUTPUT, "dynamicOutput");
  Function::type.Add(Function::DYNAMIC_OUTPUT_TRACKING, "dynamicOutputTracking");
  Function::type.Add(Function::REFLECTED_WAVE, "reflectedWave");
  Function::type.Add(Function::ABS_OUTPUT, "absOutput");
  Function::type.Add(Function::GLOBAL_DYNAMIC_COMPLIANCE, "globalDynamicCompliance");
  Function::type.Add(Function::CONJUGATE_COMPLIANCE, "conjugateCompliance");
  Function::type.Add(Function::VOLUME, "volume");
  Function::type.Add(Function::PENALIZED_VOLUME, "penalizedVolume");
  Function::type.Add(Function::GAP, "gap");
  Function::type.Add(Function::REALVOLUME, "realvolume");
  Function::type.Add(Function::TRACKING, "tracking");
  Function::type.Add(Function::ELEC_ENERGY, "elecEnergy");
  Function::type.Add(Function::ENERGY_FLUX, "energyFlux");
  Function::type.Add(Function::TEMP_TRACKING_AT_INTERFACE, "tempTrackingAtInterface");
  Function::type.Add(Function::HOM_TENSOR, "homTensor");
  Function::type.Add(Function::HOM_TRACKING, "homTracking");
  Function::type.Add(Function::HOM_FROBENIUS_PRODUCT, "homFrobeniusProduct");
  Function::type.Add(Function::POISSONS_RATIO, "poissonsRatio");
  Function::type.Add(Function::YOUNGS_MODULUS, "youngsModulus");
  Function::type.Add(Function::YOUNGS_MODULUS_E1, "youngsModulusE1");
  Function::type.Add(Function::YOUNGS_MODULUS_E2, "youngsModulusE2");
  Function::type.Add(Function::TYCHONOFF, "tychonoff");
  Function::type.Add(Function::TEMPERATURE, "temperature");
  Function::type.Add(Function::GREYNESS, "greyness");
  Function::type.Add(Function::FILTERING_GAP, "filteringGap");
  Function::type.Add(Function::GLOBAL_STRESS, "globalStress");
  Function::type.Add(Function::LOCAL_STRESS, "localStress");
  Function::type.Add(Function::ISOTROPY, "isotropy");
  Function::type.Add(Function::ISO_ORTHOTROPY, "iso-orthotropy");
  Function::type.Add(Function::ORTHOTROPY, "orthotropy");
  Function::type.Add(Function::SLOPE, "slope");
  Function::type.Add(Function::GLOBAL_SLOPE, "globalSlope");
  Function::type.Add(Function::PERIMETER, "perimeter");
  Function::type.Add(Function::MOLE, "mole");
  Function::type.Add(Function::GLOBAL_MOLE, "globalMole");
  Function::type.Add(Function::OSCILLATION, "oscillation");
  Function::type.Add(Function::GLOBAL_OSCILLATION, "globalOscillation");
  Function::type.Add(Function::JUMP, "jump");
  Function::type.Add(Function::GLOBAL_JUMP, "globalJump");
  Function::type.Add(Function::BUMP, "bump");
  Function::type.Add(Function::CURVATURE, "curvature");
  Function::type.Add(Function::GLOBAL_CURVATURE, "globalCurvature");
  Function::type.Add(Function::OVERHANG_VERT, "overhang_vert");
  Function::type.Add(Function::OVERHANG_HOR, "overhang_hor");
  Function::type.Add(Function::DISTANCE, "distance");
  Function::type.Add(Function::BENDING, "bending");
  Function::type.Add(Function::CONES, "cones");
  Function::type.Add(Function::DESIGN, "design");
  Function::type.Add(Function::GLOBAL_DESIGN, "globalDesign");
  Function::type.Add(Function::PERIODIC, "periodic");
  Function::type.Add(Function::DESIGN_TRACKING, "designTracking");
  Function::type.Add(Function::SUM_MODULI, "sumModuli");
  Function::type.Add(Function::GLOBAL_SUM_MODULI, "globalSumModuli");
  Function::type.Add(Function::TWO_SCALE_VOL, "twoScaleVolume");
  Function::type.Add(Function::GLOBAL_TWO_SCALE_VOL, "globalTwoScaleVolume");
  Function::type.Add(Function::ORTHOTROPIC_TENSOR_TRACE, "orthotropicTensorTrace");
  Function::type.Add(Function::GLOBAL_ORTHOTROPIC_TENSOR_TRACE, "globalOrthotropicTensorTrace");
  Function::type.Add(Function::TENSOR_TRACE, "tensorTrace");
  Function::type.Add(Function::GLOBAL_TENSOR_TRACE, "globalTensorTrace");
  Function::type.Add(Function::TENSOR_NORM, "tensorNorm");
  Function::type.Add(Function::PARAM_PS_POS_DEF, "parametrized-plane-stress-pos-def");
  Function::type.Add(Function::POS_DEF_DET_MINOR_1, "fmoPosDefMinor1");
  Function::type.Add(Function::POS_DEF_DET_MINOR_2, "fmoPosDefMinor2");
  Function::type.Add(Function::POS_DEF_DET_MINOR_3, "fmoPosDefMinor3");
  Function::type.Add(Function::BENSON_VANDERBEI_1, "bensonVanderbeiMinor1");
  Function::type.Add(Function::BENSON_VANDERBEI_2, "bensonVanderbeiMinor2");
  Function::type.Add(Function::BENSON_VANDERBEI_3, "bensonVanderbeiMinor3");
  Function::type.Add(Function::EIGENFREQUENCY, "eigenfrequency");
  Function::type.Add(Function::BUCKLING_LOAD_FACTOR, "bucklingLoadFactor");
  Function::type.Add(Function::LOCAL_BUCKLING_LOAD_FACTOR, "localBucklingLoadFactor");
  Function::type.Add(Function::GLOBAL_BUCKLING_LOAD_FACTOR, "globalBucklingLoadFactor");
  Function::type.Add(Function::MULTIMATERIAL_SUM, "multimaterial_sum");
  Function::type.Add(Function::SLACK, "slack");
  Function::type.Add(Function::SLACK_FNCT, "slackFunction");
  Function::type.Add(Function::BANDGAP, "bandgap");
  Function::type.Add(Function::SHAPE_INF, "shape_inf");
  Function::type.Add(Function::EXPRESSION, "expression");
  Function::type.Add(Function::PRESSURE_DROP, "pressureDrop");
  Function::type.Add(Function::HEAT_ENERGY, "heatEnergy");
  Function::type.Add(Function::SQR_MAG_FLUX_DENS_X, "sqrMagFluxDensX");
  Function::type.Add(Function::SQR_MAG_FLUX_DENS_Y,"sqrMagFluxDensY");
  Function::type.Add(Function::SQR_MAG_FLUX_DENS_RZ, "sqrMagFluxDensRZ");
  Function::type.Add(Function::LOSS_MAG_FLUX_RZ, "lossMagFluxRZ");
  Function::type.Add(Function::MAG_COUPLING,"magCoupling");
  Function::type.Add(Function::ARC_OVERLAP,"arc_overlap");
  Function::type.Add(Function::PYTHON_FUNCTION,"python");
  Function::type.Add(Function::LOCAL_PYTHON_FUNCTION,"localPython");

  Function::slackFnct.SetName("Function::SlackFnct");
  Function::slackFnct.Add(Function::NO_FUNCTION, "no_function");
  Function::slackFnct.Add(Function::ALPHA_SLACK_QUOTIENT, "a/s");
  Function::slackFnct.Add(Function::REL_BANDGAP, "(2*s)/(a-s)");
  Function::slackFnct.Add(Function::NORM_BANDGAP, "(2*s)/a");
  Function::slackFnct.Add(Function::ALPHA_MINUS_SLACK, "a-s");

  Function::multiObjType.SetName("Function::MultiObjType");
  Function::multiObjType.Add(Function::WEIGHTED_SUM, "weightedSum");
  Function::multiObjType.Add(Function::SMOOTH_MIN, "smoothMin");
  Function::multiObjType.Add(Function::SMOOTH_MAX, "smoothMax");

  Function::access.SetName("Function::Access");
  Function::access.Add(Function::PLAIN, "plain");
  Function::access.Add(Function::FILTERED, "filtered");
  Function::access.Add(Function::PHYSICAL, "physical");
  Function::access.Add(Function::DEFAULT, "default");

  Function::Local::locality.SetName("Function::Local::Locality");
  Function::Local::locality.Add(Function::Local::DEFAULT, "default");
  Function::Local::locality.Add(Function::Local::NEXT, "next");
  Function::Local::locality.Add(Function::Local::NEXT_DIAG, "next_diag");
  Function::Local::locality.Add(Function::Local::NEXT_AND_REVERSE, "next_and_reverse");
  Function::Local::locality.Add(Function::Local::PREV_NEXT, "prev_next");
  Function::Local::locality.Add(Function::Local::PREV_NEXT_AND_REVERSE, "prev_next_and_reverse");
  Function::Local::locality.Add(Function::Local::DEG_45_STAR, "45_deg_star");
  Function::Local::locality.Add(Function::Local::DEG_45_STAR_AND_REVERSE, "45_deg_star_and_reverse");
  Function::Local::locality.Add(Function::Local::BOUNDARY, "boundary");
  Function::Local::locality.Add(Function::Local::CYCLIC, "cyclic");
  Function::Local::locality.Add(Function::Local::ELEMENT, "element");
  Function::Local::locality.Add(Function::Local::MULT_DESIGNS_ELEMENT, "multiple_designs_element");
  Function::Local::locality.Add(Function::Local::MULT_DESIGNS_NEXT, "multiple_designs_next");
  Function::Local::locality.Add(Function::Local::MULT_DESIGNS_NEXT_AND_REVERSE, "multiple_designs_next_and_reverse");
  Function::Local::locality.Add(Function::Local::MULT_DESIGNS_PREV_NEXT, "multiple_designs_prev_next");
  Function::Local::locality.Add(Function::Local::MULT_DESIGNS_PREV_NEXT_AND_REVERSE, "multiple_designs_prev_next_and_reverse");
  Function::Local::locality.Add(Function::Local::SHAPE, "shape");
  Function::Local::locality.Add(Function::Local::FUNCTION_SPECIFIC, "function_specific");
  Function::Local::locality.Add(Function::Local::FUNCTION_SPECIFIC_TWO_SIGNS, "function_specific_two_signs");
  Function::Local::locality.Add(Function::Local::EXTERNALLY_DEFINED, "externally_defined");

  Function::Local::phase.SetName("Function::Local::Phase");
  Function::Local::phase.Add(Function::Local::BOTH, "both");
  Function::Local::phase.Add(Function::Local::VOID_MAT, "void");
  Function::Local::phase.Add(Function::Local::MATERIAL, "material");

  Function::stressType.SetName("Function::StressType");
  Function::stressType.Add(Function::MECH, "mech");
  Function::stressType.Add(Function::PIEZO, "piezo");
  Function::stressType.Add(Function::ONLY_COUPLING, "only_coupling");

  Objective::term.SetName("Objective::Term");
  Objective::term.Add(Objective::LINEAR, "linear");
  Objective::term.Add(Objective::PENALTY, "penalty");
  Objective::term.Add(Objective::LN1P, "ln1p");
  Objective::term.Add(Objective::POWER, "power");

  Condition::bound.SetName("Condition::Bound");
  Condition::bound.Add(Condition::EQUAL, "equal");
  Condition::bound.Add(Condition::LOWER_BOUND, "lowerBound");
  Condition::bound.Add(Condition::UPPER_BOUND, "upperBound");

  StoppingRule::type.SetName("StoppingRule::Type");
  StoppingRule::type.Add(StoppingRule::DESIGN_CHANGE, "designChange");
  StoppingRule::type.Add(StoppingRule::REL_COST_CHANGE, "relativeCostChange");
  StoppingRule::type.Add(StoppingRule::ABOVE_FUNCTION, "aboveFunction");
  StoppingRule::type.Add(StoppingRule::BELOW_FUNCTION, "belowFunction");
  StoppingRule::type.Add(StoppingRule::MAX_HOURS, "maxHours");
  StoppingRule::type.Add(StoppingRule::OSCILLATIONS, "oscillations");


  StoppingRule::condition.SetName("StoppingRule::Condition");
  StoppingRule::condition.Add(StoppingRule::SUFFICIENT, "sufficient");
  StoppingRule::condition.Add(StoppingRule::NECESSARY, "necessary");

  optimizer.SetName("Optimization::Optimizer");
  optimizer.Add(OCM_SOLVER, "ocm");
  optimizer.Add(IPOPT_SOLVER, "ipopt");
  optimizer.Add(SCPIP_SOLVER, "scpip");
  optimizer.Add(FEAS_PP_SOLVER, "feasPP");
  optimizer.Add(MMA_SOLVER, "mma");
  optimizer.Add(DUMAS_MMA, "dumas_mma");
  optimizer.Add(DUMAS_GCMMA, "dumas_gcmma");
  optimizer.Add(SGP_SOLVER, "sgp");
  optimizer.Add(SNOPT_SOLVER, "snopt");
  optimizer.Add(PYTHON_SOLVER, "python");
  optimizer.Add(SHAPE_SOLVER, "shapeOpt");
  optimizer.Add(EVALUATE_INITIAL_DESIGN, "evaluate");
  optimizer.Add(GRADIENT_CHECK, "gradientCheck");

  ErsatzMaterial::method.SetName("ErsatzMaterial::Method");
  ErsatzMaterial::method.Add(ErsatzMaterial::SIMP_METHOD, "simp");
  ErsatzMaterial::method.Add(ErsatzMaterial::PARAM_MAT, "paramMat");
  ErsatzMaterial::method.Add(ErsatzMaterial::SHAPE_GRAD, "shapeGrad");
  ErsatzMaterial::method.Add(ErsatzMaterial::SHAPE_OPT, "shapeOpt");
  ErsatzMaterial::method.Add(ErsatzMaterial::SHAPE_PARAM_MAT, "shapeParamMat");
  ErsatzMaterial::method.Add(ErsatzMaterial::SHAPE_MAP, "shapeMap");
  ErsatzMaterial::method.Add(ErsatzMaterial::SPAGHETTI, "spaghetti");
  ErsatzMaterial::method.Add(ErsatzMaterial::SPAGHETTI_PARAM_MAT, "spaghettiParamMat");
  ErsatzMaterial::method.Add(ErsatzMaterial::SPLINE_BOX, "splineBox");
  ErsatzMaterial::method.Add(ErsatzMaterial::FEATURE_MAPPING, "featureMapping");
  ErsatzMaterial::method.Add(ErsatzMaterial::FEATURE_MAPPING_PARAM_MAT, "featureMappingAniso");

  ErsatzMaterial::commitMode.SetName("ErsatzMaterial::CommitMode");
  ErsatzMaterial::commitMode.Add(ErsatzMaterial::FORWARD, "forward");
  ErsatzMaterial::commitMode.Add(ErsatzMaterial::EACH_FORWARD, "each_forward");
  ErsatzMaterial::commitMode.Add(ErsatzMaterial::ADJOINT, "adjoint");
  ErsatzMaterial::commitMode.Add(ErsatzMaterial::EACH_ADJOINT, "each_adjoint");
  ErsatzMaterial::commitMode.Add(ErsatzMaterial::BOTH, "both_cases");

  OptimizationMaterial::system.SetName("OptimizationMaterial::System");
  OptimizationMaterial::system.Add(OptimizationMaterial::PIEZOCOUPLING, "piezo");
  OptimizationMaterial::system.Add(OptimizationMaterial::MECH, "mechanic");
  OptimizationMaterial::system.Add(OptimizationMaterial::HEAT, "heat");
  OptimizationMaterial::system.Add(OptimizationMaterial::MAG, "magnetic");
  OptimizationMaterial::system.Add(OptimizationMaterial::ACOUSTIC, "acoustic");
  OptimizationMaterial::system.Add(OptimizationMaterial::LBM, "lbm");

  application.SetName("App::Type");
  application.Add(App::NO_APP, "no_app");
  application.Add(App::ACOUSTIC, "acoustic");
  application.Add(App::HEAT, "heat");
  application.Add(App::MAG, "magnetic");
  application.Add(App::LAPLACE, "laplace");
  application.Add(App::MECH, "mech");
  application.Add(App::BUCKLING, "buckling");
  application.Add(App::MASS, "mass");
  application.Add(App::ELEC, "elec");
  application.Add(App::PIEZO_COUPLING, "piezoCoupling");
  application.Add(App::PRESSURE, "pressure");
  application.Add(App::CHARGE_DENSITY, "chargeDensity");
  application.Add(App::STRESS, "stress");
  application.Add(App::LBM, "lbm");

  LevelSet::Action::type.SetName("LevelSet::Action::Type");
  LevelSet::Action::type.Add(LevelSet::Action::SIGNED_DISTANCE_FIELD, "signedDistanceField");
  LevelSet::Action::type.Add(LevelSet::Action::TRIVIAL_HOLE, "trivialHole");
  LevelSet::Action::type.Add(LevelSet::Action::DO_SHAPE_STEP, "shapeStep");

  MultipleExcitation::type.SetName("MultipleExcitation::Type");
  MultipleExcitation::type.Add(MultipleExcitation::NO_TYPE, "no_type");
  MultipleExcitation::type.Add(MultipleExcitation::FIXED_WEIGHT, "fixed_weights");
  MultipleExcitation::type.Add(MultipleExcitation::META_OBJECTIVE, "meta_objective");
  MultipleExcitation::type.Add(MultipleExcitation::HOMOGENIZATION_TEST_STRAINS, "homogenizationTestStrains");

  Transform::type.SetName("Transform::Type");
  Transform::type.Add(Transform::ROTATION, "rotate");
}

bool Optimization::IsTransient() {
  return(domain->GetDriver()->GetAnalysisType() == BasePDE::TRANSIENT);
}

double Optimization::GetStepWeight(unsigned int ts) const{
  unsigned int nts = context->GetDriver()->GetNumSteps();
  if(IsFirstTransientStepStatic()){
    if(ts == 0){
      return((1.0 - otherStepWeight));
    }else{
      return(otherStepWeight / (nts-1));
    }
  }else{
    return(1.0 / nts);
  }
}


PtrParamNode Optimization::DoStopOptimizationHelper(bool converged, const string& reason)
{
  PtrParamNode in = optInfoNode->Get(ParamNode::SUMMARY)->Get("break");

  if(reason != "")
  {
    in->Get("converged")->SetValue(converged ? "yes" : "no");
    user_break_reason = reason;
    in->Get("reason/msg")->SetValue(reason);
  }
  return in;
}

bool Optimization::DoStopOptimization()
{
  // only PythonStopOptimization() will set, we are done
  if(user_break_reason != "")
  {
    DoStopOptimizationHelper();
    return true;
  }

  // check if the HALTOPT file exists
  if(fs::exists("HALTOPT"))
  {
    bool good = fs::remove("HALTOPT");
    if(!good)
      throw new Exception("Could not remove file 'HALTOPT' after detection");
    DoStopOptimizationHelper(false, "Detected file 'HALTOPT'");
    return true;
  }

  StdVector<string> reasons;
  reasons.Reserve(objectives.stop.GetSize());

  bool not_converged = false;
  // convergence with sufficient and necessary is tricky
  // any sufficient breaks
  bool must_break = false; // a sufficient rule is triggered
  // when we have necessary, all necessary need to be true
  int n_nec = -1; // allows comparison with cnt_nec not be true for 0 necessary conditions
  for(auto& rule : objectives.stop)
    if(rule.GetCondition() == rule.NECESSARY)
      n_nec = std::max(n_nec,0) + 1;

  // this counts the active sufficient conditions
  int cnt_nec = 0;

  for(auto& rule : objectives.stop)
  {
    string reason = rule.DoStop(&objectives, &constraints, time_);

    if(rule.GetCondition() == rule.SUFFICIENT && reason != "")
      must_break = true;
    if(rule.GetCondition() == rule.NECESSARY && reason != "")
      cnt_nec++;
    if(reason != "")
    {
      reasons.Push_back(reason);
      if(rule.GetType() == rule.MAX_HOURS || rule.GetType() == rule.OSCILLATIONS)
        not_converged = true;
    }

    LOG_DBG(opt) << "DSO rule=" << rule.type.ToString(rule.GetType()) << " r='" << reason << "' cond=" << rule.condition.ToString(rule.GetCondition())
                 << " must_break=" << must_break << " n_nec=" << n_nec << " cnt_nec=" << cnt_nec;
  }

  // n_nec is -1 but not 0, therefore no necessary condition does not break with cnt_nec = 0
  if(objectives.stop.GetSize() > 0 && (must_break || (n_nec == cnt_nec)))
  {
    PtrParamNode in = optInfoNode->Get(ParamNode::SUMMARY)->Get("break");

    in->Get("converged")->SetValue(not_converged ? "no" : "yes");
    for(string rsn : reasons)
      in->Get("reason", ParamNode::APPEND)->Get("msg")->SetValue(rsn);
    return true;
  }
  return false;
}

StdVector<std::pair<string,string> > Optimization::GetStoppingRules() const
{
  StdVector<std::pair<string,string> > map;
  for(const auto& rule : objectives.stop)
  {
    string type = rule.type.ToString(rule.GetType());
    if(rule.GetType() == rule.BELOW_FUNCTION || rule.GetType() == rule.ABOVE_FUNCTION)
      type += "_" + rule.function;
    map.Push_back(std::make_pair(type, std::to_string(rule.value)));
  }
  return map;
}


OptimizationMaterial::System Optimization::ParseSystem()
{
  return OptimizationMaterial::system.Parse(domain->GetParamRoot()->Get("optimization/ersatzMaterial/material")->As<string>());
}

/** read only the very basic stuff */
Optimization* Optimization::CreateInstance()
{
  // set the enums we need
  Optimization::SetEnums(); // sets also ErsatzMaterial::Method
  DesignElement::SetEnums();
  DesignMaterial::SetEnums();

  PtrParamNode param = domain->GetParamRoot();

  if(!param->Has("optimization")) return NULL;

  // we assume ersatz material, currently there is nothing else.
  // note, we read method again in the ersatz material constructor.
  PtrParamNode em = param->Get("optimization/ersatzMaterial");

  ErsatzMaterial::Method method = ErsatzMaterial::method.Parse(em->Get("method")->As<string>());
  OptimizationMaterial::System material = ParseSystem();

  Optimization* opt = NULL;

  switch(method)
  {
  case ErsatzMaterial::SIMP_METHOD:
  case ErsatzMaterial::SHAPE_MAP: 
  case ErsatzMaterial::FEATURE_MAPPING:
  case ErsatzMaterial::SPAGHETTI: // we have also SPAGHETTI_PARAM_MAT
  {
    switch(material)
    {
    case OptimizationMaterial::MECH:
    case OptimizationMaterial::ELEC:
    case OptimizationMaterial::LBM:
      opt = new SIMP(); // generally single PDE!
      break;

    case OptimizationMaterial::HEAT:
      opt = new HeatSIMP(); 
      break;

    case OptimizationMaterial::ACOUSTIC:
      opt = new AcouSIMP();
      break;

    case OptimizationMaterial::PIEZOCOUPLING:
      opt = new PiezoSIMP();
      break;

    case OptimizationMaterial::MAG:
      opt = new MagSIMP();
      break;

    default:
      assert(false);
      break;
    }
    break;
  } // simp, shape map, feature mapping

  // FMO, ShapeGrad, ...
  case ErsatzMaterial::SPAGHETTI_PARAM_MAT: // for anisotroy we use ParamMat stuff
  case ErsatzMaterial::PARAM_MAT:
    if(material != OptimizationMaterial::PIEZOCOUPLING)
      opt = new ParamMat();
    else
      opt = new PiezoParamMat();
    break;
  case ErsatzMaterial::FEATURE_MAPPING_PARAM_MAT:
      opt = new FeatureMappingParamMat();
    break;
  case ErsatzMaterial::SHAPE_OPT:
  case ErsatzMaterial::SHAPE_PARAM_MAT:
    opt = new ShapeOpt();
    break;
  case ErsatzMaterial::SHAPE_GRAD:
    opt = new ShapeGrad();
    break;
  case ErsatzMaterial::SPLINE_BOX:
    opt = new SplineBoxOpt();
    break;
  default: throw Exception("Optimization not implemented");
  }

  // we have to do this, as PostInitSecond does already run CalcObjective/Gradient
  domain->SetOptimization(opt);

  return opt;
}

void Optimization::SolveProblem()
{
  // one driver is one multisequence step. We do this stuff here
  // and call the driver->StoreResults() multiple times

  ResultHandler* rh = NULL;

  if(!IsTransient()){ // transient optimization saves results in a different way
    rh = domain->GetResultHandler();
    unsigned int mss = context->GetDriver()->GetActSequenceStep();
    // max steps is high. The number is only relevant for hdf5, but there a hard limit
    rh->BeginMultiSequenceStep(mss, BasePDE::TRANSIENT, 9999);
  }

  Exception* e = NULL;
  try
  {
    baseOptimizer_->ToInfo(baseOptimizer_->GetInfoNode());
    baseOptimizer_->SolveOptimizationProblem();
    assert(baseOptimizer_->ValidateTimers());
    baseOptimizer_->ToInfo(baseOptimizer_->GetInfoNode());

    PtrParamNode fcts = optInfoNode->Get(ParamNode::SUMMARY)->Get("functions");
    for(Function* f : GetFunctions(false)) // objectives, constraints, observes
    {
      if(f->IsLocal() || f->history.GetSize() == 0)
        continue;
      PtrParamNode fin = fcts->Get(f->ToString());
      fin->Get("min")->SetValue(f->history.Min());
      fin->Get("max")->SetValue(f->history.Max());
      fin->Get("oscillations")->SetValue(f->CountOscillations());
    }

    if(optInfoNode->Get(ParamNode::SUMMARY)->Has("break"))
    {
      PtrParamNode in = optInfoNode->Get(ParamNode::SUMMARY)->Get("break");
      std::cout << "converged: " << in->Get("converged")->As<string>() << std::endl;
      for(auto rsn : in->GetList("reason"))
        std::cout << "reason: " << rsn->Get("msg")->As<string>() << std::endl;
    }
  }
  catch(Exception& ex)
  {
    e = new Exception(ex); // create new exception, don't throw now but let Bastian do his transient stuff
  }

  if(!IsTransient()){ // transient optimization saves results in a different way
    // do the finally - try to write results even if the optimizer broke down
    FinalizeStoreResults(); // when we have strides the results are written
    rh->FinishMultiSequenceStep();
    rh->Finalize();
  }
  if(e != NULL)
    throw *e;
  delete e;
}

bool Optimization::DoSolveAdjointWithState() const
{
  LOG_DBG(opt) << "DSAWS: buckling=" << context->DoBuckling() << " ms=" << context->DoMultiSequence() << " lbm=" << context->DoLBM()
    << " cplx=" << context->IsComplex() << " #ex=" << me->excitations.GetSize() << " #freq=" << me->GetUniqueFrequencies() << " rob=" << me->DoRobust(context);
  if(context->DoBuckling())
    return false;

  // easy case
  if(context->DoMultiSequence() || context->DoLBM())
    return true;

  // don't do it within forward when we can do it later
  if(context->IsComplex() && me->excitations.GetSize() > 1)
    return true;

  // we want to reuse the assembled system
  if(me->DoRobust(context))
    return true;  // the excitation needs to be checked externally

  else
    return false;
}


void Optimization::SolveStateProblem(Excitation* excite)
{
  assert(baseOptimizer_);
  assert(baseOptimizer_->ValidateTimers());

  // do not add the time solving the system to eval_[grad]_obj/constr_timer -> performance.py
  shared_ptr<Timer> eval_timer = baseOptimizer_ != NULL ? baseOptimizer_->GetRunningEvalTimer() : shared_ptr<Timer>();
  if(eval_timer)
    eval_timer->Stop();

  AnalysisID& id = context->driver->GetAnalysisId();
  id.iteration = currentIteration;

  assert(excite != NULL);
  assert(!(!me->IsEnabled() && excite->label == ""));

  if(excite->reassemble)
    context->pde->GetAssemble()->ResetMatrixReassembly();

  id.excite = me->IsEnabled() ? excite->GetFullLabel() : "";
  id.adjoint = false;

  if(IsTransient() && problemSolvedCounter > 0){ // transient optimization always has a mech pde
    SinglePDE* mech = context->ToPDE(App::MECH);
    assert(false);
    // FIXME mech->ReReadResults();
    design->AppendOptimizationResults(mech, true);
    assert(false);
    // FIXME mech->GetSolveStep()->ReInit();
  }

  // Do not store the results. This is to be done in CommitIteration
  if(context->IsHarmonic() && excite != NULL)
  {
    LOG_DBG(opt) << "SSP: harmonic step=" << excite->f_link->step << " f=" << excite->f_link->freq;
    context->GetHarmonicDriver()->ComputeFrequencyStep(*(excite->f_link));
  }
  else if(context->DoBloch())
  {
    LOG_DBG(opt) << "SSP: bloch step=" << excite->wave_vector.ToString() << " ex=" << excite->index << " wn=" << excite->GetWaveNumber() << " seq=" << context->sequence;
    context->GetEigenFrequencyDriver()->ComputeBlochWaveVector(excite->GetWaveNumber());
  }
  else
  {
    assert(!context->DoBloch() || !context->IsHarmonic());
    context->driver->SolveProblem();
      // FIXME driver->SolveProblem(IsTransient(), analysis_id, NULL); // static and transient optimization
  }

  if(eval_timer)
    eval_timer->Start();

  problemSolvedCounter++;
  problemWithinIteration++;
}

void Optimization::SolveAdjointProblems(Excitation* excite)
{
  assert(baseOptimizer_->ValidateTimers());

  // solve for objectives and constraints
  // do not solve for observe
  StdVector<Function*> ff = GetFunctions(false);

  for(unsigned int i = 0; i < ff.GetSize(); ++i)
  {
    Function* f = ff[i];
    assert(f != NULL);
    if(f->ctxt == context && f->IsAdjointBased() && f->DoEvaluate(excite) && !f->IsLocal() && (f->IsObjective() || !dynamic_cast<Condition*>(f)->IsObservation()))
      SolveAdjointProblem(excite, f); // virtual! calls ErsatzMaterial implementation
  }
}


StdVector<Function*> Optimization::GetFunctions(bool only_active) const
{
  StdVector<Function*> result;

  const unsigned int cn = objectives.data.GetSize();
  const unsigned int gn = constraints.active.GetSize();
  const unsigned int on = only_active ? 0 : constraints.observe.GetSize();

  result.Reserve(cn + gn + on);
  result.Resize(0); // To allow push back

  for(unsigned int i = 0; i < cn; i++)
  {
    result.Push_back(objectives.data[i]);
    LOG_DBG2(opt) << "GAF: o=" << result.Last()->ToString();
  }
  for(unsigned int i = 0; i < gn; i++)
  {
    result.Push_back(constraints.active[i]);
    LOG_DBG2(opt) << "GAF: g=" << result.Last()->ToString();
  }

  for(unsigned int i = 0; i < on; i++)
  {
    result.Push_back(constraints.observe[i]);
    LOG_DBG2(opt) << "GAF: g_observe=" << result.Last()->ToString();
  }

  return result;
}

Function* Optimization::GetFunction(const std::string& name, bool throw_exception)
{
  Function* f = objectives.Get(name, false); // no exception
  if(f == nullptr)
    f = (Function*) constraints.Get(name, throw_exception); // now exception if we don't have it

  if(f == nullptr && throw_exception)
    throw Exception("unknown function name '" + name + "'");
  return f;
}

ParamNodeList Optimization::GetMultipleExcitionsNodes(){
  ParamNodeList nodes;
  PtrParamNode pn = optParamNode->Get("costFunction");
  if(pn->Has("multipleExcitation/excitations"))
    nodes = pn->Get("multipleExcitation/excitations")->GetChildren();
  return nodes;
}

Tune* Optimization::SearchTune(Tune::Usage usage, bool silent)
{
  for(Tune* t : tunes)
    if(t->GetUsage() == usage)
      return t;

  if(!silent)
    throw Exception("none of the " + std::to_string(tunes.GetSize()) + " registered tunes is of usage " + Tune::usage.ToString(usage));
  return nullptr;
}

double Optimization::CalcSymmetry(DesignElement::Type de, DesignElement::ValueSpecifier vs, DesignElement::Access access)
{
  // the symmetry works only for squared models with a horizontal symmetry axis

  // our special result index
  int res_idx = design->GetSpecialResultIndex(de, vs, DesignElement::SYMMETRY, access);

  // plausibility check for squared
  int edge = (int) sqrt(design->data.GetSize());
  if(edge * edge != (int) design->data.GetSize())
    throw Exception("Symmetry not possible as models seems to be not squared");

  int max = (int) design->data.GetSize();
  double sum = 0;
  // we assume the first element (1) be in the lower left corner and then a lexical
  // ordering to right and then the rows up.
  for(int i = 0; i < edge/2; i++)
  {
    for(int j = 0; j < edge; j++)
    {
      // our data index
      int idx = i*edge + j;
      double idx_val = design->data[idx].GetValue(vs, access);

      // our counterpart
      int cntr = max - (i*edge)-(edge-j);
      double cntr_val = design->data[cntr].GetValue(vs, access);

      double err = (idx_val-cntr_val) / idx_val;
      LOG_DBG3(opt) << "Symmetry " << DesignElement::valueSpecifier.ToString(vs)
                    << "(" << DesignElement::access.ToString(access) << "): "
                    << idx << " (" << idx_val << ") : " << cntr << " (" << cntr_val << ") -> " << err;

      if(res_idx >= 0) design->data[idx].specialResult[res_idx] = err;
      sum += err;
    }
  }

  return sum / (double) max;
}

double Optimization::CalcObjective(Excitation* ev_only_excite)
{
  // in objective.value_ we store the sum over all excitations w/o penalty but with normalization
  // in excitation.cost we store the sum over all objectives with penalty but w/o normalization

  // reset the objective values such that we can sum up normalized but unpenalized values
  for(unsigned int o = 0; o < objectives.data.GetSize(); o++)
    objectives.data[o]->ResetValue();

  double result = 0.0;

  // term = linear -> function value, penalty max(value - parameter)^2, ...
  // the smooth max (from Daniel) want both approaches
  Vector<double> s_t_ov;   //  scale * term(objective value)
  Vector<double> w_s_t_ov; // weight * scale * term(objective value)

  // the multiple excitation case is a special case - for all other cases this is executed once
  for(unsigned int e = 0; e < (ev_only_excite != NULL ? 1 : me->excitations.GetSize()); e++)
  {
    Excitation& excite = ev_only_excite != NULL ? *ev_only_excite : me->excitations[e];
    excite.Apply(true); // sets the corresponding context
    excite.cost = 0.0;

    if(isMultiObjective_ && multiObjectiveType_ != Objective::WEIGHTED_SUM)
    {
      s_t_ov.Resize(objectives.data.GetSize());
      w_s_t_ov.Resize(objectives.data.GetSize());
      s_t_ov.Init();
      w_s_t_ov.Init();
    }

    for(unsigned int oi = 0; oi < objectives.data.GetSize(); oi++)
    {
      Objective* o = objectives.data[oi];

      // some objectives are only to be evaluated for the last excitation
      if(!o->DoEvaluate(&excite))
        continue;

      double ov = CalcFunction(excite, o, false); // this is virtual!
      o->SetOrgValue(ov);
      double tov = 0.0;
      switch(o->GetTerm())
      {
        case Objective::LINEAR:
          tov = ov;
          break;
        case Objective::PENALTY:
          tov = o->CalcPenalty(ov);
          break;
        case Objective::LN1P:
          tov = o->CalcLn1p(ov);
          break;
        case Objective::POWER:
          tov = std::pow(ov, o->GetParameter());
          break;
      }
      excite.cost += o->scale * tov;

      // we ignore the weight if the evaluation happens only once! TODO why not omega*omega? - Fabian
      double weight = !o->DoEvaluateAlways(excite.sequence) ? 1.0 : excite.normalized_weight;

      o->AddValue(tov * weight);

      result += weight * o->scale * tov;

      // if multiobjective, store function values
      if(isMultiObjective_ && (multiObjectiveType_ == Objective::SMOOTH_MAX || multiObjectiveType_ == Objective::SMOOTH_MIN))
      {
        s_t_ov[oi] = o->scale * tov;
        w_s_t_ov[oi] = weight * s_t_ov[oi];
      }

      LOG_DBG(opt) << "CO: ex=" << e << " obj=" << o->type.ToString(o->GetType()) << " ov=" << ov << " tov=" << tov << " param=" << o->GetParameter()
          << " term=" << o->term.ToString(o->GetTerm()) <<" scale=" << o->scale << " ex.cost=" << excite.cost << " nw=" << excite.normalized_weight
          << " weight=" << weight << " f->val=" << o->GetValue() << " result=" << result;
    }

    // if multiobjective, combine stored function values and overwrite results
    // case WEIGHTED_SUM has already been handled -> no overwrite necessary
    if(isMultiObjective_ && multiObjectiveType_ == Objective::SMOOTH_MIN)
    {
      excite.cost = SmoothMin(s_t_ov, multiObjectiveBeta_);
      result = SmoothMin(w_s_t_ov, multiObjectiveBeta_); // whyever we overwrite and do not sum up?! - Fabian
    }
    if(isMultiObjective_ && multiObjectiveType_ == Objective::SMOOTH_MAX)
    {
      excite.cost = SmoothMax(s_t_ov, multiObjectiveBeta_);
      result = SmoothMax(w_s_t_ov, multiObjectiveBeta_);
    }
  }
  calcObjIteration_ = this->GetCurrentIteration();

  return result;
}

void Optimization::CalcObjectiveGradient(StdVector<double>* grad_out, Excitation* ev_only_excite)
{
  assert(!baseOptimizer_ || baseOptimizer_->GetRunningEvalTimer()->IsRunning());
  // reset the cost gradients in the design elements and sum them up in a weighted way
  // to perform multiple loads
  design->Reset(DesignElement::COST_GRADIENT);

  for(unsigned int obj = 0; obj < objectives.data.GetSize(); obj++)
  {
    Objective* cost = objectives.data[obj];
    // the multiple excitation case is a special case - for all other cases this is executed once
    for(unsigned int idx = 0; idx < (ev_only_excite != NULL ? 1 : cost->ctxt->excitations.GetSize()); idx++)
    {
      Excitation* excite = ev_only_excite != NULL ? ev_only_excite : cost->ctxt->excitations[idx];

      // some objectives are only to be evaluated for the last excitation
      if(!cost->DoEvaluate(excite))
        continue;
      excite->Apply(true); // set the correct context

      CalcFunction(*excite, cost, true); // calls BaseDesignElement::AddGradient() which handles scale and term
    }
  }

  if(grad_out != NULL)
  {
    design->WriteGradientToExtern(*grad_out, DesignElement::COST_GRADIENT, DesignElement::SMART,  objectives.data[0]); // use the first such that we know about the robust index
    if(progOpts->DoDetailedInfo())
      design->WriteGradientFile(); // if constraints are not calculated yet will be overwritten later with the good data for this iterations
  }
}

double Optimization::CalcConstraint(Condition* g, Excitation* ev_only_excite)
{
  LOG_DBG2(opt) << " CC g=" << ( g != NULL ? g->ToString() : "null") <<"  eoe=" << ( ev_only_excite != NULL ? ev_only_excite->label : "null");

  // assume when we have only one constraint which is not explicitly given, this is not the stress constraint!
  assert((g == NULL && constraints.active.GetSize() == 1 && constraints.active[0]->DoEvaluateAlways(1) && !context->DoMultiSequence()) || g != NULL); // DoEvaluateAlways(): there is only one sequence

  if(g == NULL)
    g = constraints.active[0];

  double result = 0.0;

  for(unsigned int e = 0; e < (ev_only_excite != NULL ? 1 : me->excitations.GetSize()); e++)
  {
    Excitation& excite = ev_only_excite != NULL ? *ev_only_excite : me->excitations[e];
    excite.Apply(true); // switch context too for stuff like robust
    // in the evaluate once case only the last excitation
    double v = g->DoEvaluate(&excite) ? CalcFunction(excite, g, false) : 0.0;
    double w = g->DoEvaluateAlways(excite.sequence) ? excite.GetWeightedFactor(g) : 1.0;
    result += v * w;
    LOG_DBG2(opt) << "CC ex=" << e << " eval=" << g->DoEvaluate(&excite) << " v=" << v << " alw=" << g->DoEvaluateAlways(excite.sequence) << " w=" << w << " -> " << result;
  }

  g->SetValue(result);
  return result;
}

void Optimization::CalcConstraintGradient(Condition* g, StdVector<double>* grad_out, Excitation* ev_only_excite)
{
  // assume when we have only one constraint which is not explicitly given, this is not the stress constraint!
  // TODO: disable this assert as multi sequence cannot ruled out for every case
//  assert((g == NULL && constraints.active.GetSize() == 1 && !constraints.active[0]->DoEvaluateAlways(1) && !context->DoMultiSequence()) || g != NULL);

  if(g == NULL)
    g = constraints.active[0];

  for(unsigned int i = 0; i < (ev_only_excite != NULL ? 1 : g->ctxt->excitations.GetSize()); i++)
  {
    Excitation* ex = ev_only_excite != NULL ? ev_only_excite : g->ctxt->excitations[i];
    if(g->DoEvaluate(ex))
    {
      ex->Apply(true); // switch context if necessary

      if(context->DoBuckling())
      {
        unsigned int linElaExIndex = ex->index - (me->DoHomogenization() ? me->GetNumberHomogenization(context->ToApp()) : 1);
        assert(manager.GetContext(&(me->excitations[linElaExIndex])).driver->GetAnalysisType() == BasePDE::STATIC);
        ex->SetStressCoefFct( ex->GetStressCoefFctFromExcitation(linElaExIndex) );
      }

      CalcFunction(*ex, g, true);
    }
  }

  // copies from the design element gradient data to a memory array for external optimizers
  if(grad_out != NULL)
    design->WriteGradientToExtern(*grad_out, DesignElement::CONSTRAINT_GRADIENT, DesignElement::SMART, g);

  // check if we have constraint gradient as output
  // <result value="constraintGradient" detail="volume" access="plain" id="optResult_1"/>
  int res_idx = -1;
  // ALL_DESIGN becomes all
  int n    = design->FindDesign(g->GetDesignType(), false) == -1 ? design->data.GetSize() : design->GetNumberOfElements();
  int base = std::max(design->FindDesign(g->GetDesignType(), false), 0); // shift -1 to 0 for ALL_DESIGN
  if(DesignElement::detail.IsValid(g->type.ToString(g->GetType()))) // is current condition defined in the schema as detail
  {
    DesignElement::Detail detail = DesignElement::detail.Parse(g->type.ToString(g->GetType()));
    res_idx = design->GetSpecialResultIndex(g->GetDesignType(), DesignElement::CONSTRAINT_GRADIENT, detail, DesignElement::PLAIN); // TODO: excitation?
    if(res_idx > -1)
      for(int i = n * base; i < n * (base + 1); i++)
        design->data[i].specialResult[res_idx] = design->data[i].GetPlainGradient(g);
    res_idx = design->GetSpecialResultIndex(g->GetDesignType(), DesignElement::CONSTRAINT_GRADIENT, detail, DesignElement::SMART);
    if(res_idx > -1)
      for(int i = n * base; i < n * (base + 1); i++)
        design->data[i].specialResult[res_idx] = design->data[i].GetValue(DesignElement::CONSTRAINT_GRADIENT, DesignElement::SMART, g);
  }
}

void Optimization::EvaluateSpecialResults()
{
  for(unsigned int i = 0; i < design->resultDescriptions.GetSize(); i++)
  {
    const ResultDescription& rd = design->resultDescriptions[i];

    if(rd.detail == DesignElement::SYMMETRY)
      CalcSymmetry(rd.design, rd.value, rd.access);
  }
}


void Optimization::StoreResults(double step_val)
{
  // For PiezoSIMP we can do storing there and this method is overwritten
  // and might do nothing

  unsigned int internalWriteCounter;

  // this will write the CFS result and history file
  if(!IsTransient())
  { // transient optimization saves results in a different way
    if(step_val == -1)
      internalWriteCounter = context->GetDriver()->StoreResults(writeCounter_, currentIteration);
    else
      internalWriteCounter = context->GetDriver()->StoreResults(writeCounter_, step_val);

    if (!context->GetDriver()->GetResultHandler()->streamOnly)
      writeCounter_ = internalWriteCounter + 1;
  }
}

void Optimization::FinalizeStoreResults()
{
  // after the last CommitIteration the iteration counter was incremented
  bool store = (int) currentIteration-1 != lastStoredResult_ && currentIteration > 1;
  LOG_DBG(opt) << "CheckFinalStoreResults: currentIteration=" << currentIteration << " lastStoredResult="
               << lastStoredResult_ << " store=" << store << " writeCounter:" << writeCounter_;
  if(store)
    StoreResults(currentIteration-1);
}


PtrParamNode Optimization::CommitIteration()
{
  assert(!baseOptimizer_ || baseOptimizer_->ValidateTimers());
  // store the real cost -> not a scaled one
  objectives.PushBackHistory();

  // store the current design and calculate the design change!
  objectives.PushBackDesign(design);

  assert(time_.GetSize() == currentIteration);
  time_.Push_back(cfs_timer_->GetWallTime());
  LOG_DBG(opt) << "CI: ci=" << currentIteration << " wt=" << cfs_timer_->GetWallTime() << " -> " << time_.ToString();

  // eventually set special result
  EvaluateSpecialResults();

  // also log to info node, append the iteration
  PtrParamNode iteration = optInfoNode->Get(ParamNode::PROCESS)->Get("iteration", ParamNode::APPEND);

  // write the header only once - we might keep the iteration number
  if(log.file && objectives.GetHistorySize() == 1)
  {
    for(auto& prop: aux_log)
      log.AddToHeader(prop.first);
    *log.file << log.fileHeader << endl;
  }

  // write the current logging information
  LogFileLine(log.file, iteration); // also ParamNode is to be written
  baseOptimizer_->LogFileLine(log.file, iteration);
  for(auto& prop: aux_log)
  {
    if(log.file)
      *log.file << "\t " << prop.second;
    iteration->Get(prop.first)->SetValue(prop.second);
  }
  if(log.file)
    *log.file << endl;

  // this late to have observations evaluated
  constraints.PushBackHistory();

  // option gradplot for some FeaturedDesign, other (e.g. SIMP) ignore
  design->WriteGradientFile();

  // option hessexport for FeatureMappingDesign, others ignore
  design->WriteHessExportFile();

  // this writes the most current solved forward problem via the driver to gid or whatever
  // keep "commitStride == 1 || " for readability!
  bool store = currentIteration == 0 || commitStride == 1 || ((commitStride > 0) && currentIteration % commitStride == 0);
  LOG_TRACE2(opt) << "CI: " << currentIteration << " objective=" << objectives.GetHistoryValue() << " store=" << store;

  if(store)
  {
    StoreResults();
    lastStoredResult_ = currentIteration;
    // see FinalizeStoreResults() !
  }
  else
  {
    for(unsigned int e = 0; e < me->excitations.GetSize(); e++)
    {
      Excitation& ex = me->excitations[e];
      Context& ctxt = manager.GetContext(&ex);
      ctxt.GetDriver()->GetResultHandler()->streamOnly = true;
    }

    StoreResults();

    for(unsigned int e = 0; e < me->excitations.GetSize(); e++)
    {
      Excitation& ex = me->excitations[e];
      Context& ctxt = manager.GetContext(&ex);
      ctxt.GetDriver()->GetResultHandler()->streamOnly = false;
    }
  }

  // IPOPT does own logging -> otherwise show the user we are alive
  string f = GetIterationFrequency();
  if(optimizer_ != IPOPT_SOLVER && optimizer_ != SNOPT_SOLVER)
  {
    cout << "iteration " << (currentIteration);
    if(f != "") cout << " f = " << f << " Hz";
    cout << " -> cost = " << objectives.GetHistoryValue() << endl;
  }

  // if a python function is registered, call it.
  python->CallHook(PythonKernel::OPT_POST_ITER);

  // update possible tunes (they self register)
  for(Tune* tune : tunes)
    tune->Update(currentIteration);

  currentIteration++;
  problemWithinIteration = 0;

  return iteration;
}

/** call not later than PostInit2(). Add the property to file and info.xml iteration log. */
void Optimization::RegisterAuxLogValue(const std::string& name, const std::string& initial)
{
  for(auto& n : aux_log)
    if(n.first == name)
      throw Exception("name '" + name + "' already registered in Optimization::aux_log");

  aux_log.Push_back(std::make_pair(name, initial));
}

/** update the value to be used by next CommitIteration(). The name should be already set by RegisterAuxLogValue() */
void Optimization::SetAuxLogValue(const std::string& name, const std::string& value)
{
  for(auto& n : aux_log)
    if(n.first == name)
    {
      n.second = value;
      return;
    }

  throw Exception("name '" + name + "' not registered in Optimization::aux_log");
}

void Optimization::LogFileLine(ofstream* out, PtrParamNode iteration)
{
  double duration = time_.Last() - (time_.GetSize() > 1 ? time_[time_.GetSize() - 2] : 0.0);

  if(out)
  {
    *out << currentIteration;
    *out << std::defaultfloat << std::setprecision(6); // uses scientific only when needed

    if(isMultiObjective_)
      *out << " \t" << baseOptimizer_->GetObjectiveValue();

    if(context->IsHarmonic())
      *out << " \t" << GetIterationFrequency();

    for(unsigned int i = 0; i < objectives.data.GetSize(); i++)
    {
      Objective* c = objectives.data[i];
      *out << " \t" << c->GetOrgValue();
      if(c->GetType() == Function::BANDGAP)
      {
        // we search with the wave vectors for minimun and maximum
        *out << " \t" << c->bandgap.lower.col;
        *out << " \t" << c->bandgap.upper.col;
      }
      if(c->GetTerm() == Objective::PENALTY)
        *out << " \t" << c->CalcPenalty(c->GetOrgValue());
      if(c->GetTerm() == Objective::LN1P)
        *out << " \t" << c->CalcLn1p(c->GetOrgValue());
      if(c->GetTerm() == Objective::POWER)
        *out << " \t" << std::pow(c->GetOrgValue(),c->GetParameter());
    }

    // more details are written to .info.xml in ErsatzMaterial::LogFileLine()
    if(me->IsEnabled() && me->excitations.GetSize() > 1)
      for(Excitation& ex : me->excitations)
        *out << " \t" << ex.cost;

    *out << " \t" << duration;
    if(design->HasAlphaVariable())
      *out << " \t" << design->GetAlphaVariable();
    if(design->HasSlackVariable() && !objectives.Has(Function::SLACK))
      *out << " \t" << design->GetSlackVariable();
  }

  iteration->Get("number")->SetValue(currentIteration);

  if(context->IsHarmonic())
    iteration->Get("frequency")->SetValue(GetIterationFrequency());

  iteration->Get("duration")->SetValue(duration);

  if(design->HasAlphaVariable()) // needs to be written to the plot.dat file in ErsatzMaterial as Optimization::Optimization() knows no design yet
    iteration->Get("alpha")->SetValue(design->GetAlphaVariable());

  if(design->HasSlackVariable() && !objectives.Has(Function::SLACK))
    iteration->Get("slack")->SetValue(design->GetSlackVariable());

  if(isMultiObjective_)
  {
    std::stringstream ss;
    ss << std::setprecision(15) << baseOptimizer_->GetObjectiveValue();
    iteration->Get("multiObjectiveValue")->SetValue(ss.str());
  }

  for(unsigned int i = 0; i < objectives.data.GetSize(); i++)
  {
    Objective* c = objectives.data[i];

    // set the precision for the output of objective function values
    std::stringstream ss;
    ss << std::setprecision(15) << c->GetOrgValue();

    iteration->Get(c->ToString())->SetValue(ss.str());
    if(c->GetType() == Function::BANDGAP)
    {
      // we search with the wave vectors for minimun and maximum
      iteration->Get("max_ef_" + to_string(c->bandgap.lower_ev) + "_wv")->SetValue(c->bandgap.lower.col);
      iteration->Get("min_ef_" + to_string(c->bandgap.upper_ev) + "_wv")->SetValue(c->bandgap.upper.col);
    }
    if(c->GetTerm() == Objective::PENALTY)
      iteration->Get(c->ToString() + "_penalty")->SetValue(c->CalcPenalty(c->GetOrgValue()));
    if(c->GetTerm() == Objective::LN1P)
      iteration->Get(c->ToString() + "_ln1p")->SetValue(c->CalcLn1p(c->GetOrgValue()));
    if(c->GetTerm() == Objective::POWER)
      iteration->Get(c->ToString() + "_power")->SetValue(std::pow(c->GetOrgValue(),c->GetParameter()));
  }

  // we might have bloch information calculated int ErsatzMaterial::CommitIteration()
  for(unsigned int i = 0; i < log.bloch_info.GetSize(); i++)
    *out << " \t" << std::get<1>(log.bloch_info[i]);

  // For iteration 0 we want also the constraint values but they were not evaluated.
  // For any iteration we need to evaluate the observe constraints
  // A problem are the slope constraints, they need to be evaluated in local mode
  // and Done() forms the global result
  for(unsigned int i = 0, m = constraints.view->GetNumberOfTotalConstraints(); i < m; i++)
  {
    Condition* g = constraints.view->Get(i); // traverse in local mode
    if(g->GetValue() == -1.0 || g->IsObservation())
      CalcConstraint(g);
  }
  constraints.view->Done();
  double max = -1.;
  for(unsigned int i = 0; i < constraints.all.GetSize(); i++)
  {
    Condition* g = constraints.all[i]; // Now traverse in global mode
    if(g->GetType() == Function::SHAPE_INF)
      continue; //TODO: MaxValue does not correctly set indexes in view

    // Calculate the max value of multiple displacement constraints
    if((g->GetType() == Function::OUTPUT || g->GetType() == Function::SQUARED_OUTPUT) && g->output_multiple_nodes > 0) {
      max = std::max(max,std::abs(g->GetValue()));
    }
    if(g->IsLocalCondition())
    {
      LocalCondition* local = dynamic_cast<LocalCondition*>(g);
      double minmax  = local->CalcMinMaxAbsValue();
      int    inf_cnt = local->CountInfeasibles();
      double mean    = progOpts->DoDetailedInfo() ? local->CalcMeanAbsValue() : -1.0;
      if(out) {
        *out << " \t" << minmax;
        if(progOpts->DoDetailedInfo())
          *out << " \t" << mean;
        *out << " \t" << inf_cnt;
      }

      iteration->Get((g->GetBound() != g->LOWER_BOUND ? "max_abs_" : "min_abs_") + g->ToString())->SetValue(minmax);
      if(progOpts->DoDetailedInfo())
        iteration->Get("mean_abs_" + g->ToString())->SetValue(mean);
      iteration->Get("infeas_count_" + g->ToString())->SetValue(inf_cnt);
    }

    else
    {
      double value = g->GetValue();
      if(g->delta_logging)
        value = value - g->GetBoundValue();
      if(out && (g->GetType() != Function::EIGENFREQUENCY || log.plot_ev)) // don't spoil
        *out << " \t" << value;
      // excitation sensitive constraints are printed in the excitation list if there is one (ErsatzMaterial::CommitIteration())
      if(!g->IsExcitationSensitive() || g->ctxt->excitations.GetSize() < 2) {
        // set the precision for the output of objective function values
        std::stringstream ss;
        ss << std::setprecision(15) << value;
        iteration->Get(g->ToString())->SetValue(ss.str());
      }
    }
  }
  // max output_constraint value
  if (out && max > -1.0)
    *out << " \t output_max = " << max;

  if(out && log.design)
  {
    StdVector<double> d;
    d.Resize(design->GetNumberOfVariables());
    design->WriteDesignToExtern(d.GetPointer(), false);
    for(unsigned int i = 0; i < design->GetNumberOfVariables(); ++i){
      *out << " \t" << d[i];
    }
  }

  if(out && log.designGradient)
  {
    for(unsigned int i = 0; i < objectives.data.GetSize(); ++i)
    {
      Objective* f = objectives.data[i];
      StdVector<double> d;
      d.Resize(design->GetNumberOfVariables());
      d.window.Set(d);
      design->WriteGradientToExtern(d, DesignElement::COST_GRADIENT, DesignElement::PLAIN, f, false);
      for(unsigned int j = 0; j < design->GetNumberOfVariables(); ++j)
        *out << " \t" << setprecision(15) << d[j];
    }
  }

  if(out && log.designConstraintGradients)
  {
    for(unsigned int i = 0; i < constraints.all.GetSize(); ++i)
    {
      Condition* g = constraints.all[i]; // Now traverse in global mode
      if(g->GetType() == Function::SHAPE_INF) continue; //TODO: MaxValue does not correctly set indexes in view

      if(g->IsLocalCondition()) continue; // this would be huge
      else
      {
        StdVector<double> d;
        d.Resize(design->GetNumberOfVariables());
        d.window.Set(d);
        design->WriteGradientToExtern(d, DesignElement::CONSTRAINT_GRADIENT, DesignElement::PLAIN, g, false);
        for(unsigned int j = 0; j < design->GetNumberOfVariables(); ++j)
          *out << " \t" << setprecision(15) << d[j];
      }
    }
  }
  if(out)
    out->flush();
}



DesignElement::Type Optimization::ToDesign(const SinglePDE* pde) const
{
  if(pde->GetName() == "electrostatic") return DesignElement::POLARIZATION;
  if(pde->GetName() == "LatticeBoltzmann") return DesignElement::DENSITY;
  if(pde->GetName() == "mechanic") return DesignElement::DENSITY;
  if(pde->GetName() == "acoustic") return DesignElement::DENSITY;

  throw Exception("invalid");
}

//App::Type Optimization::ToApp(DesignElement::Type dt)
//{
//  switch(dt)
//  {
//  case DesignElement::DENSITY:
//      return App::MECH; // wrong for buckling
//  case DesignElement::DENSITY:
//    return App::ACOUSTIC;
//  case DesignElement::POLARIZATION:
//    return App::ELEC;
//  default:
//    EXCEPTION("DesignType " << DesignElement::type.ToString(dt) << " doesn't map to App::Type");
//  }
//}


Optimization::Log::Log()
{
  this->columns_ = 0;
  this->design = false;
  this->designGradient = false;
  this->designConstraintGradients = false;
  this->gradNorm = progOpts->DoDetailedInfo();
  this->plot_ev = true;
  this->file = nullptr;
  this->fileHeader = "";
}

void Optimization::Log::Init(Optimization* opt, const string& log_name, PtrParamNode pn_log)
{
  if(log_name != "false")
  {
    string name = log_name == "[problem]" || log_name == "true" ? (progOpts->GetSimName() + ".plot.dat") : log_name;
    file = new ofstream(name.c_str());
    if(file == NULL)
      throw Exception("cannot open log file " + name + " for writing");

    if(pn_log != NULL)
    {
      design = pn_log->Get("design")->As<bool>();
      designGradient = pn_log->Get("designGradient")->As<bool>();
      designConstraintGradients = pn_log->Get("designConstraintGradients")->As<bool>();
    }
  }

  StdVector<Condition*> ev = opt->constraints.GetList(Condition::EIGENFREQUENCY);
  if(!ev.IsEmpty() && ev.First()->GetExcitation()->DoBloch() && ev.First()->DoFullBloch())
  {
    plot_ev = progOpts->DoDetailedInfo();

    // see ErsatzMaterial::CommitInteration()
    bloch_info.Push_back(std::make_tuple("bandgap", -1.0));

    StdVector<string> found;
    for(unsigned int i = 0; i < ev.GetSize(); i++) {
      std::string key = "ev_" + std::to_string(ev[i]->GetEigenValueID()) + (ev[i]->GetBound() == Condition::LOWER_BOUND ? "_min" : "_max");
      // we have wave vector times each ev, add only one!
      if(!found.Contains(key)) {
        bloch_info.Push_back(std::make_tuple(key, -1.0));
        found.Push_back(key);
      }
    }
  }
}

 Optimization::Log::~Log()
 {
   // if write to file close it
   if(file != NULL)
   {
     file->close();
     delete file;
     file = NULL;
   }
 }

void Optimization::Log::AddToHeader(const string& label)
{
  fileHeader += columns_ == 0 ? "#" : "\t";

  columns_++;

  fileHeader += boost::lexical_cast<string>(columns_) + ":" + label;
}

