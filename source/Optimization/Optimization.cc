#include <assert.h>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/programOptions.hh"
#include "DataInOut/resultHandler.hh"
#include "Domain/domain.hh"
#include "Driver/assemble.hh"
#include "Driver/baseSolveStep.hh"
#include "Driver/basedriver.hh"
#include "Driver/formsContext.hh"
#include "Driver/harmonicDriver.hh"
#include "Forms/linearForm.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignStructure.hh"
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
#include "Optimization/ParamMat.hh"
#include "Optimization/PiezoSIMP.hh"
#include "Optimization/PiezoParamMat.hh"
#include "Optimization/SIMP.hh"
#include "Optimization/ShapeGrad.hh"
#include "Optimization/ShapeOpt.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/basePDE.hh"
#include "Utils/tools.hh"
#include "boost/filesystem.hpp"
#include "def_use_ipopt.hh"
#include "def_use_knitro.hh"
#include "def_use_scpip.hh"
#include "def_use_snopt.hh"
#include "def_use_feas_scp.hh"

// IPOPT, SCPIP and SnOpt are not necessarily linked
#ifdef USE_IPOPT
  #include "Optimization/Optimizer/IPOPTHolder.hh"
  #include "Optimization/Optimizer/FeasPP.hh"
#endif
#ifdef USE_SCPIP
  #include "Optimization/Optimizer/SCPIP.hh"
#endif
#ifdef USE_FEAS_SCP
  #include "Optimization/Optimizer/FeasSCP.hh"
#endif
#ifdef USE_SNOPT
  #include "Optimization/Optimizer/SnOpt.hh"
#endif
#ifdef USE_KNITRO
  #include "Optimization/Optimizer/KNITRO.hh"
#endif

using namespace CoupledField;
using namespace std;
namespace fs = boost::filesystem;



DECLARE_LOG(opt)
DEFINE_LOG(opt, "opt")

// instantiation of the static elements
Enum<Optimization::Optimizer>        Optimization::optimizer;
Enum<Optimization::Application>      Optimization::application;
Enum<Optimization::CommitMode>       Optimization::commitMode;



Optimization::Optimization()
{
  this->lastStoredResult_ = -1;
  this->design = NULL;
  this->baseOptimizer_ = NULL;
  this->harmonic = BasePDE::IsComplex(domain->GetDriver()->GetAnalysisType()); // no multi-sequence
  this->currentIteration = 0; // a 1 or 0 can make a lot of difference! 0 is initial design!
  this->writeCounter_ = 0;
  this->problemSolvedCounter = 0;
  this->problemWithinIteration = 0;
  this->grid = domain->GetGrid();

  // inject the driver and tell him that we do optimization
  BaseDriver* driver = domain->GetDriver();
  if(driver->GetDriverClass() != BaseDriver::SINGLE_DRIVER)
    throw Exception("optimization not implemented for driver " + driver->GetDriverClass());

  optInfoNode = info->Get("optimization");   // store our info results here
  PtrParamNode pn = param->Get("optimization"); // read our parameters from the xml file
  
  // in transient optimization one can specify the initial value as a solution to a static problem and a weight for it (just in tracking)
  firstStepStatic = pn->Has("firstStepStatic");
  if(firstStepStatic){
    otherStepWeight = 1.0 - pn->Get("firstStepStatic")->Get("weight")->As<Double>();
  }else{
    otherStepWeight = 1.0;
  }

  // the tool to solve the optimization problem
  optimizer_ = optimizer.Parse(pn->Get("optimizer")->Get("type")->As<string>());
  maxIterations = pn->Get("optimizer")->Get("maxIterations")->As<Integer>();

  // might read a multiObjective problem
  objectives.Read(pn->Get("costFunction"));
  objectives.ToInfo(optInfoNode->Get(ParamNode::HEADER)->Get("objective"));

  SetPDEs(OptimizationMaterial::system.Parse(pn->Get("ersatzMaterial")->Get("material")->As<string>()));

  this->assemble_ = pde->getPDE_assemble();

  // constraints to be added later -- it is so much easier with the ParamNodes
  this->log.fileHeader = harmonic ? "#iter\tfreq" : "#iter";
  for(unsigned int i = 0; i < objectives.data.GetSize(); i++)
    this->log.fileHeader += "\t" + objectives.data[i]->GetName();
  this->log.fileHeader += "\tchange\tproblems";

  // multiple excitations are are toggled via attribute. Only if enabled we read the optional element
  // actually part of costFunction - but we store in Optimization itself!
  bool dme = pn->Get("costFunction")->Get("multiple_excitation")->As<bool>();
  this->me = new MultipleExcitation(dme, dme ? pn->Get("costFunction")->Get("multipleExcitation", ParamNode::PASS) : PtrParamNode());
  if(dme) this->me->ToInfo(optInfoNode->Get(ParamNode::HEADER)->Get("multipleExcitations"));

  // slope constraints to be processed in SIMP -> Constraints::PostProc
  ParamNodeList list = pn->GetList("constraint");
  constraints.Read(list);
  PtrParamNode in = optInfoNode->Get(ParamNode::HEADER)->Get("constraints");

  // constraints.ToInfo() is called in PostInitSecond()

  for(unsigned int i = 0; i < constraints.all.GetSize(); i++)
  {
    Condition* g = constraints.all[i];
    if(!g->IsLocalCondition())
      this->log.fileHeader += "\t" + g->ToString();
    else
    {
      this->log.fileHeader += "\tmax_" + g->ToString();
      this->log.fileHeader += "\tmean_" + g->ToString();
    }
  }

  log.Init(pn->Get("log")->As<string>(), pn->Get("logging", ParamNode::PASS)); // is fail save

  // the commit stuff
  string cm = pn->Has("commit") ? pn->Get("commit")->Get("mode")->As<string>() : "forward";
  this->commitMode_ = commitMode.Parse(cm);
  this->commitStride = pn->Has("commit") ? pn->Get("commit")->Get("stride")->As<Integer>() : 1;
  optInfoNode->Get("commit")->Get("mode")->SetValue(cm);
  optInfoNode->Get("commit")->Get("stride")->SetValue(commitStride);
  
  // write out the directory where the HALTOPT file will be searched for
  optInfoNode->Get("haltopt_directory")->SetValue(fs::current_path().string());

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


void Optimization::PostInitSecond()
{
  PtrParamNode opt = param->Get("optimization")->Get("optimizer");

  switch(optimizer_)
  {
    case IPOPT_SOLVER:
         #ifdef USE_IPOPT
           baseOptimizer_ = new IPOPTHolder(this, opt);
         #else
           throw Exception("CFS++ was compiled w/o IPOPT!");
         #endif
         break;

    case SCPIP_SOLVER:
         #ifdef USE_SCPIP
           #ifdef USE_FEAS_SCP
             throw Exception("CFS++ was compiled with FeasSCP - this disables SCPIP");
           #endif
           baseOptimizer_ = new SCPIP(this, opt);
         #else
           throw Exception("CFS++ was compiled w/o SCPIP");
         #endif
         break;

    case FEAS_PP_SOLVER:
         #ifndef USE_IPOPT
           throw Exception("CFS++ needs to be compiled with IPOPT to use feasPP");
         #else
           baseOptimizer_ = new FeasPP(this, opt);
         #endif
         break;

    case FEAS_SCP_SOLVER:
         #ifdef USE_FEAS_SCP
           baseOptimizer_ = new FeasSCP(this, opt);
         #else
           throw Exception("CFS++ was compiled w/o FeasSCP");
         #endif
         break;
         
    case SNOPT_SOLVER:
         #ifdef USE_SNOPT
           baseOptimizer_ = new SnOpt(this, opt);
         #else
           throw Exception("CFS++ was compiled w/o SnOpt");
         #endif
      break;

    case KNITRO_SOLVER:
         #ifdef USE_KNITRO
           baseOptimizer_ = new KNITRO(this, opt);
         #else
           throw Exception("CFS++ was compiled w/o KNITRO");
         #endif
      break;

    case OPTIMALITY_CONDITION:
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

  constraints.ToInfo(optInfoNode->Get(ParamNode::HEADER)->Get("constraints"), GetMultipleExcitation());

  unsigned int n = design->GetNumberOfVariables();
  if (this->log.design) {
    for (unsigned int i = 0; i < n; ++i) {
      this->log.fileHeader += "\t";
      this->log.fileHeader += "design";
    }
  }
  if (this->log.designGradient) {
    for (unsigned int i = 0; i < n; ++i) {
      this->log.fileHeader += "\t";
      this->log.fileHeader += "designGradient";
    }
  }
  design->SetOptimizer(baseOptimizer_);
  // add plot logging of the optimizer
  this->log.fileHeader += baseOptimizer_->LogFileHeader();
}


void Optimization::SetEnums()
{
  Function::type.SetName("Function::Type");
  Function::type.Add(Function::MULTI_OBJECTIVE, "multiObjective");
  Function::type.Add(Function::COMPLIANCE, "compliance");
  Function::type.Add(Function::OUTPUT, "output");
  Function::type.Add(Function::DYNAMIC_OUTPUT, "dynamicOutput");
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
  Function::type.Add(Function::HOM_TENSOR, "homTensor");
  Function::type.Add(Function::MAXWELL_HOM_TENSOR, "maxwellHomTensor");
  Function::type.Add(Function::HOM_TRACKING, "homTracking");
  Function::type.Add(Function::HOM_FROBENIUS_PRODUCT, "homFrobeniusProduct");
  Function::type.Add(Function::MAXWELL_HOM_TRACKING, "maxwellHomTracking");
  Function::type.Add(Function::BITENSOR, "bitensor");
  Function::type.Add(Function::POISSONS_RATIO, "poissonsRatio");
  Function::type.Add(Function::YOUNGS_MODULUS, "youngsModulus");
  Function::type.Add(Function::YOUNGS_MODULUS_E1, "youngsModulusE1");
  Function::type.Add(Function::YOUNGS_MODULUS_E2, "youngsModulusE2");
  Function::type.Add(Function::TYCHONOFF, "tychonoff");
  Function::type.Add(Function::TEMPERATURE, "temperature");
  Function::type.Add(Function::PROJECTION, "projection");
  Function::type.Add(Function::GREYNESS, "greyness");
  Function::type.Add(Function::STRESS, "stress");
  Function::type.Add(Function::STRESS_DENSITY, "stressDensity");
  Function::type.Add(Function::ISOTROPY, "isotropy");
  Function::type.Add(Function::MAXWELL_ISOTROPY, "maxwellIsotropy");
  Function::type.Add(Function::BIISOTROPY, "biIsotropy");
  Function::type.Add(Function::ISO_ORTHOTROPY, "iso-orthotropy");
  Function::type.Add(Function::ORTHOTROPY, "orthotropy");
  Function::type.Add(Function::SLOPE, "slope");
  Function::type.Add(Function::GLOBAL_SLOPE, "globalSlope");
  Function::type.Add(Function::MOLE, "mole");
  Function::type.Add(Function::GLOBAL_MOLE, "globalMole");
  Function::type.Add(Function::OSCILLATION, "oscillation");
  Function::type.Add(Function::GLOBAL_OSCILLATION, "globalOscillation");
  Function::type.Add(Function::JUMP, "jump");
  Function::type.Add(Function::GLOBAL_JUMP, "globalJump");
  Function::type.Add(Function::BUMP, "bump");
  Function::type.Add(Function::DESIGN_TRACKING, "designTracking");
  Function::type.Add(Function::SUM_MODULI, "sumModuli");
  Function::type.Add(Function::GLOBAL_SUM_MODULI, "globalSumModuli");
  Function::type.Add(Function::LAMINATES_VOL, "laminatesVolume");
  Function::type.Add(Function::GLOBAL_LAMINATES_VOL, "globalLaminatesVolume");
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
  Function::type.Add(Function::DETERMINANT_MATRIX, "determinantMatrix");
  Function::type.Add(Function::ROTATIONAL_MATRIX_1, "rotationalMatrix1");
  Function::type.Add(Function::ROTATIONAL_MATRIX_2, "rotationalMatrix2");
  Function::type.Add(Function::DESIGN_BOUND, "designBound");
  Function::type.Add(Function::MULTIMATERIAL_SUM, "multimaterial_sum");
  Function::type.Add(Function::SLACK, "slack");

  Function::Local::locality.SetName("Function::Local::Locality");
  Function::Local::locality.Add(Function::Local::DEFAULT, "default");
  Function::Local::locality.Add(Function::Local::NEXT, "next");
  Function::Local::locality.Add(Function::Local::NEXT_AND_REVERSE, "next_and_reverse");
  Function::Local::locality.Add(Function::Local::PREV_NEXT, "prev_next");
  Function::Local::locality.Add(Function::Local::PREV_NEXT_AND_REVERSE, "prev_next_and_reverse");
  Function::Local::locality.Add(Function::Local::DEG_45_STAR, "45_deg_star");
  Function::Local::locality.Add(Function::Local::DEG_45_STAR_AND_REVERSE, "45_deg_star_and_reverse");
  Function::Local::locality.Add(Function::Local::BOUNDARY, "boundary");
  Function::Local::locality.Add(Function::Local::ELEMENT, "element");
  Function::Local::locality.Add(Function::Local::MULT_DESIGNS_ELEMENT, "multiple_designs_element");
  Function::Local::locality.Add(Function::Local::MULT_DESIGNS_NEXT, "multiple_designs_next");
  Function::Local::locality.Add(Function::Local::MULT_DESIGNS_NEXT_AND_REVERSE, "multiple_designs_next_and_reverse");
  Function::Local::locality.Add(Function::Local::MULT_DESIGNS_PREV_NEXT, "multiple_designs_prev_next");
  Function::Local::locality.Add(Function::Local::MULT_DESIGNS_PREV_NEXT_AND_REVERSE, "multiple_designs_prev_next_and_reverse");

  Function::Local::phase.SetName("Function::Local::Phase");
  Function::Local::phase.Add(Function::Local::BOTH, "both");
  Function::Local::phase.Add(Function::Local::VOID, "void");
  Function::Local::phase.Add(Function::Local::MATERIAL, "material");

  Function::stressType.SetName("Function::StressType");
  Function::stressType.Add(Function::MECH, "mech");
  Function::stressType.Add(Function::PIEZO, "piezo");
  Function::stressType.Add(Function::ONLY_COUPLING, "only_coupling");

  Condition::bound.SetName("Condition::Bound");
  Condition::bound.Add(Condition::EQUAL, "equal");
  Condition::bound.Add(Condition::LOWER_BOUND, "lowerBound");
  Condition::bound.Add(Condition::UPPER_BOUND, "upperBound");


  ObjectiveContainer::StoppingRule::type.SetName("ObjectiveContainer::StoppingRule::Type");
  ObjectiveContainer::StoppingRule::type.Add(ObjectiveContainer::StoppingRule::DESIGN_CHANGE, "designChange");
  ObjectiveContainer::StoppingRule::type.Add(ObjectiveContainer::StoppingRule::REL_COST_CHANGE, "relativeCostChange");

  DesignStructure::filterSpace.SetName("DesignStructure::FilterSpace");
  DesignStructure::filterSpace.Add(DesignStructure::RADIUS, "radius");
  DesignStructure::filterSpace.Add(DesignStructure::VOLUME_RADIUS, "volumeRadius");
  DesignStructure::filterSpace.Add(DesignStructure::MAX_EDGE, "maxEdge");

  optimizer.SetName("Optimization::Optimizer");
  optimizer.Add(OPTIMALITY_CONDITION, "optimalityCondition");
  optimizer.Add(IPOPT_SOLVER, "ipopt");
  optimizer.Add(SCPIP_SOLVER, "scpip");
  optimizer.Add(FEAS_SCP_SOLVER, "feasSCP");
  optimizer.Add(FEAS_PP_SOLVER, "feasPP");
  optimizer.Add(SNOPT_SOLVER, "snopt");
  optimizer.Add(KNITRO_SOLVER, "knitro");
  optimizer.Add(SHAPE_SOLVER, "shapeOpt");
  optimizer.Add(EVALUATE_INITIAL_DESIGN, "evaluate");
  optimizer.Add(GRADIENT_CHECK, "gradientCheck");

  ErsatzMaterial::method.SetName("ErsatzMaterial::Method");
  ErsatzMaterial::method.Add(ErsatzMaterial::SIMP_METHOD, "simp");
  ErsatzMaterial::method.Add(ErsatzMaterial::PARAM_MAT, "paramMat");
  ErsatzMaterial::method.Add(ErsatzMaterial::SHAPE_GRAD, "shapeGrad");
  ErsatzMaterial::method.Add(ErsatzMaterial::SHAPE_OPT, "shapeOpt");
  ErsatzMaterial::method.Add(ErsatzMaterial::SHAPE_PARAM_MAT, "shapeParamMat");
  
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
  OptimizationMaterial::system.Add(OptimizationMaterial::ACOUSTIC, "acoustic");
  OptimizationMaterial::system.Add(OptimizationMaterial::ELEC, "maxwellHom");

  application.SetName("Optimization::Application");
  application.Add(NO_APP, "no_app");
  application.Add(ACOUSTIC, "acoustic");
  application.Add(HEAT, "heat");
  application.Add(LAPLACE, "laplace");
  application.Add(MECH, "mech");
  application.Add(MASS, "mass");
  application.Add(ELEC, "elec");
  application.Add(PIEZO_COUPLING, "piezoCoupling");
  application.Add(PRESSURE, "pressure");
  application.Add(CHARGE_DENSITY, "chargeDensity");
  application.Add(STRESS, "stress");

  LevelSet::Action::type.SetName("LevelSet::Action::Type");
  LevelSet::Action::type.Add(LevelSet::Action::SIGNED_DISTANCE_FIELD, "signedDistanceField");
  LevelSet::Action::type.Add(LevelSet::Action::TRIVIAL_HOLE, "trivialHole");
  LevelSet::Action::type.Add(LevelSet::Action::DO_SHAPE_STEP, "shapeStep");

  MultipleExcitation::type.SetName("MultipleExcitation::Type");
  MultipleExcitation::type.Add(MultipleExcitation::NO_TYPE, "no_type");
  MultipleExcitation::type.Add(MultipleExcitation::FIXED_WEIGHT, "fixed_weights");
  MultipleExcitation::type.Add(MultipleExcitation::META_OBJECTIVE, "meta_objective");
  MultipleExcitation::type.Add(MultipleExcitation::HOMOGENIZATION_TEST_STRAINS, "homogenizationTestStrains");
  MultipleExcitation::type.Add(MultipleExcitation::MAXWELL_HOMOGENIZATION_TEST_STRAINS, "maxwellHomogenizationTestStrains");
}

bool Optimization::IsTransient() {
  return(domain->GetDriver()->GetAnalysisType() == BasePDE::TRANSIENT);
}

double Optimization::GetStepWeight(unsigned int ts) const{
  unsigned int nts = domain->GetDriver()->GetNumSteps();
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

bool Optimization::DoStopOptimization()
{
  PtrParamNode in = optInfoNode->Get(ParamNode::SUMMARY)->Get("break");
  // check if the HALTOPT file exists
  if(fs::exists("HALTOPT"))
  {
    bool good = fs::remove("HALTOPT");
    if(!good) throw new Exception("Could not remove file 'HALTOPT' after detection");
    in->Get("converged")->SetValue("no");
    in->Get("reason/msg")->SetValue("Detected file 'HALTOPT'");
    return true;
  }
  
  ObjectiveContainer::StoppingRule& stop = objectives.stop;

  // we need a minimum number of iterations to be sure we are in a minimum
  unsigned int hs = objectives.GetHistorySize();
  if(hs <= stop.queue) return false;

  if(stop.GetType() == ObjectiveContainer::StoppingRule::REL_COST_CHANGE)
  {
    for(unsigned int i = hs-1; i >= (hs - stop.queue); i--)
    {
      double delta = objectives.GetHistoryValue(true, i) - objectives.GetHistoryValue(true, i-1);
      double rel = abs(delta / objectives.GetHistoryValue(true, i));
      if(rel > stop.value) return false;
    }
  }
  else // ObjectiveContainer::StoppingRule::DESIGN_CHANGE
  {
    for(unsigned int n = objectives.design_change.GetSize() - 1, i = n; i >= max((unsigned int) 0, n - stop.queue); i--)
      if(objectives.design_change[i] > stop.value)
        return false;
  }

  // the relative values for the whole queue are smaller than the requirement -> we are done! :)
  in->Get("converged")->SetValue("practically");

  if(stop.GetType() == ObjectiveContainer::StoppingRule::REL_COST_CHANGE)
    in->Get("reason/msg")->SetValue("Too small relative change in objective function");
  else
    in->Get("reason/msg")->SetValue("Too small change in design");

  in->Get("reason/value")->SetValue(stop.value);
  in->Get("reason/queue")->SetValue(stop.queue);
  return true;
}


/** read only the very basic stuff */
Optimization* Optimization::CreateInstance()
{
  // set the enums we need
  Optimization::SetEnums(); // sets also ErsatzMaterial::Method
  DesignElement::SetEnums();
  DesignMaterial::SetEnums();

  if(!param->Has("optimization")) return NULL;

  // we assume ersatz material, currently there is nothing else.
  // note, we read method again in the ersatz material constructor.
  PtrParamNode em = param->Get("optimization")->Get("ersatzMaterial");

  ErsatzMaterial::Method method = ErsatzMaterial::method.Parse(em->Get("method")->As<string>());
  OptimizationMaterial::System material = OptimizationMaterial::system.Parse(em->Get("material")->As<string>());
  
  Optimization* opt = NULL;
  
  switch(method)
  {
  case ErsatzMaterial::SIMP_METHOD:
    switch(material)
    {
    case OptimizationMaterial::MECH:
    case OptimizationMaterial::ACOUSTIC:
    case OptimizationMaterial::HEAT:
    case OptimizationMaterial::ELEC:
      opt = new SIMP(); // generally single PDE!
      break;
      
    case OptimizationMaterial::PIEZOCOUPLING:
      opt = new PiezoSIMP();
      break;
      
    default:
      assert(false);
      break;
    }
    break;
    
  // FMO, ShapeGrad, ...
  case ErsatzMaterial::PARAM_MAT:
    if(material == OptimizationMaterial::PIEZOCOUPLING)
      opt = new PiezoParamMat();
    else
      opt = new ParamMat();
    break;
  case ErsatzMaterial::SHAPE_OPT:
  case ErsatzMaterial::SHAPE_PARAM_MAT:
    opt = new ShapeOpt();
    break;
  case ErsatzMaterial::SHAPE_GRAD:
  {
    opt = new ShapeGrad();
    break;
  }
  default: throw Exception("Optimization not implemented");
  }
  
  // we have to do this, as PostInitSecond does already run CalcObjective/Gradient
  domain->SetOptimization(opt);
  
  // second initialization phase, constructs material
  opt->PostInit();
  // third initialization phase, constructs optimizer
  opt->PostInitSecond();

  return opt;
}

void Optimization::SolveProblem()
{
  // one driver is one multisequence step. We do this stuff here
  // and call the driver->StoreResults() multiple times f
  
  ResultHandler* rh = NULL;

  if(!IsTransient()){ // transient optimization saves results in a different way
    rh = domain->GetResultHandler();
    unsigned int mss = domain->GetDriver()->GetActSequenceStep();
    // max steps is high. The number is only relevant for hdf5, but there a hard limit
    rh->BeginMultiSequenceStep(mss, BasePDE::TRANSIENT, 9999);
  }

  Exception* e = NULL;
  try
  {
    PtrParamNode in = optInfoNode->Get(ParamNode::HEADER)->Get(optimizer.ToString(optimizer_));
    baseOptimizer_->ToInfo(in);
    baseOptimizer_->SolveOptimizationProblem();
    baseOptimizer_->ToInfo(in);
  }
  catch(Exception& ex)
  {
    e = new Exception(ex);
  }

  if(!IsTransient()){ // transient optimization saves results in a different way
    // do the finally - try to write results even if the optimizer broke down
    FinalizeStoreResults(); // when we have strides the results are written
    rh->FinishMultiSequenceStep();
    rh->Finalize();
  }
  if(e != NULL) throw *e;
  delete e;
}
PtrParamNode Optimization::CreateAdjointAnalysisIdNode(string child_name, Excitation* excite)
{
  BaseDriver* driver = domain->GetDriver();
  PtrParamNode base = driver->GetAnalysisId();
  string temp = "excite";
  PtrParamNode in = excite == NULL ? driver->CreateAnalysisIdChild(base, child_name) : driver->CreateAnalysisIdChild(base, "excite", excite->index, child_name);
  
  return in;
}


void Optimization::SolveStateProblem(Excitation* excite)
{
  BaseDriver* driver = domain->GetDriver();

  // this is the forward problem. Store the analysis_id in the driver such that
  // we can aquire it for the adjoint problem
  PtrParamNode analysis_id = excite == NULL ? driver->CreateAnalysisId("iter", currentIteration)
                                         : driver->CreateAnalysisId("iter", currentIteration, "excite", excite->index);
  
  if(IsTransient() && problemSolvedCounter > 0){ // transient optimization always has a mech pde
    SinglePDE* mech = domain->GetSinglePDE("mechanic");
    mech->ReReadResults();
    design->AppendOptimizationResults(mech);
    mech->GetSolveStep()->ReInit();
  }
                                         
  // Do not store the results. This is to be done in CommitIteration
  if(!harmonic || excite == NULL) 
    driver->SolveProblem(IsTransient(), analysis_id, NULL); // static and transient optimization
  else
  {
    LOG_DBG(opt) << "SSP: harmonic step=" << excite->f_link->step << " f=" << excite->f_link->freq;
    dynamic_cast<HarmonicDriver*>(driver)->ComputeFrequencyStep(excite->f_link->step, analysis_id);
  }

  problemSolvedCounter++;
  problemWithinIteration++;
}

void Optimization::SolveAdjointProblem(Excitation* excite, Function* f){
  // does almost the same as SolveStateProblem now, but passing, that we want the adjoint to be solved
  BaseDriver* driver = domain->GetDriver();
  
  AdjointParameters adjointParams(f, excite);
  
  if(IsTransient()){
    SinglePDE* mech = domain->GetSinglePDE("mechanic");
    mech->GetSolveStep()->ReInit();
  }

  // Do not store the results. This is adjoint.
  if(!harmonic) 
    driver->SolveProblem(false, CreateAdjointAnalysisIdNode("adjoint", excite), &adjointParams); // static and transient optimization
  else
    EXCEPTION("Harmonic adjoint not implemented!");
}

void Optimization::SolveAdjointProblems(Excitation* excite)
{
  // solve for objectives and constraints
  StdVector<Function*> ff = GetActiveFunctions();

  for(unsigned int i = 0; i < ff.GetSize(); ++i)
  {
    Function* f = ff[i];
    if(f->IsAdjointBased() && f->DoEvaluate(excite))
      SolveAdjointProblem(excite, f); // virtual! calls ErsatzMaterial implementation
  }
}

StdVector<Function*> Optimization::GetActiveFunctions() const
{
  StdVector<Function*> result;

  const unsigned int cn = objectives.data.GetSize();
  const unsigned int gn = constraints.active.GetSize();

  result.Resize(cn + gn);

  for(unsigned int i = 0; i < cn; i++)
  {
    result[i] = objectives.data[i];
    LOG_DBG2(opt) << "GAF: o=" << result[i]->ToString();
  }
  for(unsigned int i = 0; i < gn; i++)
  {
    result[cn + i] = constraints.active[i];
    LOG_DBG2(opt) << "GAF: g=" << result[cn + i]->ToString();
  }

  return result;
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

double Optimization::CalcObjective()
{
  // in objective.value_ we store the sum over all excitations w/o penalty but with normalization
  // in excitation.cost we store the sum over all objectives with penalty but w/o normalization

  // reset the objective values such that we can sum up normalized but unpenalized values
  for(unsigned int o = 0; o < objectives.data.GetSize(); o++)
    objectives.data[o]->ResetValue();

  double result = 0.0;

  // the multiple excitation case is a special case - for all other cases this is executed once
  for(unsigned int e = 0; e < me->excitations.GetSize(); e++)
  {
    Excitation& excite = me->excitations[e];
    excite.cost = 0.0;

    for(unsigned int o = 0; o < objectives.data.GetSize(); o++)
    {
      Objective* f = objectives.data[o];

      // some objectives are only to be evaluated for the last excitation
      if(!f->DoEvaluate(&excite)) continue;

      //double ov = CalcObjective(excite, f); // this is virtual!
      double ov = CalcFunction(excite, f, false); // this is virtual!
      excite.cost += ov * f->GetPenalty();

      // we ignore the weight if the evaluation happens only once! TODO why not omega*omega? - Fabian
      double weight = !f->DoEvaluateAlways() ? 1.0 : excite.normalized_weight;

      f->AddValue(ov * weight);

      result += ov * f->GetPenalty() * weight;
      LOG_DBG(opt) << "CalcObjective: ex=" << e << " obj=" << f->type.ToString(f->GetType()) << " ov=" << ov
          << " penalty" << f->GetPenalty() << " ex.cost=" << excite.cost << " nw=" << excite.normalized_weight
          << " wei=" << weight << " f->val=" << f->GetValue() << " result=" << result;
    }
  }

  return result;
}

void Optimization::CalcObjectiveGradient(StdVector<double>* grad_out)
{
  // reset the cost gradients in the design elements and sum them up in a weighted way
  // to perform multiple loads
  design->Reset(DesignElement::COST_GRADIENT);

  for(unsigned int obj = 0; obj < objectives.data.GetSize(); obj++)
  {
    Objective* cost = objectives.data[obj];
    // the multiple excitation case is a special case - for all other cases this is executed once
    for(unsigned int idx = 0; idx < me->excitations.GetSize(); idx++)
    {
      Excitation& excite = me->excitations[idx];
      // some objectives are only to be evaluated for the last excitation
      if(!cost->DoEvaluate(&excite)) continue;

      CalcFunction(excite, cost, true);
    }
  }

  if(grad_out != NULL)
    design->WriteGradientToExtern(*grad_out, DesignElement::COST_GRADIENT, DesignElement::SMART);
}



double Optimization::CalcConstraint(Condition* g)
{
  // assume when we have only one constraint which is not explicitly given, this is not the stress constraint!
  assert((g == NULL && constraints.active.GetSize() == 1 && constraints.active[0]->DoEvaluateAlways()) || g != NULL);

  if(g == NULL)
    g = constraints.active[0];

  double result = 0.0;

  for(unsigned int e = 0; e < me->excitations.GetSize(); e++)
  {
    Excitation& excite = me->excitations[e];
    // in the evaluate once case only the last excitation
    double v = g->DoEvaluate(&excite) ? CalcFunction(excite, g, false) : 0.0;
    double w = g->DoEvaluateAlways() ? excite.GetWeightedFactor(g) : 1.0;
    result += v * w;
    LOG_DBG(opt) << "CC ex=" << e << " eval=" << g->DoEvaluate(&excite) << " v=" << v << " alw=" << g->DoEvaluateAlways() << " w=" << w << " -> " << result;
  }

  g->SetValue(result);
  return result;
}

void Optimization::CalcConstraintGradient(Condition* g, StdVector<double>* grad_out)
{
  // assume when we have only one constraint which is not explicitly given, this is not the stress constraint!
  assert((g == NULL && constraints.active.GetSize() == 1 && !constraints.active[0]->DoEvaluateAlways()) || g != NULL);

  if(g == NULL)
    g = constraints.active[0];

  for(unsigned int i = 0; i < me->excitations.GetSize(); i++)
  {
    if(g->DoEvaluate(&me->excitations[i]))
      CalcFunction(me->excitations[i], g, true);
  }

  // copies from the design element gradient data to a memory array for external optimizers
  if(grad_out != NULL)
    design->WriteGradientToExtern(*grad_out, DesignElement::CONSTRAINT_GRADIENT, DesignElement::SMART, g);

  // if there is a <result ... value="constraintGradient" detail="penalizedVolume/*"
  if(g->special_result_idx != -1)
  {
    int base = design->FindDesign(g->GetDesignType());
    int n    = design->GetNumberOfElements();
    for(int i = n * base; i < n * (base + 1); i++) // TODO add access!
      design->data[i].specialResult[g->special_result_idx] = design->data[i].GetPlainGradient(NULL, g);
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

  // this will write the CFS result and history file
  if(!IsTransient())
  { // transient optimization saves results in a different way
    if(step_val == -1) {
      domain->GetDriver()->StoreResults(writeCounter_, currentIteration);
    } else {
      domain->GetDriver()->StoreResults(writeCounter_, step_val);
    }
    writeCounter_++;
  }
}

void Optimization::FinalizeStoreResults()
{
  // after the last CommitIteration the iteration counter was incremented
  bool store = currentIteration-1 != lastStoredResult_ && currentIteration > 1;
  LOG_DBG(opt) << "CheckFinalStoreResults: currentIteration=" << currentIteration << " lastStoredResult="
               << lastStoredResult_ << " store=" << store;
  if(store)
    StoreResults(currentIteration-1);
}


PtrParamNode Optimization::CommitIteration(bool keep_iteration_number)
{
  // store the real cost -> not a scaled one
  objectives.PushBackHistory();

  // store the current design and calculate the design change!
  objectives.PushBackDesign(design);

  // eventually set special result
  EvaluateSpecialResults();

  // also log to info node, append the iteration
  PtrParamNode iteration = optInfoNode->Get(ParamNode::PROCESS)->Get("iteration", ParamNode::APPEND);

  // write the header only once - we might keep the iteration number
  if(log.file) if(objectives.GetHistorySize() == 1) *log.file << log.fileHeader << endl;
  LogFileLine(log.file, iteration); // also ParamNode is to be written
  baseOptimizer_->LogFileLine(log.file, iteration);
  if(log.file) *log.file << endl;

  // this writes the most current solved forward problem via the driver to gid or whatever
  bool store = currentIteration == 0 || commitStride == 1 || (commitStride > 0 && currentIteration % commitStride == 0);
  LOG_TRACE2(opt) << "CommitIteration " << currentIteration << " objective=" << objectives.GetHistoryValue() << " store=" << store;
  if(store)
  {
    StoreResults();
    lastStoredResult_ = currentIteration;
    // see FinalizeStoreResults() !
  }

  // IPOPT does own logging -> otherwise show the user we are alive
  string f = GetIterationFrequency();
  if(optimizer_ != IPOPT_SOLVER && optimizer_ != SNOPT_SOLVER && optimizer_ != KNITRO_SOLVER)
  {
    cout << "iteration " << (currentIteration);
    if(f != "") cout << " f = " << f << "Hz";
    cout << " -> cost = " << objectives.GetHistoryValue() << endl;
  }

  if(!keep_iteration_number)
  {
    currentIteration++;
    problemWithinIteration = 0;
  }

  return iteration;
}

void Optimization::LogFileLine(ofstream* out, PtrParamNode iteration)
{
  // the current design as handy vector
  double change = objectives.design_change.Last();

  if(out)
  {
    *out << currentIteration;
    if(harmonic) *out << "\t" << GetIterationFrequency();

    for(unsigned int i = 0; i < objectives.data.GetSize(); i++)
      *out << "\t" << objectives.data[i]->GetValue();

    *out << "\t" << change << "\t" << problemSolvedCounter;
  }

  iteration->Get("number")->SetValue(currentIteration);
  if(harmonic) iteration->Get("frequency")->SetValue(GetIterationFrequency());

  for(unsigned int i = 0; i < objectives.data.GetSize(); i++)
  {
    Function* f = objectives.data[i];
    iteration->Get(f->type.ToString(f->GetType()))->SetValue(f->GetValue());
    if(f->GetLocal() != NULL)
      iteration->Get("infeasible_" + f->type.ToString(f->GetType()))->SetValue(f->GetLocal()->infeasible);
  }

  iteration->Get("change")->SetValue(change);
  // iteration->Get("problemsSolved")->SetValue(problemSolvedCounter);

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

  for(unsigned int i = 0; i < constraints.all.GetSize(); i++)
  {
    Condition* g = constraints.all[i]; // Now traverse in global mode

    if(g->IsLocalCondition())
    {
      LocalCondition* local = dynamic_cast<LocalCondition*>(g);
      double max  = local->CalcMaxValue();
      double mean = local->CalcMeanValue();
      if(out)
        *out << "\t" << max << "\t" << mean;
      iteration->Get("max_" + g->ToString())->SetValue(max);
      iteration->Get("mean_" + g->ToString())->SetValue(mean);
    }
    else
    {
      double value = g->GetValue();
      if(g->delta_logging) value = value - g->GetBoundValue();
      if(out) *out << "\t" << value;
      // excitation sensitive constraints are printed in the excitation list if there is one
      if(!g->IsExcitationSensitive() || me->excitations.GetSize() < 2)
      {
        iteration->Get(g->ToString(me))->SetValue(value);
        // don't report for local, they should be almost always feasible for MMA, ...
        if(g->GetLocal() != NULL )
          iteration->Get("infeasible_" + g->ToString())->SetValue(g->GetLocal()->infeasible);
      }
    }
  }

  if(out && log.design){
    StdVector<double> d;
    d.Resize(design->GetNumberOfVariables());
    design->WriteDesignToExtern(d.GetPointer(), false);
    for(unsigned int i = 0; i < design->GetNumberOfVariables(); i++){
      *out << "\t" << d[i];
    }
  }
  
  if(out && log.designGradient){
    StdVector<double> d;
    d.Resize(design->GetNumberOfVariables());
    d.window.Set(d);
    design->WriteGradientToExtern(d, DesignElement::COST_GRADIENT, DesignElement::PLAIN, NULL, false);
    for(unsigned int i = 0; i < design->GetNumberOfVariables(); i++){
      *out << "\t" << d[i];
    }
  }

  if(out) out->flush();
}

LinearFormContext* Optimization::GetLinearFormContext(const RegionIdType regionId, StdPDE* pde, const string& integrator, Global::ComplexPart entryType)
{
  Assemble* ass = domain->GetBasePDE()->getPDE_assemble();
  return(ass->GetLinearForm(regionId, pde, integrator, false, entryType));
}


BiLinFormContext* Optimization::GetFormContext(const RegionIdType regionId, StdPDE* pde1, StdPDE* pde2, const string& integrator, bool throw_exception, Global::ComplexPart entryType)
{
  Assemble* ass = domain->GetBasePDE()->getPDE_assemble();
  return(ass->GetBiLinForm(regionId, pde1, pde2, integrator, !throw_exception, entryType));
}

BaseForm* Optimization::GetForm(const RegionIdType regionId, StdPDE* pde1, StdPDE* pde2, const string& integrator, bool throw_exception, Global::ComplexPart entryType)
{
  if (pde2 != NULL){
    BiLinFormContext* bc = GetFormContext(regionId, pde1, pde2, integrator, throw_exception, entryType);
    return bc != NULL ? bc->GetIntegrator() : NULL;
  }
  else {
    LinearFormContext* bc = GetLinearFormContext(regionId, pde1, integrator, entryType);
    return bc != NULL ? bc->GetIntegrator() : NULL;
  }
}

BaseForm* Optimization::GetForm(const RegionIdType reg, Application app1, Application app2, bool throw_exception)
{
  Application a1, a2;
  string integrator = "";

  if(app1 == MECH && (app2 == MECH || app2 == NO_APP))
  {
    a1 = a2 = MECH;
    integrator = "linElastInt";
  }
  if(app1 == ELEC && (app2 == ELEC || app2 == NO_APP))
  {
    a1 = a2 = ELEC;
    integrator = "linGradBDBInt";
  }
  if((app1 == MECH && app2 == ELEC) || (app1 == ELEC && app2 == MECH) || (app1 == PIEZO_COUPLING && app2 == NO_APP))
  {
    a1 = MECH;
    a2 = ELEC;
    integrator = "linPiezoCoupling";
  }
  if(app1 == MASS && (app2 == MASS || app2 == NO_APP))
  {
    a1 = a2 = MASS;
    integrator = "MassInt";
  }

  assert(integrator != "");

  SinglePDE* pde1 = ToPDE(a1, throw_exception);
  SinglePDE* pde2 = ToPDE(a2, throw_exception);

  if(pde1 == NULL || pde2 == NULL)
  {
    if(!throw_exception)
      return NULL;
    else
      EXCEPTION("No PDE for application " << a1 << " resp. " << a2);
  }

  return GetForm(reg, pde1, pde2, integrator, throw_exception);
}


void Optimization::SetPDEs(OptimizationMaterial::System sys)
{
  switch(sys)
  {
  case OptimizationMaterial::MECH:
  case OptimizationMaterial::PIEZOCOUPLING:
    pde = domain->GetSinglePDE("mechanic");
    pdes[MECH] = pde;
    break;

  case OptimizationMaterial::HEAT:
    pde = domain->GetSinglePDE("heatConduction");
    pdes[HEAT] = pde;
    break;

  case OptimizationMaterial::ACOUSTIC:
    pde = domain->GetSinglePDE("acoustic", true);
    pdes[ACOUSTIC] = pde;
    break;

  case OptimizationMaterial::ELEC:
    pde = domain->GetSinglePDE("electrostatic", true);
    pdes[ELEC] = pde;
    break;

  default:
    assert(false);
  }

  // make it more smart when using energy flux for other pdes
  if(objectives.Has(Function::ENERGY_FLUX))
    pdes[ACOUSTIC] = domain->GetSinglePDE("acoustic", true);

  // ELEC is set in PiezoSIMP()
}


SinglePDE* Optimization::ToPDE(Application app, bool throw_exception) const
{
  map<Application, SinglePDE*>::const_iterator it = pdes.find(app);
  if(it != pdes.end())
    return it->second;

  // nothing found
  if(throw_exception)
    EXCEPTION("No PDE '" << app << "' stored");

  return NULL;
}


Optimization::Application Optimization::ToApp(DesignElement::Type dt)
{
  switch(dt)
  {
  case DesignElement::DENSITY:
    return MECH;
  case DesignElement::ACOU_DENSITY:
    return ACOUSTIC;
  case DesignElement::POLARIZATION:
    return ELEC;
  default:
    EXCEPTION("DesignType " << DesignElement::type.ToString(dt) << " doesn't map to Application");
  }
}

DesignElement::Type Optimization::ToDesign(const SinglePDE* pde) const
{
  if(pde->GetName() == "electrostatic") return DesignElement::POLARIZATION;
  if(pde->GetName() == "mechanic") return DesignElement::DENSITY;
  if(pde->GetName() == "acoustic") return DesignElement::ACOU_DENSITY;

  throw Exception("invalid");
}

Optimization::Application Optimization::ToApp(const SinglePDE* pde) const
{
  if(pde->GetName() == "electrostatic") return ELEC;
  if(pde->GetName() == "mechanic") return MECH;
  if(pde->GetName() == "heatConduction") return HEAT;
  if(pde->GetName() == "acoustic") return ACOUSTIC;

  throw Exception("invalid");
}



Optimization::Log::Log()
{
  this->design = false;
  this->designGradient = false;
  this->file = NULL;
  this->fileHeader = "";
}

void Optimization::Log::Init(const string& log_name, PtrParamNode pn_log)
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
