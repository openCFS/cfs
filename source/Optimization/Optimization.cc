#include <def_use_ipopt.hh>
#include <def_use_scpip.hh>
#include <def_use_snopt.hh>
#include <map>

#include "Optimization/Optimization.hh"
#include "Optimization/SIMP.hh"
#include "Optimization/PiezoSIMP.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/OptimalityCondition.hh"
#include "Optimization/LevelSet.hh"
#include "Optimization/EvaluateOnly.hh"
#include "Optimization/ShapeOptimizer.hh"
#include "Optimization/ShapeGrad.hh"
#include "Optimization/GradientCheck.hh"
#include "Optimization/ParamMat.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Driver/assemble.hh"
#include "Driver/assemble.hh"
#include "Driver/basedriver.hh"
#include "Driver/harmonicDriver.hh"
#include "Driver/singleDriver.hh"
#include "PDE/StdPDE.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/mechPDE.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/resultHandler.hh"
#include "DataInOut/programOptions.hh"
#include "General/exception.hh"
#include <boost/filesystem.hpp>
#include "Optimization/ShapeOpt.hh"

// IPOPT, SCPIP and SnOpt are not necessarily linked
#ifdef USE_IPOPT
  #include "Optimization/IPOPTHolder.hh"
#endif
#ifdef USE_SCPIP
  #include "Optimization/SCPIP.hh"
#endif
#ifdef USE_SNOPT
  #include "Optimization/SnOpt.hh"
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

Enum<Optimization::MultipleExcitation::Type>  Optimization::MultipleExcitation::type;

Optimization::Optimization()
{
  this->lastStoredResult_ = -1;
  this->design = NULL;
  this->baseOptimizer_ = NULL;
  this->harmonic = BasePDE::IsComplex(domain->GetDriver()->GetAnalysisType());
  this->currentIteration = 0; // a 1 or 0 can make a lot of difference! 0 is initial design!
  this->problemSolvedCounter = 0;
  this->problemWithinIteration = 0;

  // inject the driver and tell him that we do optimization
  BaseDriver* driver = domain->GetDriver();
  if(driver->GetDriverClass() != BaseDriver::SINGLE_DRIVER)
    throw Exception("optimization not implemented for driver " + driver->GetDriverClass());

  optInfoNode = info->Get("optimization");   // store our info results here
  ParamNode* pn = param->Get("optimization"); // read our parameters from the xml file

  // the tool to solve the optimization problem
  optimizer_ = optimizer.Parse(pn->Get("optimizer")->Get("type"));
  maxIterations = pn->Get("optimizer")->Get("maxIterations")->AsInt();

  // might read a multiObjective problem
  objectives.Read(pn->Get("costFunction")); //
  objectives.ToInfo(optInfoNode->Get(InfoNode::HEADER)->Get("objective"));

  // constraints to be added later -- it is so much easier with the InfoNodes
  this->log.fileHeader = harmonic ? "#iter\tfreq" : "#iter";
  for(unsigned int i = 0; i < objectives.data.GetSize(); i++)
    this->log.fileHeader += "\t" + Objective::type.ToString(objectives.data[i]->GetType());
  this->log.fileHeader += "\tchange\tproblems";


  // multiple excitations are are toggled via attribute. Only if enabled we read the optional element
  // actually part of costFunction - but we store in Optimization itself!
  bool me = pn->Get("costFunction")->Get("multiple_excitation")->AsBool();
  this->multiple_excitation = new MultipleExcitation(me, me ? pn->Get("costFunction")->Get("multipleExcitation", false) : NULL);
  if(me) this->multiple_excitation->ToInfo(optInfoNode->Get(InfoNode::HEADER)->Get("multipleExcitations"));

  // the commit stuff
  string cm = pn->Has("commit") ? pn->Get("commit")->Get("mode")->AsString() : "forward";
  this->commitMode_ = commitMode.Parse(cm);
  this->commitStride = pn->Has("commit") ? pn->Get("commit")->Get("stride")->AsInt() : 1;
  optInfoNode->Get("commit")->Get("mode")->SetValue(cm);
  optInfoNode->Get("commit")->Get("stride")->SetValue(commitStride);
  
  // write out the directory where the HALTOPT file will be searched for
  optInfoNode->Get("haltopt_directory")->SetValue(fs::current_path().directory_string());

  // the constraints are optional and might not be real constraints!
  StdVector<ParamNode*> list = pn->GetList("constraint");
  InfoNode* in = optInfoNode->Get(InfoNode::HEADER)->Get("constraints");
  for(unsigned int i = 0; i < list.GetSize(); i++)
  {
    // the constraint is either added to constraints or outputs if mode is observation
    // homogenization constraints are "super constraints" which might "blow" up to 4 or 9 constraints
     Condition::AddCondition(list[i], constraints, outputs);
  }

  // first the constraints then the outputs!
  for(unsigned int i = 0; i < constraints.GetSize(); i++)
  {
    // if in the 'logging' element deltaConstraints is enabled, we modify the output!
    string delta = constraints[i].delta_logging ? "delta_" : "";
    this->log.fileHeader += "\t" + delta + constraints[i].ToString();
    constraints[i].ToInfo(in->Get("constraint", InfoNode::APPEND));
  }

  for(unsigned int i = 0; i < outputs.GetSize(); i++)
  {
    string delta = outputs[i].delta_logging ? "delta_" : "";
    this->log.fileHeader += "\t" + delta + outputs[i].ToString();
    outputs[i].ToInfo(in->Get("observation", InfoNode::APPEND));
  }

  log.Init(pn->Get("log")->AsString(), pn->Get("logging", false)); // is fail save

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
  delete multiple_excitation; multiple_excitation = NULL;
}

void Optimization::PostInit()
{
  if(param->Has("loadErsatzMaterial")) { // if loadErsatzMaterial is used with optimization specifying a starting point, we have to load it here, before scaling ist done
    domain->ReadErsatzMaterial(param->Get("loadErsatzMaterial"));
  }
}

void Optimization::PostInitSecond()
{
  ParamNode* opt = param->Get("optimization")->Get("optimizer");

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
           baseOptimizer_ = new SCPIP(this, opt);
         #else
           throw Exception("CFS++ was compiled w/o SCPIP");
         #endif
         break;
         
    case SNOPT_SOLVER:
         #ifdef USE_SNOPT
           baseOptimizer_ = new SnOpt(this, opt);
         #else
           throw Exception("CFS++ was compiled w/o SnOpt");
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

    default: throw Exception("optimizer not implemented");
  }
  if (this->log.design) {
    for (unsigned int i = 0; i < design->GetNumberOfVariables(); i++) {
      this->log.fileHeader += "\t";
      this->log.fileHeader += "design";
    }
  }
  if (this->log.designGradient) {
    for (unsigned int i = 0; i < design->GetNumberOfVariables(); i++) {
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
  Objective::type.SetName("Objective::Type");
  Objective::type.Add(Objective::MULTI_OBJECTIVE, "multiObjective");
  Objective::type.Add(Objective::COMPLIANCE, "compliance");
  Objective::type.Add(Objective::OUTPUT, "output");
  Objective::type.Add(Objective::DYNAMIC_OUTPUT, "dynamicOutput");
  Objective::type.Add(Objective::ABS_DYN_OUTPUT_SQUARED, "absDynamicOutputSquared");
  Objective::type.Add(Objective::GLOBAL_DYNAMIC_COMPLIANCE, "globalDynamicCompliance");
  Objective::type.Add(Objective::CONJUGATE_COMPLIANCE, "conjugateCompliance");
  Objective::type.Add(Objective::VOLUME, "volume");
  Objective::type.Add(Objective::TRACKING, "tracking");
  Objective::type.Add(Objective::ELEC_ENERGY, "elecEnergy");
  Objective::type.Add(Objective::HOMOGENIZATION_TENSOR, "homTensor");
  Objective::type.Add(Objective::HOMOGENIZATION_E11, "homE11");
  Objective::type.Add(Objective::HOMOGENIZATION_TRACKING, "homTracking");
  Objective::type.Add(Objective::POISSONS_RATIO, "poissonsRatio");
  Objective::type.Add(Objective::YOUNGS_MODULUS, "homYoungsModulus");
  Objective::type.Add(Objective::TYCHONOFF, "tychonoff");
  Objective::type.Add(Objective::TEMPERATURE, "temperature");

  Condition::name.SetName("Constraint::Name");
  Condition::name.Add(Condition::VOLUME, "volume");
  Condition::name.Add(Condition::COMPLIANCE, "compliance");
  Condition::name.Add(Condition::GREYNESS, "greyness");
  Condition::name.Add(Condition::GAUSS_GREYNESS, "gaussGreyness");
  Condition::name.Add(Condition::TRACKING, "tracking");
  Condition::name.Add(Condition::REALVOLUME, "realvolume");
  Condition::name.Add(Condition::HOMOGENIZATION_TENSOR, "homTensor");
  Condition::name.Add(Condition::HOMOGENIZATION_TRACKING, "homTracking");
  Condition::name.Add(Condition::ISOTROPY, "isotropy");

  Condition::type.SetName("Constraint::Type");
  Condition::type.Add(Condition::EQUAL, "equal");
  Condition::type.Add(Condition::LOWER_BOUND, "lowerBound");
  Condition::type.Add(Condition::UPPER_BOUND, "upperBound");

  optimizer.SetName("Optimization::Optimizer");
  optimizer.Add(OPTIMALITY_CONDITION, "optimalityCondition");
  optimizer.Add(IPOPT_SOLVER, "ipopt");
  optimizer.Add(SCPIP_SOLVER, "scpip");
  optimizer.Add(SNOPT_SOLVER, "snopt");
  optimizer.Add(SHAPE_SOLVER, "shapeOpt");
  optimizer.Add(EVALUATE_INITIAL_DESIGN, "evaluateInitialDesign");
  optimizer.Add(GRADIENT_CHECK, "gradientCheck");

  ErsatzMaterial::method.SetName("ErsatzMaterial::Method");
  ErsatzMaterial::method.Add(ErsatzMaterial::SIMP_METHOD, "simp");
  ErsatzMaterial::method.Add(ErsatzMaterial::PARAM_MAT, "paramMat");
  ErsatzMaterial::method.Add(ErsatzMaterial::SHAPE_GRAD, "shapeGrad");
  ErsatzMaterial::method.Add(ErsatzMaterial::SHAPE_OPT, "shapeOpt");
  ErsatzMaterial::method.Add(ErsatzMaterial::SHAPE_PARAM_MAT, "shapeParamMat");
  
  ErsatzMaterial::commitMode.SetName("ErsatzMaterial::CommitMode");
  ErsatzMaterial::commitMode.Add(ErsatzMaterial::FORWARD, "forward");
  ErsatzMaterial::commitMode.Add(ErsatzMaterial::ADJOINT, "adjoint");
  ErsatzMaterial::commitMode.Add(ErsatzMaterial::BOTH, "both_cases");

  OptimizationMaterial::system.SetName("OptimizationMaterial::System");
  OptimizationMaterial::system.Add(OptimizationMaterial::PIEZO, "piezo");
  OptimizationMaterial::system.Add(OptimizationMaterial::MECHANIC, "mechanic");
  OptimizationMaterial::system.Add(OptimizationMaterial::HEAT, "heat");

  application.SetName("Optimization::Application");
  application.Add(NO_APP, "no_app");
  application.Add(MECH, "mech");
  application.Add(MASS, "mass");
  application.Add(ELEC, "elec");
  application.Add(PIEZO_COUPLING, "piezoCoupling");
  application.Add(PRESSURE, "pressure");
  application.Add(CHARGE_DENSITY, "chargeDensity");

  LevelSet::Action::type.SetName("LevelSet::Action::Type");
  LevelSet::Action::type.Add(LevelSet::Action::SIGNED_DISTANCE_FIELD, "signedDistanceField");
  LevelSet::Action::type.Add(LevelSet::Action::TRIVIAL_HOLE, "trivialHole");
  LevelSet::Action::type.Add(LevelSet::Action::DO_SHAPE_STEP, "shapeStep");

  MultipleExcitation::type.SetName("Optimization::MultipleExcitation::Type");
  MultipleExcitation::type.Add(MultipleExcitation::NO_TYPE, "no_type:q"
      "");
  MultipleExcitation::type.Add(MultipleExcitation::FIXED_WEIGHT, "fixed_weights");
  MultipleExcitation::type.Add(MultipleExcitation::META_OBJECTIVE, "meta_objective");
  MultipleExcitation::type.Add(MultipleExcitation::HOMOGENIZATION_TEST_STRAINS, "homogenizationTestStrains");
}


bool Optimization::DoStopOptimization()
{
  InfoNode* in = optInfoNode->Get(InfoNode::SUMMARY)->Get("break");
  // check if the HALTOPT file exists
  if(fs::exists("HALTOPT"))
  {
    bool good = fs::remove("HALTOPT");
    if(!good) throw new Exception("Could not remove file 'HALTOPT' after detection");
    in->Get("converged")->SetValue("no");
    in->Get("reason")->SetValue("Detected file 'HALTOPT'");
    return true;
  }
  
  // this currently only implements relative stopping rule
  Objective* cost = objectives.data[0];

  // we need a minimum number of iterations to be sure we are in a minimum
  unsigned int hs = objectives.GetHistorySize();
  if(hs <= cost->stop.queue) return false;

  for(unsigned int i = hs-1; i >= (hs - cost->stop.queue); i--)
  {
    double delta = objectives.GetHistoryValue(true, i) - objectives.GetHistoryValue(true, i-1);
    double rel = abs(delta / objectives.GetHistoryValue(true, i));
    if(rel > cost->stop.value) return false;
  }

  // the relative values for the whole queue are smaller than the requirement -> we are done! :)
  in->Get("converged")->SetValue("practically");
  in->Get("reason")->SetValue("Too small change in objective function");
  in->Get("reason")->Get("queue")->SetValue(cost->stop.queue);
  in->Get("reason")->Get("relative")->SetValue(cost->stop.value);
  return true;
}


/** read only the very basic stuff */
Optimization* Optimization::CreateInstance()
{
  // set the enums we need
  Optimization::SetEnums(); // sets also ErsatzMaterial::Method
  DesignElement::SetEnums();

  if(!param->Has("optimization")) return NULL;

  // we assume ersatz material, currently there is nothing else.
  // note, we read method again in the ersatz material constructor.
  ParamNode* em = param->Get("optimization")->Get("ersatzMaterial");
  ErsatzMaterial::Method method = ErsatzMaterial::method.Parse(em->Get("method"));
  
  Optimization* opt = NULL;
  
  switch(method)
  {
  case ErsatzMaterial::SIMP_METHOD:
    switch(OptimizationMaterial::system.Parse(em->Get("material")))
    {
    case OptimizationMaterial::MECHANIC:
      opt = new SIMP();
      break;
      
    case OptimizationMaterial::PIEZO:
      opt = new PiezoSIMP();
      break;
      
    default: assert(false);
    }
    break;
    
  // FreeMat, ShapeGrad, ...
  case ErsatzMaterial::PARAM_MAT:
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

  ResultHandler* rh = domain->GetResultHandler();
  unsigned int mss = domain->GetDriver()->GetActSequenceStep();
  rh->BeginMultiSequenceStep(mss, BasePDE::TRANSIENT, 999); // max steps is high

  Exception* e = NULL;
  try
  {
    baseOptimizer_->SolveOptimizationProblem();
  }
  catch(Exception& ex)
  {
    e = new Exception(ex);
  }

  // do the finally - try to write results even if the optimizer broke down
  FinalizeStoreResults(); // when we have strides the results are written
  rh->FinishMultiSequenceStep();
  rh->Finalize();
  if(e != NULL) throw *e;
}


InfoNode* Optimization::CreateAdjointAnalysisIdNode()
{
  BaseDriver* driver = domain->GetDriver();
  InfoNode* base = driver->GetAnalysisId();
  InfoNode* in = driver->CreateAnalysisIdChild(base, "adjoint");
  
  return in;
}


void Optimization::SolveStateProblem(Excitation* excite)
{
  BaseDriver* driver = domain->GetDriver();

  // this is the forward problem. Store the analysis_id in the driver such that
  // we can aquire it for the adjoint problem
  InfoNode* analysis_id = excite == NULL ? driver->CreateAnalysisId("iter", currentIteration)
                                         : driver->CreateAnalysisId("iter", currentIteration, "excite", excite->index);
  
  bool reAssembleMatrices = true;
  if(excite != NULL && excite->index > 0)
  {
    reAssembleMatrices = false;
  }
                                         
  // Do not store the results. This is to be done in CommitIteration
  if(!harmonic || excite == NULL) 
    driver->SolveProblem(false, analysis_id, reAssembleMatrices);
  else
    dynamic_cast<HarmonicDriver*>(driver)->ComputeFrequencyStep(excite->f_link->step, analysis_id);

  problemSolvedCounter++;
  problemWithinIteration++;
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
  domain->GetDriver()->StoreResults(step_val == -1 ? currentIteration : step_val);
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


InfoNode* Optimization::CommitIteration(bool keep_iteration_number)
{
  // store the real cost -> not a scaled one
  objectives.PushBackHistory();

  // eventually set special result
  EvaluateSpecialResults();

  // this writes the most current solved forward problem via the driver to gid or whatever
  bool store = currentIteration == 0 || commitStride == 1 || (commitStride > 0 && currentIteration % commitStride == 0);
  LOG_TRACE2(opt) << "CommitIteration " << currentIteration << " ojective=" << objectives.GetHistoryValue() << " store=" << store;
  if(store)
  {
    StoreResults();
    lastStoredResult_ = currentIteration;
    // see FinalizeStoreResults() !
  }

  // save this iteration
  design->WriteDesignToExtern(last_iteration.GetPointer());

  // also log to info node, append the iteration
  InfoNode* iteration = optInfoNode->Get(InfoNode::PROCESS)->Get("iteration", InfoNode::APPEND);

  // write the header only once - we might keep the iteration number
  if(log.file) if(objectives.GetHistorySize() == 1) *log.file << log.fileHeader << endl;
  LogFileLine(log.file, iteration); // also InfoNode is to be written
  baseOptimizer_->LogFileLine(log.file, iteration);
  if(log.file) *log.file << endl;

  // IPOPT does own logging -> otherwise show the user we are alive
  std::string f = GetIterationFrequency();
  if(optimizer_ != IPOPT_SOLVER && optimizer_ != SNOPT_SOLVER)
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

void Optimization::LogFileLine(ofstream* out, InfoNode* iteration)
{
  // calculate the relative cost change
  int hs = objectives.GetHistorySize();
  double change = 0.0;
  if(hs >= 2)
  {
    double last = objectives.GetHistoryValue(true, hs-1);
    change = (last - objectives.GetHistoryValue(true, hs-2)) / last;
  }

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
    iteration->Get(Objective::type.ToString(objectives.data[i]->GetType()))->SetValue(objectives.data[i]->GetValue());

  iteration->Get("change")->SetValue(change);
  iteration->Get("problemsSolved")->SetValue(problemSolvedCounter);

  for(unsigned int i = 0; i < constraints.GetSize(); i++)
  {
    double value = CalcConstraint(&constraints[i]);
    if(constraints[i].delta_logging) value = value - constraints[i].value;
    if(out) *out << "\t" << value;
    iteration->Get(constraints[i].ToString())->SetValue(value);
  }

  for(unsigned int i = 0; i < outputs.GetSize(); i++)
  {
    double value = CalcConstraint(&outputs[i]);
    if(outputs[i].delta_logging) value = value - outputs[i].value;
    if(out) *out << "\t" << value;
    iteration->Get(outputs[i].ToString())->SetValue(value);
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
    design->WriteGradientToExtern(d.GetPointer(), DesignElement::COST_GRADIENT, DesignElement::PLAIN, NULL, false);
    for(unsigned int i = 0; i < design->GetNumberOfVariables(); i++){
      *out << "\t" << d[i];
    }
  }

  if(out) out->flush();
}


Optimization::MultipleExcitation::MultipleExcitation(bool multiple, ParamNode* pn)
{
  // set defaults
  this->stride = 1;
  this->damping = 1.0;
  this->max_gain = 1e4;
  this->multiple_excitation_ = multiple;
  this->type_ = NO_TYPE; // to be eventually overwritten soon

  // if disabled, we don't read anything
  if(!multiple || pn == NULL) return;

  this->type_ = type.Parse(pn->Get("type"));

  // adjust defaults
  if(pn->Has("adjustWeights")) {
    this->stride   = pn->Get("adjustWeights")->Get("stride")->AsInt();
    this->max_gain = pn->Get("adjustWeights")->Get("max_gain")->AsDouble();
    this->damping  = pn->Get("adjustWeights")->Get("damping")->AsDouble();
  }
}

void Optimization::MultipleExcitation::ToInfo(InfoNode* in) const
{
  if(!multiple_excitation_) return;

  if(type_ == META_OBJECTIVE)  {
    InfoNode* in_ = in->Get("metaObjective");
    in_->Get("type")->SetValue("best_value");
    in_->Get("damping")->SetValue(damping);
    in_->Get("stride")->SetValue(stride);
    in_->Get("max_gain")->SetValue(max_gain);
  } else
    in->Get("type")->SetValue(type.ToString(type_));
}

Condition& Optimization::GetConstraint(Condition::Name name, DesignElement::Type design)
{
  // be save and check for uniqueness!
  int count = 0;

  for(unsigned int i = 0; i < constraints.GetSize(); i++)
    if(constraints[i].GetName() == name &&
       (design != DesignElement::NO_TYPE ? constraints[i].design == design : true))
          count++;

  for(unsigned int i = 0; i < outputs.GetSize(); i++)
    if(outputs[i].GetName() == name &&
       (design != DesignElement::NO_TYPE ? outputs[i].design == design : true))
          count++;

  if(count > 1)
    throw Exception("constraint " + Condition::name.ToString(name) + "is not unique");

  for(unsigned int i = 0; i < constraints.GetSize(); i++)
    if(constraints[i].GetName() == name) return constraints[i];

  for(unsigned int i = 0; i < outputs.GetSize(); i++)
    if(outputs[i].GetName() == name) return outputs[i];

  throw Exception("no constraint " + Condition::name.ToString(name) + " found");
}


Excitation::Excitation()
{
  this->index = -1; // must be updated
  this->frequency = -1.0;
  this->f_link = NULL;
  this->weight = 1.0;
  this->normalized_weight = 1.0;
  this->apply_linForms = false;
}


void Excitation::Apply()
{
  Assemble* a = domain->GetBasePDE()->getPDE_assemble();
  if(apply_linForms){
    a->SetLoads(loads);
    a->SetLinForms(linForms);
  }else if(! loads.IsEmpty()){
    a->SetLoads(loads);
  }
  // a frequency cannot really be applied but has to be used as parameter
  // in the driver call
}

double Excitation::GetOmega()
{
  if(frequency < 0.0)
    EXCEPTION("No frequency given");
  return 2 * M_PI * frequency;
}

double Excitation::GetFactor(Objective* cost)
{
  double factor = 1.0; // default
  
  if(cost->FactorOmegaOmega())
  {
    if(frequency < 0.0) EXCEPTION("No frequency given");
    factor *= 4 * M_PI * M_PI * frequency * frequency;
  }
  
  return factor;
}

double Excitation::GetWeightedFactor(Objective* cost)
{
  return normalized_weight * GetFactor(cost);
}

void Excitation::ReadTrackings(ParamNode* ts){
  StdPDE* mech = domain->GetStdPDE("mechanic");
  trackings.Clear();
  StdVector<ParamNode*> tracking_list = ts->GetChildren();
  for(unsigned int i = 0; i < tracking_list.GetSize(); i++){
    std::string name, dof, value;
    dof = "";
    tracking_list[i]->Get("name", name);
    tracking_list[i]->Get("dof", dof, false);
    tracking_list[i]->Get("value", value);
    shared_ptr<LoadBc> actLoad( new LoadBc );
    shared_ptr<EntityList> actList = domain->GetGrid()->GetEntityList( EntityList::NODE_LIST, name, EntityList::NAMED_NODES );
    actLoad->entities = actList;
    actLoad->eqnMap = mech->GetEqnMap();
    if ( dof == "" ) {
      actLoad->dof = 1;
    } else {
      actLoad->dof = mech->GetResultInfo(MECH_DISPLACEMENT)->GetDofIndex(dof);
    }
    actLoad->value = value;
    trackings.Push_back(actLoad);
  }
}

void Excitation::ReadLoads(ParamNode* ls)
{
  apply_linForms = true;
  
  // loads
  loads.Clear();
  MechPDE* mech = (MechPDE*)domain->GetSinglePDE("mechanic");
  mech->ReadLoads(ls->GetList("load"), loads);

  linForms.clear();
  
  // pressures
  StdVector<shared_ptr<EntityList> > pressSurf;
  StdVector<std::string> pressVals;
  StdVector<std::string> pressPhase;
  mech->ReadPressureLoadsFromXML(ls, pressSurf, pressVals, pressPhase);
  mech->DefinePressureIntegrators(pressSurf, pressVals, pressPhase, &linForms);
  
  // regionLoads
  std::map<RegionIdType, SinglePDE::RegionLoad> regionLoads;
  mech->ReadRegionLoadsFromXML(ls, regionLoads);
  mech->DefineRegionLoadIntegrators(regionLoads, &linForms);
  
  // all already set linear Forms
  std::set<LinearFormContext*>* assLinForms = domain->GetBasePDE()->getPDE_assemble()->GetLinForms();
  linForms.insert(assLinForms->begin(), assLinForms->end());
}

void Excitation::ReadTestStrain(const Vector<double>& vec)
{
  apply_linForms = true;

  this->test_strain = vec;

  loads.Clear();
  linForms.clear();
  MechPDE* mech = dynamic_cast<MechPDE*>(domain->GetSinglePDE("mechanic"));
  mech->DefineTestStrainIntegrators(vec, &linForms);

  // all already set linear Forms
  std::set<LinearFormContext*>* assLinForms = domain->GetBasePDE()->getPDE_assemble()->GetLinForms();
  linForms.insert(assLinForms->begin(), assLinForms->end());
}

Optimization::Log::Log()
{
  this->design = false;
  this->designGradient = false;
  this->file = NULL;
  this->fileHeader = "";
}

void Optimization::Log::Init(const string& log_name, ParamNode* pn_log)
{
  if(log_name != "false")
  {
    string name = log_name == "[problem]" || log_name == "true" ? (progOpts->GetSimName() + ".plot.dat") : log_name;
    file = new ofstream(name.c_str());
    if(file == NULL)
      throw Exception("cannot open log file " + name + " for writing");

    if(pn_log != NULL)
    {
      design = false;
      pn_log->Get("design", design, false);
      designGradient = false;
      pn_log->Get("designGradient", designGradient, false);
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
