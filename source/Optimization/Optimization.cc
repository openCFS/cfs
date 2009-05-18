#include <def_use_ipopt.hh>
#include <def_use_scpip.hh>

#include "Optimization/Optimization.hh"
#include "Optimization/SIMP.hh"
#include "Optimization/PiezoSIMP.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/OptimalityCondition.hh"
#include "Optimization/LevelSet.hh"
#include "Optimization/EvaluateOnly.hh"
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

// IPOPT and SCPIP are not necessarily linked
#ifdef USE_IPOPT
  #include "Optimization/IPOPTHolder.hh"
#endif
#ifdef USE_SCPIP
  #include "Optimization/SCPIP.hh"
#endif

using namespace CoupledField;
using namespace std;
namespace fs = boost::filesystem;



DECLARE_LOG(opt)
DEFINE_LOG(opt, "opt")

// instantiation of the static elements
Enum<Optimization::ObjectiveType>    Optimization::objectiveType;
Enum<Optimization::Optimizer>        Optimization::optimizer;
Enum<Optimization::Application>      Optimization::application;
Enum<Optimization::CommitMode>       Optimization::commitMode;

Optimization::Objective::Objective(ParamNode* pn, bool harmonic)
{
  // multiple excitation is handled in Optimization itself!

  // the current value -> check <Get/Set>Value() when altering the presets!
  this->value_       = -1.0;
  this->pn           = pn;
  this->harmonic_    = harmonic;
  this->omega_omega_ = pn->Has("factor") ? pn->Get("factor")->Get("omega_omega")->AsBool() : false;
  if(!harmonic && omega_omega_)
    throw Exception("It makes no sense to set costFunction/factor/omega_omega in static optimization");
  
  this->tychonoff_  = pn->Has("tychonoff") ? pn->Get("tychonoff")->AsDouble() : 0.0;

  type = objectiveType.Parse(pn->Get("type"));

  if(pn->Get("stoppingRule")->AsString() != "relative")
      throw Exception("stopping rule not implemented yet");

  task = pn->Get("task")->AsString() == "minimize" ? MINIMIZE : MAXIMIZE;

  stop.value = pn->Get("value")->AsDouble();

  stop.queue = pn->Get("queue")->AsUInt();
  if(stop.queue < 1) throw Exception("minimal queue value for stopping rule is 1");
}

void Optimization::Objective::ToInfo(InfoNode* in) const
{
  in->Get("type")->SetValue(objectiveType.ToString(type));
  in->Get("task")->SetValue(task == MINIMIZE ? "minimize" : "maximize");
  if(harmonic_)
    in->Get("factor")->Get("omega_omega")->SetValue(omega_omega_);
}

double Optimization::Objective::GetValue()
{
  return value_;
}

void Optimization::Objective::SetValue(double val)
{
  value_ = val;
}



Optimization::Optimization()
{
  this->lastStoredResult_ = -1;
  this->logFile_ = NULL;
  this->design = NULL;
  this->baseOptimizer_ = NULL;
  this->harmonic = BasePDE::IsComplex(domain->GetDriver()->GetAnalysisType());
  this->currentIteration = 0; // a 1 or 0 can make a lot of difference! 0 is initial design!
  this->problemSolvedCounter = 0;
  this->problemWithinIteration = 0;
  // constraints to be added later -- it is so much easier with the InfoNodes
  this->logFileHeader = harmonic ? "#iter\tfreq" : "#iter";
  this->logFileHeader += "\tcost\tchange\tproblems";

  // inject the driver and tell him that we do optimization
  BaseDriver* driver = domain->GetDriver();
  if(driver->GetDriverClass() != BaseDriver::SINGLE_DRIVER)
    throw Exception("optimization not implemented for driver " + driver->GetDriverClass());

  optInfoNode = info->Get("optimization");   // store our info results here
  ParamNode* pn = param->Get("optimization"); // read our parameters from the xml file

  // the tool to solve the optimization problem
  optimizer_ = optimizer.Parse(pn->Get("optimizer")->Get("type"));
  maxIterations = pn->Get("optimizer")->Get("maxIterations")->AsInt();

  // the cost function is mandatory
  this->cost = new Objective(pn->Get("costFunction"), harmonic);
  this->cost->ToInfo(optInfoNode->Get(InfoNode::HEADER)->Get("objective", InfoNode::APPEND));

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

  // the constraints are optional and might not be real constraints!
  StdVector<ParamNode*> list = pn->GetList("constraint");
  InfoNode* in = optInfoNode->Get(InfoNode::HEADER)->Get("constraints");
  for(unsigned int i = 0; i < list.GetSize(); i++)
  {
     // the index is the current constraints size, it such works also if there are outputs
     Condition g(list[i], constraints.GetSize());
     g.ToInfo(in->Get("constraint", InfoNode::APPEND));
     if(g.active) constraints.Push_back(g);
             else outputs.Push_back(g);
  }
  // first the constraints then the outputs!
  for(unsigned int i = 0; i < constraints.GetSize(); i++)
    this->logFileHeader += "\t" + constraints[i].ToString();
  for(unsigned int i = 0; i < outputs.GetSize(); i++)
    this->logFileHeader += "\t" + outputs[i].ToString();


  // We almost always make a gnuplot log file
  string log_name = pn->Get("log")->AsString();
  if(log_name != "false")
  {
    if(log_name == "[problem]" || log_name == "true") log_name = progOpts->GetSimName() + ".plot.dat";
    logFile_ = new ofstream(log_name.c_str());
    if(logFile_ == NULL)
      throw Exception("cannot open log file " + pn->Get("log")->AsString() + " for writing");
    this->logDesign = pn->Get("logdesign")->AsBool();
    this->logDesignGradient = pn->Get("logdesigngradient")->AsBool();
  }
  
  // remove a stop file, if found
  std::string command = "rm -f ./HALTOPT";
  std::system( command.c_str() );
}

Optimization::~Optimization()
{
  // if write to file close it
  if(logFile_ != NULL)
  {
    logFile_->close();
    delete logFile_;
    logFile_ = NULL;
  }

  delete cost; cost = NULL;
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
  if(param->Has("loadErsatzMaterial")) { // if loadErsatzMaterial is used with optimization specifying a starting point, we have to load it here, before scaling ist done
    domain->ReadErsatzMaterial(param->Get("loadErsatzMaterial"));
  }

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

    case OPTIMALITY_CONDITION:
         baseOptimizer_ = new OptimalityCondition(this, opt);
         break;

    case LEVEL_SET:
         baseOptimizer_ = new LevelSet(this, opt);
         break;

    case EVALUATE_INITIAL_DESIGN:
         baseOptimizer_ = new EvaluateOnly(this, opt);
         break;
         
    case GRADIENT_CHECK:
         baseOptimizer_ = new GradientCheck(this, opt);
         break;

    default: throw Exception("optimizer not implemented");
  }
  if (this->logDesign) {
    for (unsigned int i = 0; i < design->GetNumberOfVariables(); i++) {
      this->logFileHeader += "\t";
      this->logFileHeader += "design";
    }
  }
  if (this->logDesignGradient) {
    for (unsigned int i = 0; i < design->GetNumberOfVariables(); i++) {
      this->logFileHeader += "\t";
      this->logFileHeader += "designGradient";
    }
  }
  design->SetOptimizer(baseOptimizer_);
  // add plot logging of the optimizer
  this->logFileHeader += baseOptimizer_->LogFileHeader();
}


void Optimization::SetEnums()
{
  objectiveType.SetName("Optimization::ObjectiveType");
  objectiveType.Add(COMPLIANCE, "compliance");
  objectiveType.Add(OUTPUT, "output");
  objectiveType.Add(DYNAMIC_OUTPUT, "dynamicOutput");
  objectiveType.Add(ABS_DYN_OUTPUT_SQUARED, "absDynamicOutputSquared");
  objectiveType.Add(GLOBAL_DYNAMIC_COMPLIANCE, "globalDynamicCompliance");
  objectiveType.Add(CONJUGATE_COMPLIANCE, "conjugateCompliance");
  objectiveType.Add(VOLUME, "volume");
  objectiveType.Add(TRACKING, "tracking");
  objectiveType.Add(ELEC_ENERGY, "elecEnergy");

  Condition::name.SetName("Contraint::Name");
  Condition::name.Add(Condition::VOLUME, "volume");
  Condition::name.Add(Condition::COMPLIANCE, "compliance");
  Condition::name.Add(Condition::GREYNESS, "greyness");
  Condition::name.Add(Condition::GAUSS_GREYNESS, "gaussGreyness");
  Condition::name.Add(Condition::TRACKING, "tracking");
  Condition::name.Add(Condition::REALVOLUME, "realvolume");

  Condition::type.SetName("Contraint::Type");
  Condition::type.Add(Condition::EQUAL, "equal");
  Condition::type.Add(Condition::LOWER_BOUND, "lowerBound");
  Condition::type.Add(Condition::UPPER_BOUND, "upperBound");

  optimizer.SetName("Optimization::Optimizer");
  optimizer.Add(OPTIMALITY_CONDITION, "optimalityCondition");
  optimizer.Add(IPOPT_SOLVER, "ipopt");
  optimizer.Add(SCPIP_SOLVER, "scpip");
  optimizer.Add(LEVEL_SET, "levelSet");
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
  LevelSet::Action::type.Add(LevelSet::Action::TOP_GRAD, "topGrad");
  LevelSet::Action::type.Add(LevelSet::Action::TRIVIAL_HOLE, "trivialHole");
  LevelSet::Action::type.Add(LevelSet::Action::FIRST_ORDER_FD, "firstOrderFiniteDifferences");
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

  // we need a minimum number of itations to be sure we are in a minimum
  if(cost->history.GetSize() <= cost->stop.queue) return false;

  for(unsigned int i = cost->history.GetSize()-1; i >= (cost->history.GetSize() - cost->stop.queue); i--)
  {
    double delta = abs(cost->history[i] - cost->history[i-1]);
    double rel = abs(delta / cost->history[i]);
    if(rel > cost->stop.value) return false;
  }

  // the relative values for the whole queue are smaller than the requirement -> we are done! :)
  in->Get("converged")->SetValue("practically");
  in->Get("reason")->SetValue("Too small change in objective function");
  in->Get("reason")->Get("queue")->SetValue(cost->stop.queue);
  in->Get("reason")->Get("realtive")->SetValue(cost->stop.value);
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
  // note, we read method again in the ersat material constructor.
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
  default: throw Exception("Optimization not implemented");
  }
  
  // second initialzation phase, constrcuts material
  opt->PostInit();
  // third initizalization phase, constructs optimizer
  opt->PostInitSecond();

  return opt;
}

void Optimization::SolveProblem()
{
  // one driver is one multisequence step. We do this stuff here
  // and call the driver->StoreResults() multiple times f

  ResultHandler* rh = domain->GetResultHandler();
  unsigned int mss = domain->GetDriver()->GetActSequenceStep();
  rh->BeginMultiSequenceStep(mss, BasePDE::TRANSIENT, 9999); // max steps is high
  baseOptimizer_->SolveProblem();
  FinalizeStoreResults(); // when we have strides the results are written
  rh->FinishMultiSequenceStep();
  rh->Finalize();
}

string Optimization::GetSolveComment(Excitation* excite)
{
  ostringstream os;
  os << "iter"<< currentIteration;
  if(excite != NULL)
  {
    if(excite->f_link > 0) os << "_f_" << excite->frequency;
    else os << "_f_load_idx" << excite->index;
  }
  if(problemWithinIteration > 0) os << "_cnt" << problemWithinIteration;
  return os.str();
}

void Optimization::SolveStateProblem(Excitation* excite)
{
  BaseDriver* driver = domain->GetDriver();

  // we solve, give a part of the filename in case we use export linear system
  // but do not store the results. This is to be done in CommitIteration
  string comment = GetSolveComment(excite);
  if(!harmonic || excite == NULL)
    driver->SolveProblem(false, comment);
  else
    dynamic_cast<HarmonicDriver*>(driver)->ComputeFrequencyStep(excite->f_link->step, comment);

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
  // we assume the first element (1) be in the lower left corner and then a lexicaly
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
  cost->history.Push_back(cost->GetValue());

  // eventually set special result
  EvaluateSpecialResults();

  // this writes the most current solved forward problem via the driver to gid or whatever
  bool store = currentIteration == 0 || commitStride == 1 || (commitStride > 0 && currentIteration % commitStride == 0);
  LOG_TRACE2(opt) << "CommitIteration " << currentIteration << " ojective=" << cost->GetValue() << " store=" << store;
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

  if(logFile_)
  {
    // write the header only once - we might keep the iteration number
    if(cost->history.GetSize() == 1) *logFile_ << logFileHeader << endl;
    LogFileLine(logFile_, iteration);
    baseOptimizer_->LogFileLine(logFile_, iteration);
    *logFile_ << endl;
  }

  // IPOPT does own logging -> otherwise show the user we are alive
  std::string f = GetIterationFrequency();
  if(optimizer_ != IPOPT_SOLVER)
  {
    cout << "iteration " << (currentIteration);
    if(f != "") cout << " f = " << f << "Hz";
    cout << " -> cost = " << cost->history.Last() << endl;
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
  double change = cost->history.GetSize() < 2 ? 0.0 : (cost->history.Last()
         - cost->history[cost->history.GetSize() - 2]) / cost->history.Last();


  *out << currentIteration;
  if(harmonic) *out << "\t" << GetIterationFrequency();
  *out << "\t" << cost->history.Last() << "\t" << change
       << "\t" << problemSolvedCounter;

  iteration->Get("number")->SetValue(currentIteration);
  if(harmonic) iteration->Get("frequency")->SetValue(GetIterationFrequency());
  iteration->Get("objective")->SetValue(cost->history.Last());
  iteration->Get("change")->SetValue(change);
  iteration->Get("problemsSolved")->SetValue(problemSolvedCounter);

  for(unsigned int i = 0; i < constraints.GetSize(); i++)
  {
    double value = CalcConstraint(&constraints[i]);
    *out << "\t" << value;
    iteration->Get(constraints[i].ToString())->SetValue(value);
  }

  for(unsigned int i = 0; i < outputs.GetSize(); i++)
  {
    double value = CalcConstraint(&outputs[i]);
    *out << "\t" << value;
    iteration->Get(outputs[i].ToString())->SetValue(value);
  }

  for(unsigned int i = 0; i < outputs.GetSize(); i++)    
    *out << "\t" << CalcConstraint(&outputs[i]);
  
  if(logDesign){
    StdVector<double> d;
    d.Resize(design->GetNumberOfVariables());
    design->WriteDesignToExtern(d.GetPointer(), false);
    for(unsigned int i = 0; i < design->GetNumberOfVariables(); i++){
      *out << "\t" << d[i];
    }
  }
  
  if(logDesignGradient){
    StdVector<double> d;
    d.Resize(design->GetNumberOfVariables());
    design->WriteGradientToExtern(d.GetPointer(), DesignElement::COST_GRADIENT, DesignElement::PLAIN, NULL, false);
    for(unsigned int i = 0; i < design->GetNumberOfVariables(); i++){
      *out << "\t" << d[i];
    }
  }

  out->flush();
}

Optimization::MultipleExcitation::MultipleExcitation(bool multiple, ParamNode* pn)
{
  // set defaults
  this->stride = 1;
  this->damping = 1.0;
  this->max_gain = 1e4;
  this->meta_objective_ = false;
  this->multiple_excitation_ = multiple;

  // if disabled, we don't read anything
  if(!multiple || pn == NULL) return;

  // adjust defaults
  if(pn->Has("adjustWeights")) {
    this->stride   = pn->Get("adjustWeights")->Get("stride")->AsInt();
    this->max_gain = pn->Get("adjustWeights")->Get("max_gain")->AsDouble();
    this->damping  = pn->Get("adjustWeights")->Get("damping")->AsDouble();
  }

  // up to now we don't implement something else
  this->meta_objective_ = pn->Get("type")->AsString() == "meta_objective";
}

void Optimization::MultipleExcitation::ToInfo(InfoNode* in) const
{
  if(!multiple_excitation_) return;

  if(meta_objective_)  {
    InfoNode* in_ = in->Get("metaObjective");
    in_->Get("type")->SetValue("best_value");
    in_->Get("damping")->SetValue(damping);
    in_->Get("stride")->SetValue(stride);
    in_->Get("max_gain")->SetValue(max_gain);
  } else
    in->Get("type")->SetValue("fixed_weights");
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

double Excitation::GetFactor(Optimization::Objective* cost)
{
  double factor = 1.0; // default
  
  if(cost->FactorOmegaOmega())
  {
    if(frequency < 0.0) EXCEPTION("No frequency given");
    factor *= 4 * M_PI * M_PI * frequency * frequency;
  }
  
  return factor;
}

double Excitation::GetWeightedFactor(Optimization::Objective* cost)
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

void Excitation::ReadLoads(ParamNode* ls){
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
