#include <assert.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ResultHandler.hh"
#include "Domain/Domain.hh"
#include "Driver/Assemble.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Function.hh"
#include "Optimization/Objective.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "Optimization/Tune.hh"
#include "PDE/BasePDE.hh"
#include "Utils/StdVector.hh"
#include "Utils/Timer.hh"
#include "Utils/tools.hh"
#include "Utils/PythonKernel.hh"

using namespace CoupledField;
using std::abs;

DEFINE_LOG(optimizer, "optimizer")

const string BaseOptimizer::Tuned::MOVE_LIMIT = "move_limit";

BaseOptimizer::Scale::Scale(BaseOptimizer* base, PtrParamNode autoscale, double manual_scale, bool no_autoscale)
 : target(0.0),
   tol(0.0),
   logscale(false),
   manual(0.0),
   opt_scaling(DesignMemory(-1, 0.0)),
   scaling(DesignMemory(-1, 0.0)),
   current(DesignMemory(-1, 0.0)),
   autoscale_(false),
   base_(base)
{
  // the vectors are zero size by default!
  if(autoscale != NULL && !no_autoscale)
  {
    target = autoscale->Get("target")->As<Double>();
    
    if(manual_scale != 0.0 && manual_scale != 1.0) 
      throw Exception("Don't give explicit objective scaling with autoscale");

    tol    = autoscale->Get("tolerance")->As<Double>();

    logscale = autoscale->Get("logscale")->As<bool>();

    manual   = autoscale->Get("manual")->As<Double>();

    // Cannot CalcAutoscale() has to be done in PostInit()
  }
  if(target == 0.0) // there might have been an autoscale element with 0 target
  {
    this->scaling.value = manual_scale;
    this->scaling.value *= base_->optimization->objectives.DoMaximize() ? -1.0 : 1.0;
  }
  LOG_TRACE(optimizer) << "Scale::Scale() target=" << target << " tol=" << tol << " manual_scale=" 
                       << manual_scale << " -> PostInit() pending!";
}

void BaseOptimizer::Scale::PostInit()
{
  // calculate the autoscale values (more or less expensive, never cheap)
  CalcAutoscale();

  LOG_TRACE(optimizer) << "Scale::PostInit() target=" << target << " tol=" << tol << " -> scale=" << scaling.value;
}


void BaseOptimizer::Scale::CalcAutoscale()
{
  if(target == 0.0 && manual == 0.0)
  {
    autoscale_ = false;
    return;
  }

  LOG_DBG(optimizer) << "Scale:CalcAutoscale() enter";
  if (manual == 0.0) {
    // save the current tolerance and temporarily set to to 0 so we can call common methods
    double tol_save = tol;
    tol = 0.0;
  
    // We need a double design and gradient array. We use temporary ones

    // evaluate with the current (initial) design. Use this temporary gradient space
    StdVector<double> grad(base_->optimization->GetDesign()->GetNumberOfVariables());
    // make a temporary design as a copy from design space to copy it back there :) - Fabian: what is copied back??
    StdVector<double> data(grad.GetSize());
    // copy the design to our temporary space
    int design_id = base_->optimization->GetDesign()->WriteDesignToExtern(data.GetPointer());
  
    // evaluate the the gradient -> will be cheap in restart case
    bool good = base_->EvalGradObjective(grad.GetSize(), data.GetPointer(), false, grad);
    if(!good) EXCEPTION("internal error"); // needs to be good as tol = set to 0.0;
    assert(opt_scaling.value != 0.0);
    // reset the tolerance
    tol = tol_save;
    // our new scaling is the optimal scaling for now!
    scaling.design_id = opt_scaling.design_id;
    scaling.value = opt_scaling.value;
    LOG_TRACE(optimizer) << "Scale::CalcAutoscale(): scale=" << scaling.value << " design="
                         << scaling.design_id << " needed_eval=" << (opt_scaling.design_id != design_id);
  } else {
    scaling.design_id = opt_scaling.design_id;
    scaling.value = manual;
    opt_scaling.value = manual;
    LOG_TRACE(optimizer) << "Scale::CalcAutoscale(): manual scaling. Scale=" << scaling.value << " design="
                         << scaling.design_id;
  }
  assert(opt_scaling.value != 0.0);
}


/** Checks the current gradients in the design space against the tolerance.
 * @return true if not active or no tolerance or within tolerance */
bool BaseOptimizer::Scale::CheckScaling(int n, StdVector<double>& grad)
{
  assert(base_->design_.design_id != -1);

  current.design_id     = base_->design_.design_id;

  // analyse the gradient. Is always an plot output
  current.value = 0.0;  
  for(int i = 0; i < n; i++)
    current.value = std::max(current.value, abs(grad[i]));
  
  // we cannot do anything w/o target
  if(target == 0.0) return true;
  
  // first calc the optimal scaling
  if(base_->design_.design_id != opt_scaling.design_id)
  {
    // we assume the current settings are correct!
    opt_scaling.design_id = base_->design_.design_id;

    opt_scaling.value = target / current.value; 
    opt_scaling.value *= base_->optimization->objectives.DoMaximize() ? -1.0 : 1.0;
    
    LOG_TRACE(optimizer) << "Scale::check::set opt scaling: design=" << opt_scaling.design_id << " target=" << target 
                         << " opt_scaling=" << opt_scaling.value << " maximize="
                         << base_->optimization->objectives.DoMaximize();
  }

  assert(base_->design_.design_id == current.design_id);
  // now we have a valid optimal scaling
  double act = abs(current.value * scaling.value);
  bool good = act <= abs(target) * (1.0 + tol) && act >= abs(target) * (1.0 - tol);

  // the total result:
  // do we have to check the current scaling for rescaling? If not it is alwys good by definition :)!
  bool result = tol != 0.0 ? good : true; 
  
  LOG_DBG2(optimizer) << "Scale::check tolerance:"
                      << " design: " << base_->design_.design_id << "/" << base_->design_.value 
                      << " scaling: " << scaling.design_id << "/" << scaling.value 
                      << " opt_scaling: " << opt_scaling.design_id << "/" << opt_scaling.value 
                      << " current: " << current.design_id << "/" << current.value
                      << ": " << act << " in [" << (abs(target) * (1.0 - tol)) << ";" 
                      << (abs(target) * (1.0 + tol)) << "] -> " << good << " ->> " << result;
  return result;
}


std::string BaseOptimizer::Scale::ToString()
{
  // scaling.value = -1.0 when no target but maximization
  if(scaling.value == 0.0 || target == 0.0)
  {
    return "No scaling";
  }

  std::ostringstream os;
  
  if(target != 0.0)
  {
    os << "Autosale target=" << target;
    if(tol != 0.0) os << " tolerance=" << tol;
    os << " scaling=" << scaling.value;
  }
  
  return os.str();
}


BaseOptimizer::Tuned::Tuned(PtrParamNode pn, double* value, double max, double divider, BaseOptimizer* base)
{
  assert(value != nullptr);
  this->value = value;
  this->max = max; // has a default
  this->divider = divider;
  if(pn)
  {
    pn->GetValue("max", this->max, ParamNode::PASS);
    pn->GetValue("min", this->min, ParamNode::PASS);
    pn->GetValue("transition_divider", this->divider, ParamNode::PASS);
  }
  for(Tune* tune : base->optimization->tunes)
    if(tune->GetUsage() == tune->BETA)
      this->tune = tune;
  if(this->tune == nullptr)
    throw Exception("'tuned' move_limit required 'tune' for the density projection beta value.");

  this->base = base;
  base->tuned = this; // to have Update called within CommitIteration
  base->optimization->RegisterAuxLogValue("move_limit", 0);
}

void BaseOptimizer::Tuned::ToInfo(PtrParamNode in)
{
  in->Get("tune")->SetValue(tune->usage.ToString(tune->BETA));
  in->Get("max")->SetValue(max);
  in->Get("min")->SetValue(min);
  in->Get("divider")->SetValue(divider);
}

void BaseOptimizer::Tuned::Update()
{
  double start = this->tune->CalcTransistionZone(.1);
  double end = this->tune->CalcTransistionZone(1-.1);
  double tz = end-start;
  *value = std::min(tz/divider, max);
  base->optimization->Optimization::SetAuxLogValue("move_limit", *value);
}


BaseOptimizer::BaseOptimizer(Optimization* opt, PtrParamNode pn, Optimization::Optimizer type) :
  optimization(opt),
  type_(type),
  info_(opt->optInfoNode->Get("optimizer")->Get(Optimization::optimizer.ToString(type))),
  objective(NULL),
  restart_requested(false),
  design_(DesignMemory(-1, 0.0))
{
  assert(pn != NULL);
  gen_opt_pn_  =  pn;
  this_opt_pn_ =  pn->Get(Optimization::optimizer.ToString(optimization->GetOptimizerType()), ParamNode::PASS);

  LOG_DBG(optimizer) << "BO: gen_opt_pn_=" << gen_opt_pn_->GetName();
  LOG_DBG(optimizer) << "BO: this_opt_pn_=" << (this_opt_pn_ != NULL ? this_opt_pn_->GetName() : "null") ;

  // evaluate  timer makes no sense and is also not implemented in a clean way
  if(type == Optimization::EVALUATE_INITIAL_DESIGN)
    optimizer_timer_ = boost::shared_ptr<Timer>(new Timer());
  else
    optimizer_timer_ = info_->Get(ParamNode::SUMMARY)->Get(Optimization::optimizer.ToString(type))->Get("timer")->AsTimer();

  optimizer_timer_->SetNesting(false);
  optimizer_timer_->Start(); // stop after specialized constrctor or in PostInit implementation

  eval_obj_timer_        = optimization->optInfoNode->Get(ParamNode::SUMMARY)->Get("eval_objective/timer")->AsTimer();
  eval_const_timer_      = optimization->optInfoNode->Get(ParamNode::SUMMARY)->Get("eval_constraints/timer")->AsTimer();
  eval_grad_obj_timer_   = optimization->optInfoNode->Get(ParamNode::SUMMARY)->Get("eval_grad_objective/timer")->AsTimer();
  eval_grad_const_timer_ = optimization->optInfoNode->Get(ParamNode::SUMMARY)->Get("eval_grad_constraints/timer")->AsTimer();

  assert(!GetRunningEvalTimer());
  assert(GetRunningEvalTimer() == nullptr);
}

BaseOptimizer::~BaseOptimizer()
{ 
  if(objective) { delete objective; objective = nullptr; }
  if(tuned) { delete tuned; tuned = nullptr; }
}

void BaseOptimizer::PostInit()
{
  if(tuned)
    tuned->Update(); // log initial value
}

boost::shared_ptr<Timer> BaseOptimizer::GetRunningEvalTimer()
{
  // only one may run
  assert(!(eval_obj_timer_->IsRunning() && (eval_const_timer_->IsRunning() || eval_grad_obj_timer_->IsRunning() || eval_grad_const_timer_->IsRunning() )));
  assert(!(eval_const_timer_->IsRunning() && (eval_obj_timer_->IsRunning() || eval_grad_obj_timer_->IsRunning() || eval_grad_const_timer_->IsRunning() )));
  assert(!(eval_grad_obj_timer_->IsRunning() && (eval_const_timer_->IsRunning() || eval_obj_timer_->IsRunning() || eval_grad_const_timer_->IsRunning() )));
  assert(!(eval_grad_const_timer_->IsRunning() && (eval_const_timer_->IsRunning() || eval_grad_obj_timer_->IsRunning() || eval_obj_timer_->IsRunning() )));

  if(eval_obj_timer_->IsRunning())
    return eval_obj_timer_;
  if(eval_const_timer_->IsRunning())
    return eval_const_timer_;
  if(eval_grad_obj_timer_->IsRunning())
    return eval_grad_obj_timer_;
  if(eval_grad_const_timer_->IsRunning())
    return eval_grad_const_timer_;
  // Nothing is running when we do a direct eval constraint call
  //std::cout << "BO:GRET -> NULL\n";
  return boost::shared_ptr<Timer>(); // http://stackoverflow.com/questions/16229401/initialize-a-boostshared-ptr-to-null
}

bool BaseOptimizer::ValidateTimers()
{
  assert(!optimizer_timer_->IsRunning());
  assert(!eval_obj_timer_->IsRunning());
  assert(!eval_const_timer_->IsRunning());
  assert(!eval_grad_obj_timer_->IsRunning());
  assert(!eval_grad_const_timer_->IsRunning());
  return true;
}

void  BaseOptimizer::CommitIteration()
{
  bool restart = optimizer_timer_->IsRunning();
  if(restart)
    optimizer_timer_->Stop();
  assert(ValidateTimers());
  optimization->CommitIteration();
  if(restart)
    optimizer_timer_->Start();
  if(tuned)
    tuned->Update(); // read beta and write move_limit
}

void BaseOptimizer::PostInitScale(double manual_scaling, bool no_autoscale)
{
  PtrParamNode as = gen_opt_pn_->Get("autoscale", ParamNode::PASS);
  objective = new Scale(this, as, manual_scaling, no_autoscale);
  objective->PostInit();
}

void BaseOptimizer::SolveOptimizationProblem()
{
  optimizer_timer_->Start(); // we already might have spent time in the optimzer in Init()
  SolveProblem();
  if(optimizer_timer_->IsRunning())
    optimizer_timer_->Stop();
}

void BaseOptimizer::LogFileHeader(Optimization::Log& log)
{
  if(log.gradNorm)
    log.AddToHeader("max_f_grad");

  if(objective->target != 0.0) {
    log.AddToHeader("scale");
    log.AddToHeader("opt_scale");
  }

  // give the norm of the nonlinear gradients (ignore blown up)
  for(unsigned int i = 0; i < optimization->constraints.all.GetSize(); i++)
  {
    Condition* g = optimization->constraints.all[i];
    if(!g->IsLinear() && log.gradNorm)
      log.AddToHeader("max_" + g->ToString() + "_grad");
  }
}


void BaseOptimizer::LogFileLine(std::ofstream* out, PtrParamNode iteration)
{
  if(out && optimization->log.gradNorm)
    *out << " \t" << objective->current.value;
  
  if(optimization->log.gradNorm)
    iteration->Get("max_f_grad")->SetValue(objective->current.value);

  if(objective->target != 0.0)
  {
    if(out)
      *out << " \t" << objective->scaling.value
           << " \t" << objective->opt_scaling.value;

    iteration->Get("scale")->SetValue(objective->scaling.value);
    iteration->Get("opt_scale")->SetValue(objective->opt_scaling.value);
    if(GetCurrenActiveSetSize() >= 0)
      iteration->Get("active_set")->SetValue(GetCurrenActiveSetSize());
  }

  // give the norm of the nonlinear gradients (ignore blown up)
  for(unsigned int i = 0; i < optimization->constraints.all.GetSize(); i++)
  {
    Condition* g = optimization->constraints.all[i];
    if(!g->IsLinear() && this->optimization->log.gradNorm)
    {
      double mv = 0.0;
      const StdVector<DesignElement>& data = optimization->GetDesign()->data;

      for(int i = 0, n = data.GetSize(); i < n; i++)
        mv = std::max(mv, abs(data[i].GetValue(DesignElement::CONSTRAINT_GRADIENT, DesignElement::SMART, g)));

      iteration->Get("max_" + g->ToString() + "_grad")->SetValue(mv);
      if(out)
        *out << " \t" << mv;
    }
  }
}


double BaseOptimizer::EvalObjective(int n, const double* x, bool cfs_scale)
{

  assert(optimization->GetDesign()->GetNumberOfVariables() == (unsigned int) n);
  // we might come from another eval, then the optimizer is already stopped and we must not restart it
  bool restart_timer = optimizer_timer_->IsRunning();
  if(restart_timer)
    optimizer_timer_->Stop();

  python->CallHook(PythonKernel::OPT_EVAL_FUNC);

  assert(!GetRunningEvalTimer()); // no currently running timer!

  eval_obj_timer_->Start();

  // set the design and see if it is a new one
  int new_design = optimization->GetDesign()->ReadDesignFromExtern(x);
  LOG_DBG(optimizer) << " set new design: avg " <<  Average(x, n)  << " std_dev = " << StandardDeviation(x, n) << " -> " << new_design;
  bool need_eval;
  
  if(new_design != design_.design_id)
  {
    design_.design_id = new_design;
    need_eval = true;
    
    // TODO: handle case for non-design sensitive objective!

    // tell assemble, the design has changed
    for(unsigned int c = 0; c < optimization->manager.context.GetSize(); c++)
      optimization->manager.context[c].pde->GetAssemble()->SetAllReassemble();

    // does a lot of work.
    eval_obj_timer_->Stop();
    optimization->SolveStateProblem();
    eval_obj_timer_->Start();


    // take the fruit of the work
    design_.value = optimization->CalcObjective();
  }
  else
  {
    // the design did not change
    // in the debug case we evaluate the objective to be sure!
    assert(close(design_.value, optimization->CalcObjective()));
    need_eval = false;
  }

  if(objective->logscale && design_.value <= 0.0) EXCEPTION("Cannot do logarithmic autoscale for objective value " << design_.value);

  double sov = objective->logscale ? std::log(design_.value) : design_.value;

  double ret = cfs_scale ? objective->scaling.value * sov : sov;

  LOG_DBG(optimizer) << "EvalObjective: x_avg=" << Average(x, n) 
                     << " std_dev=" << StandardDeviation(x, n) 
                     << " is_new=" << need_eval << " -> "
                     << " ov=" << design_.value << " sov=" << sov
                     << " scaled=" << ret;
  
  LOG_DBG3(optimizer) << "x=" << StdVector<double>::ToString(n, x);

  eval_obj_timer_->Stop();
  if(restart_timer)
    optimizer_timer_->Start();

  return ret;
}

bool BaseOptimizer::SolveAdjointProblemsIfNeeded(int n, const double* x, bool cfs_scale)
{
  // we need to take care about measurement! We don't want to count times double!
  assert(!optimizer_timer_->IsRunning());
  // we eval objective next so there must not be anything running at the moment!
  assert(!eval_grad_obj_timer_->IsRunning());
  assert(!eval_grad_const_timer_->IsRunning());
  // The function has to be evaluated before the gradient can be computed
  // This is true most times, as usually the rhs of the adjoint problem depends on the solution
  // On the other hand, it is has usually been called by the Optimizer before and so generates no cost. (But one cannot rely on this behavior.)
  EvalObjective(n, x, cfs_scale);
  
  bool need_eval = design_.design_id != design_.gradient_design_id; 
  
  if(need_eval){ // the adjoints have to be recalculated

    optimization->SolveAdjointProblems();
    design_.gradient_design_id = design_.design_id;
  }

  return(need_eval);  
}

bool BaseOptimizer::EvalGradObjective(int n, const double* x, bool cfs_scale, StdVector<double>& grad_f)
{
//  assert(optimization->CalcObjectiveCalled());
  Optimization* opt = domain->GetOptimization();
  if(!opt->CalcObjectiveCalled())
    opt->CalcObjective();

  bool opt_run = optimizer_timer_->IsRunning(); // in the scale case we have to optimization_timer
  if(opt_run)
    optimizer_timer_->Stop();

  python->CallHook(PythonKernel::OPT_EVAL_GRAD);

  // triggers EvalObjective, so start timer afterwards
  bool need_eval = SolveAdjointProblemsIfNeeded(n, x, cfs_scale);

  eval_grad_obj_timer_->Start();

  LOG_DBG2(optimizer) << "EvalGradObjective: call CalcObjectiveGradient()";
  // calc our gradient - it is not stored anywhere
    
  assert(n <= (int) grad_f.GetSize()); // FIXME
  // grad_f.window.Set(grad_f);
  grad_f.window.Set(0, n);
  optimization->CalcObjectiveGradient(&grad_f);

  if(objective->logscale)
  {
    double ov = EvalObjective(n, x, false);
    for(int i = 0; i < n; i++)
      grad_f[i] /= ov;
    LOG_DBG2(optimizer) << "EGO: ov=" << ov;
  }

  // check the scaling - harmless if not restarted autoscale
  restart_requested = false;
  LOG_DBG2(optimizer) << "EvalGradObjective: call CheckScaling()";
  if(!objective->CheckScaling(n, grad_f))
  {
    LOG_DBG(optimizer) << "EvalGradObjective: scaling violated";
    assert(design_.design_id == objective->opt_scaling.design_id);
    restart_requested = true;
  }
  
  // eventually perform the scaling
  assert(!(cfs_scale && objective->scaling.value == 0.0));
  if(cfs_scale)
    for(int i = 0; i < n; i++)
      grad_f[i] *= objective->scaling.value;

  LOG_DBG(optimizer) << "EvalGradObjective: cfs_scale=" << cfs_scale << " design=" << design_.design_id
                     << " need_eval=" << need_eval << " -> grad.scale=" << objective->scaling.value
                     << " grad.opt_scaling=" << objective->opt_scaling.value;

  if(opt_run)
    optimizer_timer_->Start();
  eval_grad_obj_timer_->Stop();
  
  return !restart_requested;
}

void BaseOptimizer::EvalConstraints(int n, const double* x, int m, bool cfs_scale, double* g_val, bool normalize)
{
  optimizer_timer_->Stop();
  python->CallHook(PythonKernel::OPT_EVAL_FUNC);

  ConditionContainer& cc = optimization->constraints;

  assert(m == cc.view->GetNumberOfActiveConstraints());

  // Before the constraints can be calculated it might be the case, that the forward problem needs recalculation
  // if it does not, this does not cost more than reading the design
  if(cc.IsAllStateDependent())
    EvalObjective(n, x, cfs_scale);
  else
    optimization->GetDesign()->ReadDesignFromExtern(x);
  
  eval_const_timer_->Start(); // After EvalObjective();
  
  // iterate over all constraints
  for(int i = 0; i < m; i++)
  {
    Condition* g = cc.view->Get(i);

    double val = EvalConstraint(g, cfs_scale, normalize, false); // no direct call

    g_val[i] = val;
  }
  cc.view->Done(); // reset local constraint to global mode

  eval_const_timer_->Stop();
  optimizer_timer_->Start();
}

double BaseOptimizer::EvalConstraint(Condition* g, bool cfs_scale, bool normalize, bool direct_call, Excitation* ev_only_excite)
{
  // for a proper time measurement we have to know if we are called by EvalConstraints() or directly (FeasPP)
  assert(!(!direct_call && !eval_const_timer_->IsRunning()));
  if(direct_call) {
    assert(optimizer_timer_->IsRunning());
    assert(!GetRunningEvalTimer());
    optimizer_timer_->Stop();
    eval_const_timer_->Start();
  }
  assert(!optimizer_timer_->IsRunning());


  // do a complicated detection of local conditions handle the Local::active counter for logging
  double manual_scaling = 1.0;
  double objective_scaling = 1.0;
  if(cfs_scale)
  {
    if(g->DoObjectiveScaling())
      objective_scaling = objective->scaling.value;
    else
      manual_scaling = g->manual_scaling_value;
  }
  double org = optimization->CalcConstraint(g, ev_only_excite);
  double base = org;
  if(g->HasGeneralSlackBound())
  {
    if(g->IsGeneralSlackBound(Condition::SLACK_VALUE))
      base -= optimization->GetDesign()->GetSlackVariable();
    if(g->IsGeneralSlackBound(Condition::ALPHA_VALUE))
      base -= optimization->GetDesign()->GetAlphaVariable();
    else if(g->IsGeneralSlackBound(Condition::ALPHA_PLUS_SLACK_VALUE))
      base -= optimization->GetDesign()->GetAlphaVariable() + optimization->GetDesign()->GetSlackVariable();
    else if(g->IsGeneralSlackBound(Condition::ALPHA_MINUS_SLACK_VALUE))
      base -= optimization->GetDesign()->GetAlphaVariable() - optimization->GetDesign()->GetSlackVariable();
  }

  double scaled = base * manual_scaling * objective_scaling;
  double val = normalize ? (g->GetBound() == Condition::LOWER_BOUND ? -1.0 : 1.0) * (scaled - g->GetBoundValue()) : scaled;

  LOG_DBG2(optimizer) << "BO:EC g=" << g->type.ToString(g->GetType()) << " org=" << org
      << " base=" << base << " ms=" << manual_scaling
      << " os=" << objective_scaling << " scaled=" << scaled << " -> " << val;

  if(direct_call) {
    eval_const_timer_->Stop();
    optimizer_timer_->Start();
  }

  return val;
}



void BaseOptimizer::EvalGradConstraints(int n, const double* x, int m, int nentries, bool cfs_scale, bool normalize,
      StdVector<double>& values, GradientType grtype)
{
  // Attention! there is a copy and paste clone in FeasPP::SolveSubProblem()!
  optimizer_timer_->Stop();

  python->CallHook(PythonKernel::OPT_EVAL_GRAD);

  if(grtype != LINEAR)
    // triggers EvalObjective, so start timer afterwards
    SolveAdjointProblemsIfNeeded(n, x, cfs_scale);
  else
    EvalObjective(n, x, cfs_scale);

  eval_grad_const_timer_->Start();

  // note, that we have dense gradients!
  // iterate over the gradients

  // set the window
  int start = 0;
  
  // reset values of the constraint gradients before the loop
  // as it also contains a loop over all the design elements
  optimization->GetDesign()->Reset(DesignElement::CONSTRAINT_GRADIENT, DesignElement::DEFAULT);
  
  for(int c = 0; c < m; c++)
  {
    Condition* g = optimization->constraints.view->Get(c);

    if(grtype == ALL || (grtype == LINEAR && g->IsLinear()) || (grtype == NONLINEAR && !g->IsLinear()))
    {
      int tmp = EvalGradConstraint(optimization->constraints.view->Get(c), start, cfs_scale, normalize, values, false); // no direct call
      LOG_DBG3(optimizer) << "EGCs: g=" << g->ToString() << " co=" << c << " scaled val=" << values.ToString(true);
      start += tmp;
    }
  }
  optimization->constraints.view->Done(); // reset slope constraint to global mode

  assert(start == nentries);

  eval_grad_const_timer_->Stop();
  optimizer_timer_->Start();
}

int BaseOptimizer::EvalGradConstraint(Condition* g, int start, bool cfs_scale, bool normalize, StdVector<double>& values, bool direct_call)
{
  // for a proper time measurement we have to know if we are called by EvalGradConstraints() or directly (FeasPP)
  assert(!(!direct_call && !eval_grad_const_timer_->IsRunning()));
  if(direct_call) {
    assert(optimizer_timer_->IsRunning());
    assert(!GetRunningEvalTimer());
    optimizer_timer_->Stop();
    eval_grad_const_timer_->Start();
  }
  // always initialize the window!!
  int nnz = 0;
  values.window.Set(start, nnz);

  LOG_DBG(optimizer) << "EGC g=" << g->ToString() << " start=" << start << " nnz=" << g->GetSparsityPattern().GetSize();
  assert(g->GetSparsityPatternSize() == g->GetSparsityPattern().GetSize());
  nnz = g->GetSparsityPatternSize();
  values.window.Set(start, nnz);

  // evaluate
  optimization->CalcConstraintGradient(g, &values);

  // we might ignore that value
  double scaling = g->DoObjectiveScaling() ? objective->scaling.value : g->manual_scaling_value;

  for(int p = 0; cfs_scale && p < nnz; p++)
  {
    values[start + p] *= scaling;
    assert(values.InWindow(start + p));
  }

  // apply flip_sign only when we normalize
  double flip_sign = g->GetBound() == Condition::LOWER_BOUND ? -1.0 : 1.0;
  for(int p = 0; normalize && p < nnz; p++)
    values[start + p] *= flip_sign;

  if(direct_call) {
     eval_grad_const_timer_->Stop();
     optimizer_timer_->Start();
   }
  return nnz;
}

void BaseOptimizer::GetBounds(int n, double* x_l, double* x_u, int m, double* g_l, double* g_u)
{
  assert(n == (int) optimization->GetDesign()->GetNumberOfVariables());

  bool restart_timer = optimizer_timer_->IsRunning();
  if(restart_timer)
    optimizer_timer_->Stop(); // makes not much sense for EvaluateOnly!
  
  optimization->GetDesign()->WriteBoundsToExtern(x_l,x_u);

  GetConstraintsBounds(m, g_l, g_u);

  if(restart_timer)
    optimizer_timer_->Start();
}

void BaseOptimizer::GetConstraintsBounds(int m, double* g_l, double* g_u)
{
  assert(m == (int) optimization->constraints.view->GetNumberOfActiveConstraints());

  bool restart_timer = optimizer_timer_->IsRunning();
  if(restart_timer)
    optimizer_timer_->Stop(); // makes not much sense for EvaluateOnly!

  // normalization to =0 and <=0 constraints is done SCPIPBase   
  for(int i = 0; i < m; i++)
  {
    Condition* g = optimization->constraints.view->Get(i);

    // this is the equal case
    g_u[i] = g->GetBoundValue();
    g_l[i] = g_u[i];

    // Snopt alone is able to handle bounding boxes for constraints. This needs to
    // be reflected when the number of constraints is determined. Here Function::Local::Local()
    // only when we have NO box constraints we need NEXT_AND_REVERSE
    if(g->IsDoubleBounded() && g->GetBound() != Condition::EQUAL) {// checks the locality!!
      // the snopt solver was necessary to set to a non-REVERSE locality. This check needs to be done when the number of constraints is to
      // be determined.
      assert(type_ == Optimization::SNOPT_SOLVER);
      g_u[i] =  g->GetBoundValue();
      g_l[i] = -g->GetBoundValue();
    }
    else
    {
      if(g->GetBound() == Condition::LOWER_BOUND) {
        g_u[i] = GetInfBound();
        g_l[i] = g->GetBoundValue();
      }
      if(g->GetBound() == Condition::UPPER_BOUND) {
        g_u[i] = g->GetBoundValue();
        g_l[i] = -GetInfBound();
      }
    }
    LOG_DBG2(optimizer) << "BO::GB i=" << i << " g=" << g->ToString() << " bv=" << g->GetBoundValue() << " DB=" << g->IsDoubleBounded() << " l=" << g_l[i] << " u=" << g_u[i];
  }
  optimization->constraints.view->Done(); // reset slope constraint to global mode

  if(restart_timer)
    optimizer_timer_->Start();
}
