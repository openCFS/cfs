#include "Optimization/Optimization.hh"
#include "Optimization/BaseOptimizer.hh"
#include "Optimization/DesignSpace.hh"
#include "General/exception.hh"
#include "Utils/StdVector.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "Utils/Timer.hh"

#include <sstream>
#include <cmath>
#include <fstream>

using namespace CoupledField;
using std::abs;

DECLARE_LOG(optimizer)
DEFINE_LOG(optimizer, "optimizer")

BaseOptimizer::Scale::Scale(BaseOptimizer* base, PtrParamNode autoscale, double manual_scale, bool no_autoscale)
 : target(0.0),
   tol(0.0),
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
    
    if(target == 0.0) 
      throw Exception("A target of 0.0 disabled autoscaling");
    if(manual_scale != 0.0 && manual_scale != 1.0) 
      throw Exception("Don't give explicit objective scaling with autoscale");

    tol    = autoscale->Get("tolerance")->As<Double>();
    // Cannot CalcAutoscale() has to be done in PostInit()
  }
  else
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
  if(target == 0.0) return;
  
  LOG_DBG(optimizer) << "Scale:CalcAutoscale() enter";
  
  // save the current tolerance and temporarily set to to 0 so we can call common methods
  double tol_save = tol;
  tol = 0.0;

  // We need a double design and gradient array. We use temporay ones  
 
  // evalue with the current (initial) design. Use this temporary gradient space
  StdVector<double> grad(base_->optimization->GetDesign()->GetNumberOfVariables());
  // make a temporary design as a copy from design space to copy it back there :)
  StdVector<double> data(grad.GetSize());
  // copy the design to our temporary space
  int design_id = base_->optimization->GetDesign()->WriteDesignToExtern(data.GetPointer());

  // evaluate the the gradient -> will be cheap in restart case
  bool good = base_->EvalGradObjective(grad.GetSize(), data.GetPointer(), false, grad);
  if(!good) EXCEPTION("internal error"); // needs to be good as tol = set to 0.0;
  assert(opt_scaling.value != 0.0);
  
  // our new scaling is the optimal scaling for now!
  scaling.design_id = opt_scaling.design_id;
  scaling.value = opt_scaling.value;
  
  // reset the tolerance
  tol = tol_save;
  autoscale_ = true;

  LOG_TRACE(optimizer) << "Scale::CalcAutoscale(): scale=" << scaling.value << " desing=" 
                       << scaling.design_id << " needed_eval=" << (opt_scaling.design_id != design_id);
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

BaseOptimizer::BaseOptimizer(Optimization* opt, PtrParamNode pn) :
  optimization(opt),
  info_(info->Get("optimization")->Get("optimizer")),
  objective(NULL),
  restart_requested(false),
  design_(DesignMemory(-1, 0.0)),
  optimizer_pn_(pn)
{
  this->timer_ =  new Timer();
  info_->Get(ParamNode::SUMMARY)->Get("timer")->SetValue(this->timer_ );
}

BaseOptimizer::~BaseOptimizer()
{ 
  if(objective != NULL) { delete objective; objective = NULL; }
}

void BaseOptimizer::PostInit(double manual_scaling, bool no_autoscale)
{
  PtrParamNode as = optimizer_pn_->Get("autoscale", ParamNode::PASS);
  objective = new Scale(this, as, manual_scaling, no_autoscale);
  objective->PostInit();
}

void BaseOptimizer::SolveOptimizationProblem()
{
  timer_->Start();
  SolveProblem();
  timer_->Stop();
}

std::string BaseOptimizer::LogFileHeader()
{
  std::string tmp = objective->target != 0.0 ? "\tmax_f_grad\tscale\topt_scale" : "\tmax_f_grad";
  if(GetCurrenActiveSetSize() >= 0)
    tmp += "\tactive_set";

  return tmp;
}


void BaseOptimizer::LogFileLine(std::ofstream* out, PtrParamNode iteration)
{
  if(out) *out << "\t" << objective->current.value;
  
  iteration->Get("max_f_grad")->SetValue(objective->current.value);

  if(objective->target != 0.0)
  {
    if(out) *out << "\t" << objective->scaling.value
                 << "\t" << objective->opt_scaling.value;
    
    iteration->Get("scale")->SetValue(objective->scaling.value);
    iteration->Get("opt_scale")->SetValue(objective->opt_scaling.value);
    if(GetCurrenActiveSetSize() >= 0)
      iteration->Get("active_set")->SetValue(GetCurrenActiveSetSize());
  }
}

double BaseOptimizer::EvalObjective(int n, const double* x, bool cfs_scale)
{
  assert(optimization->GetDesign()->GetNumberOfVariables() == (unsigned int) n);

  timer_->Stop();

  // set the design and see if it is a new one
  int new_design = optimization->GetDesign()->ReadDesignFromExtern(x);

  bool need_eval;
  
  if(new_design != design_.design_id)
  {
    design_.design_id = new_design;
    need_eval = true;

    // does a lot of work.
    optimization->SolveStateProblem();

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

  timer_->Start();

  double ret = cfs_scale ? objective->scaling.value * design_.value : design_.value;

  LOG_DBG(optimizer) << "EvalObjective: x_avg=" << Average(x, n) 
                     << " std_dev=" << StandardDeviation(x, n) 
                     << " is_new=" << need_eval << " -> " << design_.value
                     << " scaled=" << ret;
  return ret;
}

bool BaseOptimizer::EvalGradObjective(int n, const double* x, bool cfs_scale, StdVector<double>& grad_f)
{
  assert(optimization->GetDesign()->GetNumberOfVariables() == (unsigned int) n);

  timer_->Start();

  // set the design and see if it is a new one
  int new_design = optimization->GetDesign()->ReadDesignFromExtern(x);
  
  LOG_DBG2(optimizer) << "EvalGradObjective() external design_id=" << new_design << " design.design_id=" << design_.design_id;
  
  // do we have a valid value?
  bool need_eval = true;

  // new evaluation necessary?
  if(new_design != design_.design_id)
  {
    need_eval = true;

    // does a lot of work.
    optimization->SolveStateProblem();

    // set a proper objective - should be cheap and might be read from EvalObjective()!
    design_.value = optimization->CalcObjective();
    design_.design_id = new_design;
  }    
  else
  {
    need_eval = false;
    // only done in debug
    assert(close(design_.value, optimization->CalcObjective()));
  }
  
  if(new_design != design_.gradient_design_id){ // the adjoints have to be recalculated
    optimization->SolveAdjointProblems();
    design_.gradient_design_id = new_design;
  }
  
  LOG_DBG2(optimizer) << "EvalGradObjective: call CalcObjectiveGradient()";
  // calc our gradient - it is not stored anywhere
    
  assert(n <= (int) grad_f.GetSize()); // FIXME
  grad_f.window.Set(grad_f);
  optimization->CalcObjectiveGradient(&grad_f);

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
  if(cfs_scale)
    for(int i = 0; i < n; i++)
      grad_f[i] *= objective->scaling.value;

  timer_->Stop();

  LOG_DBG(optimizer) << "EvalGradObjective: cfs_scale=" << cfs_scale << " design=" << design_.design_id
                     << " need_eval=" << need_eval << " -> grad.scale=" << objective->scaling.value
                     << " grad.opt_scaling=" << objective->opt_scaling.value;
  
  return !restart_requested;
}

void BaseOptimizer::EvalConstraints(int n, const double* x, int m, bool cfs_scale, double* g_val)
{
  assert(m == optimization->constraints.view->GetNumberOfActiveConstraints());
  timer_->Start();

  // we overwrite the design space, but we do this all the time - especially before eval_f
  optimization->GetDesign()->ReadDesignFromExtern(x); 
     
  // iterate over all constraints
  for(int i = 0; i < m; i++)
  {
    Condition* g = optimization->constraints.view->Get(i);
    double org = optimization->CalcConstraint(g);
    double manual_scaling = 1.0;
    double objective_scaling = 1.0;
    if(cfs_scale)
    {
      if(g->DoObjectiveScaling())
        objective_scaling = objective->scaling.value;
      else
        manual_scaling = g->manual_scaling_value;
    }
    double val = org * manual_scaling * objective_scaling;

    LOG_DBG2(optimizer) << "EvalConstraints: g=" << g->type.ToString(g->GetType()) << " org=" << org
                        << " manual_scale=" << manual_scaling << " objective_scale=" << objective_scaling
                        << " ->" << val;
    g_val[i] = val;
  }
  optimization->constraints.view->Done(); // reset slope constraint to global mode

  timer_->Stop();
}

int BaseOptimizer::EvalGradConstraints(Condition* g, int start, bool cfs_scale,
    StdVector<double>& values, GradientType grtype)
{
  // always initialize the window!!
  int nnz(0); 
  values.window.Set(start, nnz);
  
  if(grtype == ALL || (grtype == LINEAR && g->IsLinear()) || (grtype == NONLINEAR && !g->IsLinear()))
  {
    nnz = g->GetSparsityPattern().GetSize();
    values.window.Set(start, nnz);
    
    // evaluate
    optimization->CalcConstraintGradient(g, &values);

    // we might ignore that value
    double scaling = g->DoObjectiveScaling() ? objective->scaling.value : g->manual_scaling_value;

    if(cfs_scale)
    {
      for(int p = 0; p < nnz; p++)
      {
        values[start + p] *= scaling;
        assert(values.InWindow(start + p));
      }
    }
  }
  return nnz;
}

void BaseOptimizer::EvalGradConstraints(int n, const double* x, int m, int nentries, bool cfs_scale,
      StdVector<double>& values, GradientType grtype)
{
  timer_->Start();

  // note, that we have dense gradients!
  // iterate over the gradients

  // set the window
  int start = 0;
  
  // reset values of the constraint gradients before the loop
  // as it also contains a loop over all the design elements
  optimization->GetDesign()->Reset(DesignElement::CONSTRAINT_GRADIENT, DesignElement::DEFAULT);
  
  for(int c = 0; c < m; c++)
  {
    int tmp = EvalGradConstraints(optimization->constraints.view->Get(c), start, cfs_scale, values, grtype); 
    start += tmp;
    LOG_DBG3(optimizer) << "EvalGradConstraints: co=" << c << " scaled val=" << values.ToString(true);
  }
  optimization->constraints.view->Done(); // reset slope constraint to global mode

  assert(start == nentries);

  timer_->Stop();
}

void BaseOptimizer::GetBounds(int n, double* x_l, double* x_u, int m, double* g_l, double* g_u){
  
  assert(n == (int) optimization->GetDesign()->GetNumberOfVariables());
  
  optimization->GetDesign()->WriteBoundsToExtern(x_l,x_u);

  assert(m == (int) optimization->constraints.view->GetNumberOfActiveConstraints());
    
  // normalization to =0 and <=0 constraints is done SCPIPBase   
    
  for(int i = 0; i < m; i++)
  {
    Condition* g = optimization->constraints.view->Get(i);
    // handle as in IPOPT 
    g_l[i] = g_u[i] = g->GetBoundValue();
    if(g->GetType() == Condition::SLOPE)
      g_l[i] *= -1.0;
    else
    {
      if(g->GetBound() == Condition::LOWER_BOUND) g_u[i] =  GetInfBound();
      if(g->GetBound() == Condition::UPPER_BOUND) g_l[i] = -GetInfBound();
    }
  }
  optimization->constraints.view->Done(); // reset slope constraint to global mode
  
}
