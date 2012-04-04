#include <assert.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/resultHandler.hh"
#include "Domain/domain.hh"
#include "Driver/assemble.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Function.hh"
#include "Optimization/Objective.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "PDE/basePDE.hh"
#include "Utils/StdVector.hh"
#include "Utils/Timer.hh"
#include "Utils/tools.hh"

using namespace CoupledField;
using std::abs;

DECLARE_LOG(optimizer)
DEFINE_LOG(optimizer, "optimizer")

BaseOptimizer::Scale::Scale(BaseOptimizer* base, PtrParamNode autoscale, double manual_scale, bool no_autoscale)
 : target(0.0),
   tol(0.0),
   logscale(false),
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
  if(target == 0.0)
  {
    autoscale_ = false;
    return;
  }

  LOG_DBG(optimizer) << "Scale:CalcAutoscale() enter";
  
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

BaseOptimizer::BaseOptimizer(Optimization* opt, PtrParamNode pn, Optimization::Optimizer type) :
  optimization(opt),
  type_(type),
  info_(info->Get("optimization")->Get("optimizer")),
  objective(NULL),
  restart_requested(false),
  timer_(new Timer()),
  design_(DesignMemory(-1, 0.0))
{
  optimizer_pn_ =  pn->Get(Optimization::optimizer.ToString(optimization->GetOptimizerType()), ParamNode::PASS);

  order.SetName("BaseOptimizer::Order");
  order.Add(BY_DESIGN, "by_design");
  order.Add(BY_ELEMENT, "by_element");

  // snopt, scpip and feasscp
  order_ = (optimizer_pn_ != NULL && optimizer_pn_->Has("order")) ? order.Parse(optimizer_pn_->Get("order")->As<std::string>()) : BY_DESIGN;
  SetupOrderMap(order_);

  LOG_DBG(optimizer) << "BO: optimizer_pn_=" << (optimizer_pn_ == NULL ? "null" : optimizer_pn_->GetName());
  LOG_DBG(optimizer) << "BO: order= " << order.ToString(order_);

  info_->Get(ParamNode::SUMMARY)->Get("timer")->SetValue(this->timer_ );
}

BaseOptimizer::~BaseOptimizer()
{ 
  if(objective != NULL) { delete objective; objective = NULL; }
}

void BaseOptimizer::PostInitScale(double manual_scaling, bool no_autoscale)
{
  PtrParamNode as = optimizer_pn_->Get("autoscale", ParamNode::PASS);
  objective = new Scale(this, as, manual_scaling, no_autoscale);
  objective->PostInit();
}

void BaseOptimizer::SolveOptimizationProblem()
{
  timer_->Start();
  SolveProblem();

  // dirty fix to have the final status streamed for iTop
  if(domain->GetResultHandler()->GetOutputWriter("streaming", true) != NULL && this->type_ != Optimization::EVALUATE_INITIAL_DESIGN)
    optimization->CommitIteration(true);

  timer_->Stop();
}

std::string BaseOptimizer::LogFileHeader()
{
  std::string tmp = objective->target != 0.0 ? "\tmax_f_grad\tscale\topt_scale" : "\tmax_f_grad";
  if(GetCurrenActiveSetSize() >= 0)
    tmp += "\tactive_set";

  // give the norm of the nonlinear gradients (ignore blown up)
  for(unsigned int i = 0; i < optimization->constraints.all.GetSize(); i++)
  {
    Condition* g = optimization->constraints.all[i];
    if(g->IsLinear()) continue;
    tmp += "\tmax_" + Function::type.ToString(g->GetType()) + "_grad";
  }

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

  // give the norm of the nonlinear gradients (ignore blown up)
  for(unsigned int i = 0; i < optimization->constraints.all.GetSize(); i++)
  {
    Condition* g = optimization->constraints.all[i];
    if(g->IsLinear()) continue;

    double mv = 0.0;
    const StdVector<DesignElement>& data = optimization->GetDesign()->data;

    for(int i = 0, n = data.GetSize(); i < n; i++)
      mv = std::max(mv, abs(data[i].GetValue(DesignElement::CONSTRAINT_GRADIENT, DesignElement::SMART, g)));

    iteration->Get("max_" + Function::type.ToString(g->GetType()) + "_grad")->SetValue(mv);
    if(out) *out << "\t" << mv;
  }

}

void BaseOptimizer::SetupOrderMap(Order order)
{
  if(!order_map.IsEmpty()) return;

  // we have to handle the slackvariable!
  const DesignSpace* space = optimization->GetDesign();

  order_map.Reserve(space->GetNumberOfVariables());
  assert(space->GetNumberOfVariables() == space->elements * space->design.GetSize());

  // cfs has design as outer loop, we need design as outer loop
  for(unsigned int d = 0; d < space->design.GetSize(); d++)
    for(unsigned int e = 0; e < space->elements; e++)
      order_map.Push_back(d + (space->design.GetSize() * e));

  LOG_DBG3(optimizer) << "SDEO: -> " << order_map.ToString();
}


double BaseOptimizer::EvalObjective(int n, const double* x, bool cfs_scale)
{
  assert(optimization->GetDesign()->GetNumberOfVariables() == (unsigned int) n);

  timer_->Stop();

  // set the design and see if it is a new one
  int new_design = optimization->GetDesign()->ReadDesignFromExtern(x);
  LOG_DBG(optimizer) << " set new design: avg " <<  Average(x, n)  << " std_dev = " << StandardDeviation(x, n) << " -> " << new_design;

  bool need_eval;
  
  if(new_design != design_.design_id)
  {
    design_.design_id = new_design;
    need_eval = true;
    
    // tell assemble, the design has changed
    domain->GetBasePDE()->getPDE_assemble()->SetAllReassemble();    

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

  if(objective->logscale && design_.value <= 0.0) EXCEPTION("Cannot do logarithmic autoscale for objective value " << design_.value);

  double sov = objective->logscale ? std::log(design_.value) : design_.value;

  double ret = cfs_scale ? objective->scaling.value * sov : sov;

  LOG_DBG(optimizer) << "EvalObjective: x_avg=" << Average(x, n) 
                     << " std_dev=" << StandardDeviation(x, n) 
                     << " is_new=" << need_eval << " -> "
                     << " ov=" << design_.value << " sov=" << sov
                     << " scaled=" << ret;
  
  LOG_DBG3(optimizer) << "x=" << StdVector<double>::ToString(n, x);

  timer_->Start();
  
  return ret;
}

bool BaseOptimizer::SolveAdjointProblemsIfNeeded(int n, const double* x, bool cfs_scale){
  // The function has to be evaluated before the gradient can be computed
  // This is true most times, as usually the rhs of the adjoint problem depends on the solution
  // On the other hand, it is called be the Optimizer before usually and so generates no cost. (But one cannot rely on this behaviour.)
  EvalObjective(n, x, cfs_scale);

  timer_->Stop();
  
  bool need_eval = design_.design_id != design_.gradient_design_id; 
  
  if(need_eval){ // the adjoints have to be recalculated
    optimization->SolveAdjointProblems();
    design_.gradient_design_id = design_.design_id;
  }
  
  timer_->Start();
  
  return(need_eval);  
}

bool BaseOptimizer::EvalGradObjective(int n, const double* x, bool cfs_scale, StdVector<double>& grad_f)
{
  bool need_eval = SolveAdjointProblemsIfNeeded(n, x, cfs_scale);
  
  timer_->Stop();
  
  LOG_DBG2(optimizer) << "EvalGradObjective: call CalcObjectiveGradient()";
  // calc our gradient - it is not stored anywhere
    
  assert(n <= (int) grad_f.GetSize()); // FIXME
  grad_f.window.Set(grad_f);
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

  timer_->Start();
  
  return !restart_requested;
}

void BaseOptimizer::EvalConstraints(int n, const double* x, int m, bool cfs_scale, double* g_val)
{
  assert(m == optimization->constraints.view->GetNumberOfActiveConstraints());

  // Before the constraints can be calculated it might be the case, that the forward problem needs recalculation
  // if it does not, this does not cost more than reading the design
  EvalObjective(n, x, cfs_scale);
  
  timer_->Stop();
  
  // iterate over all constraints
  for(int i = 0; i < m; i++)
  {
    Condition* g = optimization->constraints.view->Get(i);

    // do a complicated detection of local conditions handle the Local::active counter for logging
    if(g->IsLocalCondition())
    {
      assert(dynamic_cast<LocalCondition*>(g) != NULL);
      assert(g->GetLocal() != NULL);

      // reset the active counter for the first element
      if(static_cast<LocalCondition*>(g)->GetCurrentRelativePosition() == 0)
        g->GetLocal()->infeasible = 0;

      if(!g->IsFeasible())
        g->GetLocal()->infeasible++;
    }

    double manual_scaling = 1.0;
    double objective_scaling = 1.0;
    if(cfs_scale)
    {
      if(g->DoObjectiveScaling())
        objective_scaling = objective->scaling.value;
      else
        manual_scaling = g->manual_scaling_value;
    }
    double org = optimization->CalcConstraint(g);
    double base = org - (g->HasSlackBound() ? optimization->GetDesign()->GetSlackVariable() : 0.0);
    double val = base * manual_scaling * objective_scaling;

    LOG_DBG2(optimizer) << "EvalConstraints: g=" << g->type.ToString(g->GetType()) << " org=" << org
                        << " slack_corrected=" << base << " manual_scale=" << manual_scaling
                        << " objective_scale=" << objective_scaling << " ->" << val;
    g_val[i] = val;
  }
  optimization->constraints.view->Done(); // reset local constraint to global mode

  timer_->Start();
}

int BaseOptimizer::EvalGradConstraints(Condition* g, int start, bool cfs_scale,
    StdVector<double>& values, GradientType grtype)
{
  timer_->Stop();
  
  // always initialize the window!!
  int nnz(0); 
  values.window.Set(start, nnz);
  
  LOG_DBG(optimizer) << "EGC g=" << g->ToString(NULL) << " start=" << start << " gt=" << grtype << " nnz=" << g->GetSparsityPattern().GetSize();

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
  
  timer_->Start();
  
  return nnz;
}

void BaseOptimizer::EvalGradConstraints(int n, const double* x, int m, int nentries, bool cfs_scale,
      StdVector<double>& values, GradientType grtype)
{
  SolveAdjointProblemsIfNeeded(n, x, cfs_scale);
  
  timer_->Stop();

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

  timer_->Start();
}

void BaseOptimizer::ReorderDesign(int n, double* x, bool reverse)
{
  assert((int) order_map.GetSize() == n);

  StdVector<double> tmp(n);

  for(int i = 0; i < n; i++)
    reverse ? tmp[i] = x[order_map[i]] : tmp[order_map[i]] = x[i];

  for(int i = 0; i < n; i++)
    x[i] = tmp[i];
}

void BaseOptimizer::GetBounds(int n, double* x_l, double* x_u, int m, double* g_l, double* g_u)
{
  assert(n == (int) optimization->GetDesign()->GetNumberOfVariables());

  timer_->Stop();
  
  optimization->GetDesign()->WriteBoundsToExtern(x_l,x_u);

  assert(m == (int) optimization->constraints.view->GetNumberOfActiveConstraints());
    
  // normalization to =0 and <=0 constraints is done SCPIPBase   
    
  for(int i = 0; i < m; i++)
  {
    Condition* g = optimization->constraints.view->Get(i);
    // handle as in IPOPT 
    g_l[i] = g_u[i] = g->GetBoundValue();
    
    // Snopt alone is able to handle bounding boxes for constraints. This needs to
    // be reflected when the number of constraints is determined. Here Function::Local::Local()
    // only when we have NO box box constraints we need NEXT_AND_REVERSE
    if(g->GetType() == Condition::SLOPE && g->GetLocal()->GetLocality() == Function::Local::NEXT)
      g_l[i] *= -1.0;
    else
    {
      if(g->GetBound() == Condition::LOWER_BOUND) g_u[i] =  GetInfBound();
      if(g->GetBound() == Condition::UPPER_BOUND) g_l[i] = -GetInfBound();
    }
  }
  optimization->constraints.view->Done(); // reset slope constraint to global mode
  
  timer_->Start();
}
