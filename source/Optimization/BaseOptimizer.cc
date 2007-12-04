#include "Optimization/Optimization.hh"
#include "Optimization/BaseOptimizer.hh"
#include "Optimization/DesignSpace.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/exception.hh"
#include "Utils/StdVector.hh"
#include "DataInOut/Logging/cfslog.hh"

#include <sstream>
#include <cmath>

using namespace CoupledField;
using std::abs;

DECLARE_LOG(optimizer)
DEFINE_LOG(optimizer, "optimizer")


BaseOptimizer::Scale::Scale(BaseOptimizer* base, ParamNode* autoscale, double manual_scale)
{
  this->base_ = base;
  this->target = 0.0;
  this->tol    = 0.0;
  
  this->opt_scaling.design_id = -1;
  this->opt_scaling.value     = 0.0;

  this->scaling.design_id = -1;
  this->scaling.value     = 0.0;
  
  this->current.design_id = -1;
  this->current.value     = 0.0;

  // the vectors are zero size by default!
  if(autoscale != NULL)
  {
    target = autoscale->Get("target")->AsDouble();
    
    if(target == 0.0) 
      throw Exception("A target of 0.0 disabled autoscaling");
    if(manual_scale != 0.0 && manual_scale != 1.0) 
      throw Exception("Don't give explicit objective scaling with autoscale");

    tol    = autoscale->Get("tolerance")->AsDouble();
    // Cannot CalcAutoscale() has to be done in PostInit()
  }
  else
  {
    this->scaling.value = manual_scale;
    this->scaling.value *= base_->optimization->GetObjective()->task == Optimization::MAXIMIZE ? -1.0 : 1.0;
  }
  LOG_TRACE(optimizer) << "Scale::Scale() target=" << target << " tol=" << tol << " manual_scale=" 
                       << manual_scale << " -> PostInit() pending!";
}

void BaseOptimizer::Scale::PostInit()
{
  // calculate the autoscale values (more or less expensive, never cheep)
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
  StdVector<double> grad(base_->optimization->GetDesign()->data.GetSize());
  // make a temporary design as a copy from design space to copy it back there :)
  StdVector<double> data(grad.GetSize());
  // copy the design to our temporary space
  int design_id = base_->optimization->GetDesign()->WriteDesignToExtern(data.GetPointer());

  // evaluate the the gradient -> will be cheep in restart case
  bool good = base_->EvalGradObjective(grad.GetSize(), data.GetPointer(), grad.GetPointer());
  assert(good); // needs to be good as tol = set to 0.0;
  assert(opt_scaling.value != 0.0);
  
  // our new scaling is the optimal scaling for now!
  scaling.design_id = opt_scaling.design_id;
  scaling.value = opt_scaling.value;
  
  // reset the tolerance
  tol = tol_save;
  
  LOG_TRACE(optimizer) << "Scale::CalcAutoscale(): scale=" << scaling.value << " desing=" 
                       << scaling.design_id << " needed_eval=" << (opt_scaling.design_id != design_id);
}


/** Checks the current gradients in the design space against the tolerance.
 * @return true if not active or no tolerance or within tolerance */
bool BaseOptimizer::Scale::CheckScaling(int n, double* grad) 
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
    opt_scaling.value *= base_->optimization->GetObjective()->task == Optimization::MAXIMIZE ? -1.0 : 1.0;
    
    LOG_TRACE(optimizer) << "Scale::check::set opt scaling: design=" << opt_scaling.design_id << " target=" << target 
                         << " opt_scaling=" << opt_scaling.value << " maximize="
                         << (base_->optimization->GetObjective()->task == Optimization::MAXIMIZE);
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
  if(scaling.value == 0.0)
  {
    return "No scaling";
  }

  std::ostringstream os;
  
  if(target != 0.0)
  {
    os << "Autosale target=" << target;
    if(tol != 0.0) os << " tolerance=" << tol;
  }
  else
  {
    os << "Scale";
  }
  os << " scaling=" << scaling.value;
  
  return os.str();
}

BaseOptimizer::BaseOptimizer(Optimization* optimization, ParamNode* pn)
{ 
  this->optimization = optimization;
  this->objective = NULL;
  this->restart_requested = false;
  this->design_.design_id = -1;
  this->design_.value = 0.0;
  this->optimizer_pn_ = pn;
}

BaseOptimizer::~BaseOptimizer()
{ 
  if(objective != NULL) { delete objective; objective = NULL; }
}

void BaseOptimizer::PostInit(double manual_scaling)
{
  ParamNode* as = optimizer_pn_ != NULL ? optimizer_pn_->Get("autoscale", false) : NULL;
  objective = new Scale(this, as, manual_scaling);
  objective->PostInit();
}

std::string BaseOptimizer::LogFileHeader()
{
  return objective->target != 0.0 ? "\tmax_f_grad\tscale\topt_scale" : "\tmax_f_grad";
}


void BaseOptimizer::LogFileLine(std::ofstream* out)
{
  *out << "\t" << objective->current.value;

  if(objective->target != 0.0)
    *out << "\t" << objective->scaling.value 
         << "\t" << objective->opt_scaling.value;
}

double BaseOptimizer::EvalObjective(int n, const double* x)
{
  assert(optimization->GetDesign()->data.GetSize() == (unsigned int) n);

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
    assert(design_.value == optimization->CalcObjective());
    need_eval = false;
  }

  LOG_DBG(optimizer) << "EvalObjective: x_avg=" << Average(x, n) 
                     << " std_dev=" << StandardDeviation(x, n) 
                     << " is_new=" << need_eval << " -> " << design_.value;

  return design_.value;
}

bool BaseOptimizer::EvalGradObjective(int n, const double* x, double* grad_f)
{
  assert(optimization->GetDesign()->data.GetSize() == (unsigned int) n);

  // set the design and see if it is a new one
  int new_design = optimization->GetDesign()->ReadDesignFromExtern(x);
  
  LOG_DBG2(optimizer) << "EvalGradObjective() extern design_id=" << new_design << " design.design_id=" << design_.design_id;    
  
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
    assert(design_.value == optimization->CalcObjective());
  }
  
  LOG_DBG2(optimizer) << "EvalGradObjective: call CalcObjectiveGradient()";
  // calc our gradient - it is not stored anywhere
  optimization->CalcObjectiveGradient(grad_f);

  // check the scaling - harmless if not restarted autoscale
  restart_requested = false;
  LOG_DBG2(optimizer) << "EvalGradObjective: call CheckScaling()";
  if(!objective->CheckScaling(n, grad_f))
  {
    LOG_DBG(optimizer) << "EvalGradObjective: scaling violated";
    assert(design_.design_id == objective->opt_scaling.design_id);
    restart_requested = true;
  }
  
  LOG_DBG(optimizer) << "EvalGradObjective: design=" << design_.design_id << " need_eval=" << need_eval 
                     << " -> grad.scale=" << objective->scaling.value << " grad.opt_scaling="
                     << objective->opt_scaling.value;
  
  return !restart_requested;
}

