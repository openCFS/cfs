#include "Optimization/Optimizer/FeasPP.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"



DECLARE_LOG(feasPP)
DEFINE_LOG(feasPP, "feas_pp")

using namespace CoupledField;
using std::pow;


FeasPP::FeasPP(Optimization* opt, PtrParamNode pn) : BaseOptimizer(opt, pn, Optimization::FEAS_PP_SOLVER)
{
  obj = NULL;
  n = 0;
  m = 0;

  approx_linear = pn->Has("feasPP/approx_linear") ? pn->Get("feasPP/approx_linear")->As<bool>() : false;
  approx_feasibility = pn->Has("feasPP/approx_feasible") ? pn->Get("feasPP/approx_feasible")->As<bool>() : false;
  non_approx_constraints = true; // set it PostInit()

  PostInitScale(1.0);

  ipopt = new FeasSubProblem(this, pn->Get("feasPP/ipopt", ParamNode::PASS));

}

FeasPP::~FeasPP()
{
  if(obj != NULL) { delete obj; obj = NULL; }

  for(unsigned int i = 0; i < constr.GetSize(); i++)
    if(constr[i] != NULL)
      { delete constr[i]; constr[i] = NULL; }

}

void FeasPP::PostInit()
{
  assert(obj == NULL); // call once
  assert(optimization->objectives.data.GetSize() == 1); // trivial case only

  n = optimization->GetDesign()->GetNumberOfVariables();
  m = optimization->constraints.view->GetNumberOfActiveConstraints();
  x_outer.Resize(n, 0.0);

  // setup functions
  Function* f = optimization->objectives.data[0];
  obj = new Approximation(this, -1);
  obj->pattern = f->GetSparsityPattern();
  obj->outer_grad.Resize(obj->pattern.GetSize());
  obj->approximate = approx_linear || !f->IsLinear();

  constr.Resize(m);

  non_approx_constraints = false;

  for(unsigned int i = 0; i < m; i++)
  {
    Condition* g = optimization->constraints.view->Get(i);
    Approximation* func = new Approximation(this, i);
    constr[i] = func;

    func->pattern = g->GetSparsityPattern();
    func->outer_grad.Resize(func->pattern.GetSize());

    func->lower = g->GetBoundValue();
    func->upper = g->GetBoundValue();
    if(g->GetBound() == Condition::LOWER_BOUND) func->upper =  GetInfBound();
    if(g->GetBound() == Condition::UPPER_BOUND) func->lower = -GetInfBound();

    func->approximate = approx_linear || !g->IsLinear();
    if(g->IsFeasibilityConstraint()) // feasibility is more important than linear
      func->approximate = approx_feasibility;
    if(!func->approximate)
      non_approx_constraints = true;

    LOG_DBG3(feasPP) << "FP:PI i=" << i << " g=" << func->ToString() << " approx=" << func->approximate << " grad=" << func->outer_grad.GetSize() << " pattern=" << func->pattern.ToString();
  }
  optimization->constraints.view->Done();

  // setup fixed asymptotes
  L.Resize(n);
  U.Resize(n);

  assert(n == optimization->GetDesign()->data.GetSize()); // TODO include aux design!
  for(unsigned int i = 0; i < n; i++)
  {
    DesignElement& de = optimization->GetDesign()->data[i];
    L[i] = de.GetLowerBound() > 0 ? 0.0 : de.GetLowerBound() * 1.1; // be smaller!!
    U[i] = de.GetUpperBound() * 1.1;
  }
  optimization->constraints.view->Done();

}

void FeasPP::SolveProblem()
{
  // start with iteration 0 which is the initial design
  int iter = 0;
  int max_iter = optimization->GetMaxIterations();

  while(!optimization->DoStopOptimization() && iter++ <= max_iter)
  {
    try
    {
      SolveSubProblem();
    }
    catch(Exception& ec)
    {

    }

    assert(ipopt->x_final.GetSize() == n);
    optimization->GetDesign()->ReadDesignFromExtern(ipopt->x_final.GetPointer());
    optimization->CommitIteration();
  }
}

void FeasPP::SolveSubProblem()
{
  DesignSpace* space = optimization->GetDesign();
  space->WriteDesignToExtern(x_outer);

  // prepare all functions for the present design

  obj->outer_val = EvalObjective(n, x_outer.GetPointer(), true);
  EvalGradObjective(n, x_outer.GetPointer(), true, obj->outer_grad);
  LOG_DBG3(feasPP) << "FP:SSP obj=" << obj->ToString() << " outer_val=" << obj->outer_val << " grad=" << obj->outer_grad.ToString();

  StdVector<double> g_val(m);
  EvalConstraints(n, x_outer.GetPointer(), m, true, g_val.GetPointer());
  // this is copy and paste from EvalGradConstraints()
  // reset values of the constraint gradients before the loop
  // as it also contains a loop over all the design elements
  optimization->GetDesign()->Reset(DesignElement::CONSTRAINT_GRADIENT, DesignElement::DEFAULT);
  for(unsigned int i = 0; i < m; i++)
  {
    constr[i]->outer_val = g_val[i];
    EvalGradConstraint(optimization->constraints.view->Get(i), 0, true, constr[i]->outer_grad);

    LOG_DBG3(feasPP) << "FP:SSP g[" << i << "]=" << constr[i]->ToString() << " outer_val=" << constr[i]->outer_val << " grad=" << constr[i]->outer_grad.ToString();
  }

  optimization->constraints.view->Done(); // reset slope constraint to global mode

  // solve sub-problem
  ipopt->Init();
  ipopt->SolveProblem(this->info_->Get("subsolver", ParamNode::APPEND)); // KILLME
}



Approximation::Approximation(FeasPP* feas_pp, int constraint_idx)
{
  this->common = feas_pp;
  this->constraint_idx = constraint_idx;
  this->outer_val = -1.0;
}

double Approximation::Eval(const double* x_inner, bool gradient, StdVector<double>* out)
{
  assert(outer_grad.GetSize() == pattern.GetSize());
  assert((gradient && out != NULL) || (!gradient && out == NULL));
  assert((out == NULL) || (out->GetSize() == pattern.GetSize()));

  double result = 0.0;

  if(approximate)
  {
    // Svanberg (2):
    // f(xi) = f(xo) - sum(p/(U-xo) + q/(xo-L)) + sum(p/(U-xi) + q/(xi-L))
    // d f(xi)/d xi = + p/((U-xi)**2) - q/((xi-L)**2)
    result = gradient ? 0.0 : outer_val;

    for(unsigned int e = 0; e < pattern.GetSize(); e++)
    {
      unsigned int j = pattern[e];
      double grad = outer_grad[e];
      double xo   = common->x_outer[j];
      double xi   = x_inner[j];
      double U    = common->U[j];
      double L    = common->L[j];

      assert(xo < U && xi < U && xo > L && xi > L);

      double p = grad > 0 ? (U-xo)*(U-xo) *  grad : 0.0;
      double q = grad < 0 ? (xo-L)*(xo-L) * -grad : 0.0;

      LOG_DBG3(feasPP) << "A:E f=" << ToString() << " j=" << j << " xo=" << xo << " xi=" << xi << " grad=" << grad << " U=" << U << " L=" << L << " p=" << p << " q=" << q << " ov=" << outer_val;

      if(!gradient)
      {
        double r_ = p/(U-xo) + q/(xo-L);
        double s  = p/(U-xi) + q/(xi-L);
        result += -r_ + s;

        LOG_DBG3(feasPP) << "A:E f=" << ToString() << " func r_=" << r_ << " s=" << s << " -> " << result;
      }
      else
      {
        (*out)[e] = + p/((U-xi)*(U-xi)) - q/((xi-L)*(xi-L));
        LOG_DBG3(feasPP) << "A:E f=" << ToString() << " grad out[" << e << "]=" << (*out)[e];
      }
    }
  }
  else // no approximation
  {
    if(GetFunction()->IsObjective())
    {
      if(!gradient)
        result = common->EvalObjective(common->n, x_inner, true);
      else
        common->EvalGradObjective(common->n, x_inner, true, *out);
    }
    else
    {
      if(!gradient)
        result = common->EvalConstraint(GetCondition(), true);
      else
        common->EvalGradConstraint(GetCondition(), 0, true, *out);
    }
  }

  LOG_DBG3(feasPP) << "A:E f=" << ToString() << " d=" << gradient << " x=" << StdVector<double>::ToString(common->n, x_inner) << " -> " << (gradient ? out->ToString() : boost::lexical_cast<std::string>(result));
  return result;
}

std::string Approximation::ToString()
{
  return Function::type.ToString(GetFunction()->GetType());
}

Condition* Approximation::GetCondition()
{
  Function* f = GetFunction();
  assert(!f->IsObjective());
  return static_cast<Condition*>(f);
}

Function* Approximation::GetFunction()
{
  assert(constraint_idx >= -1);
  assert(constraint_idx > 0 || common->optimization->objectives.data.GetSize() == 1);

  if(constraint_idx == -1)
    return common->optimization->objectives.data[0];
  else
    return common->optimization->constraints.view->Get(constraint_idx);
}
