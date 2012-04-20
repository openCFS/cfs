#include "Optimization/Optimizer/FeasPP.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "MatVec/scrs_matrix.hh"
#include "Utils/tools.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"

#include <limits>


DECLARE_LOG(feasPP)
DEFINE_LOG(feasPP, "feas_pp")

using namespace CoupledField;
using std::pow;
using std::string;

Enum<FeasPP::Globalization> FeasPP::global;

FeasPP::FeasPP(Optimization* opt, PtrParamNode pn) : BaseOptimizer(opt, pn, Optimization::FEAS_PP_SOLVER)
{
  obj = NULL;
  n = 0;
  m = 0;
  hessian = NULL;

  global.SetName("FeasPP::Globalization");
  global.Add(NONE, "none");
  global.Add(BACKTRACKING, "backtracking");

  approx_linear = this_opt_pn_ != NULL ? this_opt_pn_->Get("approx_linear")->As<bool>() : false;
  approx_feasibility = this_opt_pn_ != NULL ? this_opt_pn_->Get("approx_feasible")->As<bool>() : false;
  non_approx_constraints = true; // set it PostInit()
  global_ = this_opt_pn_ != NULL ? global.Parse(this_opt_pn_->Get("globalize")->As<string>()) : NONE;
  convex_tau = this_opt_pn_ != NULL ? this_opt_pn_->Get("convex_tau")->As<double>() : 0.0; // switch off

  PostInitScale(1.0);

  PtrParamNode ipopt_pn;
  if(this_opt_pn_ != NULL)
    ipopt_pn = this_opt_pn_->Get("ipopt", ParamNode::PASS);
  ipopt = new FeasSubProblem(this, ipopt_pn);

}

FeasPP::~FeasPP()
{
  if(obj != NULL) { delete obj; obj = NULL; }

  for(unsigned int i = 0; i < constr.GetSize(); i++)
    if(constr[i] != NULL)
      { delete constr[i]; constr[i] = NULL; }

  if(hessian != NULL) { delete hessian; hessian = NULL; }
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
  obj = new Approximation(this, -1, approx_linear || !f->IsLinear());

  constr.Resize(m);
  non_approx_constraints = false;
  for(unsigned int i = 0; i < m; i++)
  {
    Condition* g = optimization->constraints.view->Get(i);

    bool approx = approx_linear || !g->IsLinear();
    if(g->IsFeasibilityConstraint()) // feasibility is more important than linear
      approx = approx_feasibility;
    if(!approx)
      non_approx_constraints = true;

    Approximation* func = new Approximation(this, i, approx);
    constr[i] = func;

    LOG_DBG3(feasPP) << "FP:PI i=" << i << " g=" << func->ToString() << " approx=" << func->approximate << " grad=" << func->outer_grad.GetSize() << " pattern=" << func->jac_pattern.ToString();
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

  SetupHessian();
}

void FeasPP::SolveProblem()
{
  // start with iteration 0 which is the initial design
  int max_iter = optimization->GetMaxIterations();

  // evaluate initial iteration
  UpdateToCurrentStep();
  optimization->CommitIteration();

  bool converged = false;
  while(!optimization->DoStopOptimization() && optimization->GetCurrentIteration() <= max_iter && !converged)
  {
    int iter = optimization->GetCurrentIteration();
    PtrParamNode in;
    try
    {
      // solve sub-problem
      // ipopt->Init();
      in = info_->Get(ParamNode::PROCESS)->Get("feasPP")->Get("subsolver", ParamNode::APPEND);
      in->Get("iteration")->SetValue(iter);
      ipopt->SolveProblem(in);
    }
    catch(Exception& ex)
    {
      std::string msg = "subproblem failed in major iteration " + lexical_cast<string>(iter) + ": " + ex.GetMsg();

      PtrParamNode summary = optimization->optInfoNode->Get(ParamNode::SUMMARY);
      summary->Get("break/converged")->SetValue("no");
      summary->Get("problem")->SetValue("FeasPP/IPOPT: " + msg);

      throw Exception(msg);
    }

    assert(ipopt->x_final.GetSize() == n);
    optimization->GetDesign()->ReadDesignFromExtern(ipopt->x_final.GetPointer());

    LOG_DBG2(feasPP) << "SP: it=" << iter << " c_grad = [" << obj->outer_grad.ToString() << "]";
    LOG_DBG2(feasPP) << "SP: it=" << iter << " x_old = [" << x_outer.ToString() << "]";
    LOG_DBG2(feasPP) << "SP: it=" << iter << " x_new_org = [" << ipopt->x_final.ToString() << "]";

    Vector<double> x_old(n);
    x_old.Fill(x_outer.GetPointer(), n);
    Vector<double> x_new(n);
    x_new.Fill(ipopt->x_final.GetPointer(), n);
    Vector<double> c_grad(n);
    c_grad.Fill(obj->outer_grad.GetPointer(), n);

    Vector<double> d(n);
    d = x_new - x_old; // principal direction

    // needs to be smaller zero
    in->Get("angle")->SetValue(d * c_grad);


    if(global_ == BACKTRACKING)
    {
      BTI bti = Backtracking(x_old, d, ipopt->x_final);
      in->Get("steps")->SetValue(bti.steps);
      in->Get("org_dx_norm")->SetValue(bti.org_dx);
      in->Get("curr_dx_norm")->SetValue(bti.curr_dx);
      in->Get("msg")->SetValue("backtracking returned to old point");
      if(bti.old_point_is_optimal)
        converged = true;
    }

    PtrParamNode summary = optimization->optInfoNode->Get(ParamNode::SUMMARY);
    summary->Get("break/converged")->SetValue(converged ? "yes" : "no");

    optimization->GetDesign()->WriteDesignToExtern(x_outer); // only for the LOG below, can be removed as it is done in UpdateToCurrentStep()
    LOG_DBG2(feasPP) << "SP: it=" << iter << " x_new_curr = [" << x_outer.ToString() << "]";

/*    for(unsigned int i = 0; i < optimization->GetDesign()->elements; i++)
    {
      Matrix<double> E;
      optimization->GetDesign()->GetErsatzMaterialTensor(E, PLANE_STRAIN, optimization->GetDesign()->data[i].elem, DesignElement::NO_DERIVATIVE, DesignMaterial::HILL_MANDEL);
      LOG_DBG2(feasPP) << "SP sps i=" << i << " -> " << E.ToString(2);
    }*/

    // evaluate all functions such that we have the function values for the current design
    // the subproblem is based on the approximated values only
    UpdateToCurrentStep();
    optimization->CommitIteration();
  }
}


FeasPP::BTI FeasPP::Backtracking(const Vector<double>& x_old, const Vector<double>& d, const StdVector<double>& std_x_new)
{
  assert(x_old.GetSize() == n);
  assert(std_x_new.GetSize() == n);

  BTI result;

  LOG_DBG2(feasPP) << "FP:B dx_org=" << d.ToString(0, ' ');

  result.org_dx = d.NormL2();

  Vector<double> x_new(n);

  /*for(double t = 0.0; t <= 1.0001; t += 0.002)
  {
    x_new = x_old + t * d;
    optimization->GetDesign()->ReadDesignFromExtern(x_new.GetPointer());
    double ov = EvalObjective(n, x_new.GetPointer(), true);
    // LOG_DBG2(feasPP) << "FP:B t/ov: " << t << "\t" << ov;
    std::cout << "FP:B t/ov: " << t << "\t" << ov << std::endl;
  }*/

  int eval = 0;
  // obj->outer_value is proper value for an approximation of the objective function
  // in case of non->approximation ipopt made the globalization!
  for(double t = 2.0, ov = std::numeric_limits<double>::max(); ov >= obj->outer_val && eval < 25; eval++) // enter at least once!
  {
    t *= 0.5;
    x_new = t * d; // temporary only!
    result.curr_dx = x_new.NormL2();
    LOG_DBG2(feasPP) << "FP:B dx_new=" << x_new.ToString(0, ' ');
    x_new += x_old;
    LOG_DBG2(feasPP) << "FP:B x_new=" << x_new.ToString(0, ' ');
    optimization->GetDesign()->ReadDesignFromExtern(x_new.GetPointer());
    ov = EvalObjective(n, x_new.GetPointer(), true);
    LOG_DBG2(feasPP) << "FP:B e=" << eval << " t=" << t << " ov=" << ov << " old_ov=" << obj->outer_val << " x=" << x_new.ToString(0, ' ');
  }
  assert(!(!obj->approximate && eval > 1)); // see comment above

  result.old_point_is_optimal = eval >= 25;

  result.steps = eval;
  return result;
}

void FeasPP::UpdateToCurrentStep()
{
  DesignSpace* space = optimization->GetDesign();
  space->WriteDesignToExtern(x_outer);

  LOG_DBG(feasPP) << "UTCP x=" << x_outer.ToString();

  // prepare all functions for the present design
  obj->outer_val = EvalObjective(n, x_outer.GetPointer(), true);
  EvalGradObjective(n, x_outer.GetPointer(), true, obj->outer_grad);
  LOG_DBG3(feasPP) << "FP:UTCP obj=" << obj->ToString() << " outer_val=" << obj->outer_val << " grad=" << obj->outer_grad.ToString();

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

    LOG_DBG3(feasPP) << "FP:UTCP g[" << i << "]=" << constr[i]->ToString() << " outer_val=" << constr[i]->outer_val << " grad=" << constr[i]->outer_grad.ToString();
  }
  optimization->constraints.view->Done(); // reset slope constraint to global mode
}

void FeasPP::SetupHessian()
{
  hessian = new compressed_matrix<double>(n, n);

  obj->AddHessianPattern(*hessian);

  for(unsigned int i = 0; i < m; i++)
    constr[i]->AddHessianPattern(*hessian);

  LOG_DBG3(feasPP) << "FP:SH nnz= " << hessian->nnz() << " -> " << hessian;
}

Approximation::Approximation(FeasPP* feas_pp, int constraint_idx, bool approx)
{
  this->common = feas_pp;
  this->constraint_idx = constraint_idx;
  this->outer_val = -1.0;
  this->approximate = approx;

  Function*  f = GetFunction();
  Condition* g = f->IsObjective() ? NULL : GetCondition();

  jac_pattern = f->GetSparsityPattern();
  outer_grad.Resize(jac_pattern.GetSize());

  if(approximate)
  {
    // the Hessian pattern of the approximations is diagonal
    // the Jacobian might be sparse and the Hessian has the same size as we multiply with the Jacobian
    hess_pattern.Resize(jac_pattern.GetSize(), 2);
    for(unsigned int i = 0; i < jac_pattern.GetSize(); i++)
    {
      hess_pattern(i, 0) = jac_pattern[i];
      hess_pattern(i, 1) = jac_pattern[i];
    }
  }
  else
    hess_pattern = f->GetHessianSparsityPattern();

  lower = f->IsObjective() ? -1.0 * common->GetInfBound() : g->GetBoundValue();
  upper = f->IsObjective() ?  1.0 * common->GetInfBound() : g->GetBoundValue();

  if(g != NULL)
  {
    if(g->GetBound() == Condition::LOWER_BOUND) upper =  1.0 * common->GetInfBound();
    if(g->GetBound() == Condition::UPPER_BOUND) lower = -1.0 * common->GetInfBound();
  }
}

double Approximation::Evaluate(const double* x_inner, Eval eval, StdVector<double>* out)
{
  assert(outer_grad.GetSize() == jac_pattern.GetSize());
  assert((eval == FUNC && out == NULL) || (eval != FUNC && out != NULL));
  assert(eval != GRAD || out->GetSize() == jac_pattern.GetSize());
  assert(eval != HESSIAN || out->GetSize() == hess_pattern.GetNumRows());

  double result = 0.0;

  if(approximate)
    result = EvalApproximation(x_inner, eval, out);
  else
    result = EvalDirect(x_inner, eval, out);

  LOG_DBG2(feasPP) << "A:E f=" << ToString() << " d=" << eval << " -> " << (eval != FUNC ? "..." : boost::lexical_cast<std::string>(result)) << " dx=" << NormL2(common->x_outer.GetPointer(), x_inner, common->n);
  LOG_DBG3(feasPP) << "A:E f=" << ToString() << " d=" << eval << " -> " << (eval != FUNC ? out->ToString() : boost::lexical_cast<std::string>(result));
  return result;
}

double Approximation::EvalApproximation(const double* x_inner, Eval eval, StdVector<double>* out)
{
  // Svanberg (2):
  // f(xi) = f(xo) - sum(p/(U-xo) + q/(xo-L)) + sum(p/(U-xi) + q/(xi-L))
  // d f(xi)/d xi = + p/((U-xi)**2) - q/((xi-L)**2)
  double result = eval == FUNC ? outer_val : 0.0;

  // the function evaluation is based on the (pattern of the) Jacobian
  for(unsigned int e = 0, en = jac_pattern.GetSize(); e < en; e++)
  {
    unsigned int j = jac_pattern[e]; // Hessian is diagonal with row=col
    assert((eval != HESSIAN)  || ((hess_pattern(e,0) == hess_pattern(e,1)) && (hess_pattern(e,0) == jac_pattern[e])));
    double grad = outer_grad[e];
    double xo   = common->x_outer[j];
    double xi   = x_inner[j];
    double U    = common->U[j];
    double L    = common->L[j];

    assert(xo < U && xi < U && xo > L && xi > L);

    double p = grad > 0 ? (U-xo)*(U-xo) *  grad : 0.0;
    double q = grad < 0 ? (xo-L)*(xo-L) * -grad : 0.0;

    LOG_DBG3(feasPP) << "A:E f=" << ToString() << " j=" << j << " xo=" << xo << " xi=" << xi << " grad=" << grad << " U=" << U << " L=" << L << " p=" << p << " q=" << q << " ov=" << outer_val;

    // (6) in Sonjas paper, only for the objective, this is only tau (or 0.0 if not applicable)
    double tau = constraint_idx == -1 ? common->convex_tau : 0.0;

    // maxima:
    // fp: tau*(xi-xo)**2/(U-xi);
    // fn: tau*(xi-xo)**2/(xi-L);


    switch(eval)
    {
    case FUNC:
    {
      double r_ = p/(U-xo) + q/(xo-L);
      double s  = p/(U-xi) + q/(xi-L);
      // (6) in Sonjas paper, only for the objective
      double c = tau * (grad >= 0 ? (xi-xo)*(xi-xo)/(U-xi) : (xi-xo)*(xi-xo)/(xi-L));
      result += -r_ + s + c;
      LOG_DBG3(feasPP) << "A:E f=" << ToString() << " func r_=" << r_ << " s=" << s << " c=" << c << " -> " << result;
      break;
    }
    case GRAD:
    {
      // diff(fp,xi); diff(fn,xi)
      double c = grad >= 0 ? 2.0*tau*(xi-xo)/(U-xi) + tau*(xi-xo)*(xi-xo)/((U-xi)*(U-xi))
                           : 2.0*tau*(xi-xo)/(xi-L) - tau*(xi-xo)*(xi-xo)/((xi-L)*(xi-L));
      (*out)[e] = + p/((U-xi)*(U-xi)) - q/((xi-L)*(xi-L)) + c;
      LOG_DBG3(feasPP) << "A:E f=" << ToString() << " c=" << c << " grad out[" << e << "]=" << (*out)[e];
      break;
    }
    case HESSIAN:
    {
      // diff(diff(fp,xi),xi); diff(diff(fn,xi),xi);
      double c = grad >= 0 ? 2.0*tau/(U-xi) + 4.0*tau*(xi-xo)/((U-xi)*(U-xi)) + 2.0*tau*(xi-xo)*(xi-xo)/((U-xi)*(U-xi)*(U-xi))
                           : 2.0*tau/(xi-L) - 4.0*tau*(xi-xo)/((xi-L)*(xi-L)) + 2.0*tau*(xi-xo)*(xi-xo)/((xi-L)*(xi-L)*(xi-L));
      (*out)[e] = 2.0*p/((U-xi)*(U-xi)*(U-xi)) + 2.0*q/((xi-L)*(xi-L)*(xi-L)) + c;
      LOG_DBG3(feasPP) << "A:E f=" << ToString() << " c=" << c << " hess out[" << e << "]=" << (*out)[e];
      break;
    }
    }
  }
  return result;
}

double Approximation::EvalDirect(const double* x_inner, Eval eval, StdVector<double>* out)
{
  double result = 0.0;

  if(GetFunction()->IsObjective())
  {
    switch(eval)
    {
    case FUNC:
      result = common->EvalObjective(common->n, x_inner, true);
      break;
    case GRAD:
      common->EvalGradObjective(common->n, x_inner, true, *out);
      break;
    case HESSIAN:
      GetFunction()->CalcHessian(*out);
      break;
    }
  }
  else
  {
    switch(eval)
    {
    case FUNC:
      result = common->EvalConstraint(GetCondition(), true);
      break;
    case GRAD:
      common->EvalGradConstraint(GetCondition(), 0, true, *out);
      break;
    case HESSIAN:
      LOG_DBG3(feasPP) << "A:ED g=" << ToString() << " hess: pattern= " << hess_pattern.ToString();
      LOG_DBG3(feasPP) << "A:ED g=" << ToString() << " hess: GHSP   = " << GetCondition()->GetHessianSparsityPattern().ToString();
      GetCondition()->CalcHessian(*out);
      break;
    }
  }

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

void Approximation::AddHessianPattern(compressed_matrix<double>& hessian)
{
  if(hess_pattern.GetNumRows() == 0)
    return;

  assert(hess_pattern.GetNumCols() == 2);

  for(unsigned int i = 0; i < hess_pattern.GetNumRows(); i++)
    hessian(hess_pattern(i, 0), hess_pattern(i, 1)) = 1.0;

  LOG_DBG3(feasPP) << "A:AHP f=" << ToString() << " -> " << hess_pattern.ToString(0, false);
}
