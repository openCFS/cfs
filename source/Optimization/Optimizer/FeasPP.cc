#include "Optimization/Optimizer/FeasPP.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "MatVec/scrs_matrix.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"




DECLARE_LOG(feasPP)
DEFINE_LOG(feasPP, "feas_pp")

using namespace CoupledField;
using std::pow;
using std::string;

FeasPP::FeasPP(Optimization* opt, PtrParamNode pn) : BaseOptimizer(opt, pn, Optimization::FEAS_PP_SOLVER)
{
  obj = NULL;
  n = 0;
  m = 0;
  hessian = NULL;

  approx_linear = pn->Has("feasPP/approx_linear") ? pn->Get("feasPP/approx_linear")->As<bool>() : false;
  approx_feasibility = pn->Has("feasPP/approx_feasible") ? pn->Get("feasPP/approx_feasible")->As<bool>() : false;
  non_approx_constraints = true; // set it PostInit()

  PostInitScale(1.0);

  ipopt = new FeasSubProblem(this, pn->Get("feasPP/ipopt", ParamNode::PASS));


  compressed_matrix<double> m(3,3);
  for (unsigned i = 0; i < m.size1 (); ++ i)
    for (unsigned j = 0; j < m.size2 (); ++ j)
      if(j >= i)
        m (i, j) = 3.0 * i + j;

  m.clear();
  //m.resize(3,3);
  assert(m.nnz() == 0);
  m(0,0) = 0.0;
  assert(m.nnz() == 1);

  std::cout << m << std::endl;

  unbounded_array<double> vd = m.value_data();
  for(unbounded_array<double>::iterator it = vd.begin(); it != vd.end(); ++it)
  {
  //  std::cout << "vd: " << *it->index() << ": " << *it << std::endl;
  }
  for(compressed_matrix<double>::const_iterator1 it = m.begin1(); it != m.end1(); ++it)
  {
    //std::cout << "it1: " << *it << " idx1=" << it.index1() << " idx2=" << it.index2() << std::endl;
    for(compressed_matrix<double>::const_iterator2 it2 = it.begin(); it2 != it.end(); ++it2)
    {
      // std::cout << "it2: " << *it2 << " idx1=" << it2.index1() << " idx2=" << it2.index2() << std::endl;
      //std::cout << "it2: " << *it2 << std::endl;
    }
  }
  for(compressed_matrix<double>::const_iterator2 it = m.begin2(); it != m.end2(); ++it)
  {
  //  std::cout << "it2: " << *it << std::endl;
  }

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

  while(!optimization->DoStopOptimization() && optimization->GetCurrentIteration() <= max_iter)
  {
    int iter = optimization->GetCurrentIteration();
    try
    {
      // solve sub-problem
      ipopt->Init();
      PtrParamNode in = info_->Get(ParamNode::PROCESS)->Get("feasPP")->Get("subsolver", ParamNode::APPEND);
      in->Get("iteration")->SetValue(iter);
      ipopt->SolveProblem(in);
    }
    catch(Exception& ex)
    {
      string msg = "Ignore failed subproblem in major iteration " + boost::lexical_cast<string>(iter) + " :" + ex.GetMsg();
      info_->Get(ParamNode::PROCESS)->Get("feasPP")->Get(ParamNode::WARNING)->SetValue(msg);
    }

    assert(ipopt->x_final.GetSize() == n);
    optimization->GetDesign()->ReadDesignFromExtern(ipopt->x_final.GetPointer());

    // evaluate all functions such that we have the function values for the current design
    // the subproblem is based on the approximated values only
    UpdateToCurrentStep();
    optimization->CommitIteration();
  }
}


void FeasPP::UpdateToCurrentStep()
{
  DesignSpace* space = optimization->GetDesign();
  space->WriteDesignToExtern(x_outer);

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

  LOG_DBG3(feasPP) << "A:E f=" << ToString() << " d=" << eval << " x=" << StdVector<double>::ToString(common->n, x_inner) << " -> " << (eval != FUNC ? out->ToString() : boost::lexical_cast<std::string>(result));
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

    switch(eval)
    {
    case FUNC:
    {
      double r_ = p/(U-xo) + q/(xo-L);
      double s  = p/(U-xi) + q/(xi-L);
      result += -r_ + s;
      LOG_DBG3(feasPP) << "A:E f=" << ToString() << " func r_=" << r_ << " s=" << s << " -> " << result;
      break;
    }
    case GRAD:
      (*out)[e] = + p/((U-xi)*(U-xi)) - q/((xi-L)*(xi-L));
      LOG_DBG3(feasPP) << "A:E f=" << ToString() << " grad out[" << e << "]=" << (*out)[e];
      break;
    case HESSIAN:
      (*out)[e] = 2.0*p/((U-xi)*(U-xi)*(U-xi)) + 2.0*q/((xi-L)*(xi-L)*(xi-L));
      LOG_DBG3(feasPP) << "A:E f=" << ToString() << " hess out[" << e << "]=" << (*out)[e];
      break;
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
