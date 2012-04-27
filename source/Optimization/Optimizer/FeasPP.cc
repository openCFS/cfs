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
using std::max;
using std::abs;
using std::string;

Enum<FeasPP::Globalization> FeasPP::global;

FeasPP::FeasPP(Optimization* opt, PtrParamNode pn) : BaseOptimizer(opt, pn, Optimization::FEAS_PP_SOLVER)
{
  obj = NULL;
  n = 0;
  m = 0;
  m_c = 0;
  m_e = 0;
  hessian = NULL;

  global.SetName("FeasPP::Globalization");
  global.Add(NONE, "none");
  global.Add(BACKTRACKING, "backtracking");
  global.Add(AUG_LAGRANGIAN, "augmentedLagrangian");

  approx_linear = this_opt_pn_ != NULL ? this_opt_pn_->Get("approx_linear")->As<bool>() : false;
  approx_feasibility = this_opt_pn_ != NULL ? this_opt_pn_->Get("approx_feasible")->As<bool>() : false;
  non_approx_constraints = true; // set it PostInit()
  global_ = this_opt_pn_ != NULL ? global.Parse(this_opt_pn_->Get("globalize")->As<string>()) : NONE;
  convex_tau = this_opt_pn_ != NULL ? this_opt_pn_->Get("convex_tau")->As<double>() : 0.0; // switch off
  rho_init_  = this_opt_pn_ != NULL && this_opt_pn_->Has("augmentedLagranian") ? this_opt_pn_->Get("augmentedLagranian/rho")->As<double>() : 1.0;
  decrease_  = this_opt_pn_ != NULL && this_opt_pn_->Has("augmentedLagranian") ? this_opt_pn_->Get("augmentedLagranian/decrease")->As<double>() : 0.1;
  stepwidth_ = this_opt_pn_ != NULL && this_opt_pn_->Has("augmentedLagranian") ? this_opt_pn_->Get("augmentedLagranian/stepwidth")->As<double>() : 0.5;
  // same value for backtracking and augmented lagrangian. Re-read in the latter case
  min_step_  = this_opt_pn_ != NULL && this_opt_pn_->Has("backtracking") ? this_opt_pn_->Get("backtracking/min_step")->As<double>() : 1e-9;
  if(global_ == AUG_LAGRANGIAN)
  {
    min_step_  = this_opt_pn_ != NULL && this_opt_pn_->Has("augmentedLagranian") ? this_opt_pn_->Get("augmentedLagranian/min_step")->As<double>() : 1e-9;

    if(rho_init_ <= 0.0) throw Exception("'augmentedLagrangian/rho' needs to be positive");
    if(decrease_ <= 0.0 || decrease_ >= 1.0) throw Exception("'augmentedLagrangian/decrease' needs to be in (0,1)");
    if(stepwidth_ <= 0.0 || stepwidth_ >= 1.0) throw Exception("'augmentedLagrangian/stepwidth' needs to be in (0,1)");
  }

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

  if(global_ == AUG_LAGRANGIAN)
    rho.Resize(m, rho_init_);

  // setup functions
  Function* f = optimization->objectives.data[0];
  obj = new Approximation(this, -1, approx_linear || !f->IsLinear());

  constr.Resize(m);
  non_approx_constraints = false;
  for(unsigned int i = 0; i < m; i++)
  {
    Condition* g = optimization->constraints.view->Get(i);

    bool approx = approx_linear || !g->IsLinear();
    if(g->IsFeasibilityConstraint()) { // feasibility is more important than linear
      approx = approx_feasibility;
      m_e++;
    }
    else {
      m_c++;
    }
    if(!approx)
      non_approx_constraints = true;
    Approximation* func = new Approximation(this, i, approx);
    constr[i] = func;

    LOG_DBG3(feasPP) << "FP:PI i=" << i << " g=" << func->ToString() << " approx=" << func->approximate << " grad=" << func->outer_grad.GetSize() << " pattern=" << func->jac_pattern.ToString();
  }
  optimization->constraints.view->Done();

  assert(m == m_c + m_e);

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

  PtrParamNode summary = optimization->optInfoNode->Get(ParamNode::SUMMARY);

  bool converged = false;
  LSR ls;
  ls.old_point_is_optimal = false;
  int iter = 0;
  while(!optimization->DoStopOptimization() && iter <= max_iter && !converged)
  {
    // solve sub-problem
    PtrParamNode in = info_->Get(ParamNode::PROCESS)->Get("feasPP")->Get("subsolver", ParamNode::APPEND);
    in->Get("iteration")->SetValue(iter);
    std::string err = ipopt->SolveProblem(in);

    if(err != "")
    {
      std::string msg = "subproblem failed in major iteration " + lexical_cast<string>(iter) + ": " + err;

      summary = optimization->optInfoNode->Get(ParamNode::SUMMARY);
      summary->Get("break/converged")->SetValue("no");
      summary->Get("problem")->SetValue("FeasPP/IPOPT: " + msg);

      throw Exception(msg);
    }

    assert(ipopt->x_final.GetSize() == n);
    optimization->GetDesign()->ReadDesignFromExtern(ipopt->x_final.GetPointer());

    LOG_DBG2(feasPP) << "SP: it=" << iter << " c_grad = [" << obj->outer_grad.ToString() << "]";
    LOG_DBG2(feasPP) << "SP: it=" << iter << " x_old = [" << x_outer.ToString() << "]";
    LOG_DBG2(feasPP) << "SP: it=" << iter << " x_new_org = [" << ipopt->x_final.ToString() << "]";

    Vector<double> f_grad(n);
    f_grad.Fill(obj->outer_grad.GetPointer(), n);

    Vector<double> d(n);
    d = ipopt->x_final - x_outer; // principal direction

    // needs to be smaller zero
    in->Get("angle")->SetValue(d * f_grad);

    if(global_ != NONE)
    {
      if(global_ == BACKTRACKING)
        ls = Backtracking(x_outer, d, ipopt->x_final);
      else
        ls = AugmentedLagrangianLineSearch(iter, ipopt->x_final, x_outer, ipopt->lambda, ipopt->old_lambda);

      in->Get("steps")->SetValue(ls.steps);
      in->Get("step_factor")->SetValue(ls.stepwidth);
      in->Get("org_dx_norm")->SetValue(ls.org_dx);

      if(ls.stepwidth != 1.0)
        in->Get("curr_dx_norm")->SetValue(ls.curr_dx);

      if(ls.old_point_is_optimal)
      {
        in->Get("msg")->SetValue("backtracking returned to old point");
        converged = true;
      }
    }

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

    iter = optimization->GetCurrentIteration();
  }

  summary->Get("break/converged")->SetValue(converged ? "yes" : "no");
  if(ls.old_point_is_optimal)
    summary->Get("reason/msg")->SetValue("linesearch returned to old point");
  if(iter >= max_iter)
    summary->Get("reason/msg")->SetValue("maximum iterations reached");

  if(converged)
    std::cout << "\nFeasPP finished with " << iter << " iterations" << std::endl;
  else
    std::cout << "\nFeasPP did not converge: " << summary->Get("reason/msg")->As<string>() << std::endl;
}


FeasPP::LSR FeasPP::Backtracking(const Vector<double>& x_old, const Vector<double>& d, const Vector<double>& std_x_new)
{
  assert(x_old.GetSize() == n);
  assert(std_x_new.GetSize() == n);


  LOG_DBG2(feasPP) << "FP:B dx_org=" << d.ToString(0, ' ');
  LSR result;
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

  result.steps = 0;
  double t = 2.0;
  // obj->outer_value is proper value for an approximation of the objective function
  // in case of non->approximation ipopt made the globalization!
  for(double ov = std::numeric_limits<double>::max(); ov >= obj->outer_val && t >= min_step_; result.steps++) // enter at least once!
  {
    t *= 0.5;
    x_new = t * d; // temporary only!
    result.curr_dx = x_new.NormL2();
    LOG_DBG2(feasPP) << "FP:B dx_new=" << x_new.ToString(0, ' ');
    x_new += x_old;
    LOG_DBG2(feasPP) << "FP:B x_new=" << x_new.ToString(0, ' ');
    optimization->GetDesign()->ReadDesignFromExtern(x_new.GetPointer());
    ov = EvalObjective(n, x_new.GetPointer(), true);
    LOG_DBG2(feasPP) << "FP:B e=" << result.steps << " t=" << t << " ov=" << ov << " old_ov=" << obj->outer_val << " x=" << x_new.ToString(0, ' ');
  }
  assert(!(!obj->approximate && result.steps > 1)); // see comment above

  result.stepwidth = t;
  result.old_point_is_optimal = (t <= min_step_);
  return result;
}



FeasPP::LSR FeasPP::AugmentedLagrangianLineSearch(int k, const Vector<double>& x, const Vector<double>& z, const Vector<double>& y, const Vector<double>& v)
{
  // the upper part of d
  Vector<double> dx(n);
  dx = z - x;
  // the lower part of d
  Vector<double> dy(m);
  dy = v - y;

  // Algorithm 1 in the paper
  Vector<double> d(dx, dy);

  Vector<double> grad_phi(n + m);
  CalcGradAugmentedLagrangian(x, y, rho, grad_phi);

  double grad_phi_d = grad_phi  * d;
  double phi = CalcAugmentedLagrangian(x, y, rho);

  double eta = CalcEta(convex_tau, x, z);

  // norm(z -x)^2
  double dist = pow(NormL2(z.GetPointer(), x.GetPointer(), n), 2);

  if(k > 0)
  {
    // (22) sufficient decent property
    if(grad_phi_d > -0.5 * eta * dist)
      CalcPenaltyRho(eta, dist, y, v, rho);

    CalcGradAugmentedLagrangian(x, y, rho, grad_phi);
    assert(grad_phi  * d  <= -0.5 * eta * dist);
  }
  double sigma = 1.0;
  Vector<double> x_next(n);
  x_next = x + sigma * dx;
  Vector<double> y_next(m);
  y_next = y + sigma * dy;

  // perform Armijo
  // descent condition (23)
  double sigma_phi = CalcAugmentedLagrangian(x_next, y_next, rho);

  LSR result;
  result.org_dx = dx.NormL2();
  result.steps = 1;

  while((sigma_phi > phi + decrease_ * sigma * grad_phi_d) && sigma >= min_step_)
  {
    LOG_DBG2(feasPP) << "FP:ALLS sigma_phi=" << sigma_phi << " > phi=" << phi << " dec=" << decrease_ << " sig=" << sigma << " sp=" << grad_phi_d;
    sigma *= stepwidth_;
    x_next = x + sigma * dx;
    y_next = y + sigma * dy;

    sigma_phi = CalcAugmentedLagrangian(x_next, y_next, rho);

    result.steps++;
  }

  result.stepwidth = sigma;
  result.curr_dx = sigma  * result.org_dx;
  result.old_point_is_optimal = sigma < min_step_;

  return result;
}

double FeasPP::CalcEta(double tau, const Vector<double>& x_vec, const Vector<double>& z_vec)
{
  assert(x_vec.GetSize() == n && z_vec.GetSize() == n);
  // feasibility paper (18)

  double eta = std::numeric_limits<double>::max();

  const StdVector<double>& grad_f = obj->outer_grad;

  for(unsigned int i = 0; i < n; i++)
  {
    double x = x_vec[i];
    double z = z_vec[i];
    double u = U[i];
    double l = L[i];

    double eta_i;
    if(grad_f[i] >= 0.0)
      eta_i = (grad_f[i] + tau) * (2.0*u - z - x) / ((u-x)*(u-z));
    else
      eta_i = -1.0 * (grad_f[i] - tau) * (-2.0*l + z + x) / ((z-l)*(z-l));

    eta = std::min(eta, eta_i);
  }

  return eta;
}

void FeasPP::CalcPenaltyRho(double eta, double diff, const Vector<double>& y_vec,  const Vector<double>& v_vec, StdVector<double>& rho) const
{
  assert(m == m_e + m_c);

  for(unsigned int i = 0; i < m; i++)
  {
    Condition* g = optimization->constraints.view->Get(i);
    double y = y_vec[i];
    double v = v_vec[i];

    if(!g->IsFeasibilityConstraint())
      rho[i] = max(rho[i], (40.0 * m_c / (eta * diff*diff)) * max(y * abs(v - y), (v-y) * (v-y)));
    else
      rho[i] = max(rho[i], (10.0 * m_e / (eta * diff*diff)) * max(y * abs(v - y), y * y));
  }
  optimization->constraints.view->Done();
}


double FeasPP::CalcAugmentedLagrangian(const Vector<double>& x, const Vector<double>& y, const StdVector<double>& rho)
{
  LOG_DBG3(feasPP) << "FP:CAL x=" << x.ToString(0, ' ');
  LOG_DBG3(feasPP) << "FP:CAL y=" << y.ToString(0, ' ');
  LOG_DBG3(feasPP) << "FP:CAL rho=" << rho.ToString(0, ' ');

  assert(x.GetSize() == n);
  assert(y.GetSize() == m && rho.GetSize() == m);

  // evaluate function values
  optimization->GetDesign()->ReadDesignFromExtern(x.GetPointer());
  double f = EvalObjective(n, x.GetPointer(), true);
  StdVector<double> c(m);
  EvalConstraints(n, x.GetPointer(), m, true, c.GetPointer(), true);

  double result = f;

  for(unsigned int i = 0; i < m; i++)
  {
    Condition* g = optimization->constraints.view->Get(i);
    LOG_DBG3(feasPP) << "FP:CAL i=" << i << " g=" << g->ToString() << " y=" << y[i] << " c=" << c[i] << " rho=" << rho[i]  << " feas=" << g->IsFeasibilityConstraint();

    if(-1.0 * y[i] / rho[i] <= c[i])
      result += y[i] * c[i] + 0.5 * rho[i] * c[i] * c[i];
    else
      result += y[i]*y[i] / (-2.0 * rho[i]);

    LOG_DBG3(feasPP) << "FP:CAL i=" << i << " g=" << g->ToString() << " y=" << y[i] << " c=" << c[i] << " rho=" << rho[i]  << " -> " << result;
  }
  optimization->constraints.view->Done();

  return result;
}

void FeasPP::CalcGradAugmentedLagrangian(const Vector<double>& x, const Vector<double>& y_vec, const StdVector<double>& rho_vec, Vector<double>& grad)
{
  LOG_DBG3(feasPP) << "FP:CGAL x=" << x.ToString(0, ' ');
  LOG_DBG3(feasPP) << "FP:CGAL y=" << y_vec.ToString(0, ' ');
  LOG_DBG3(feasPP) << "FP:CGAL rho=" << rho_vec.ToString(0, ' ');

  assert(x.GetSize() == n);
  assert(y_vec.GetSize() == m && rho_vec.GetSize() == m);
  assert(grad.GetSize() == n + m);

  // we need the constraint values to decide about the case
  StdVector<double> c_vec(m);
  EvalConstraints(n, x.GetPointer(), m, true, c_vec.GetPointer(), true);

  // for the  the x-part of grad, we need function gradients
  // evaluate function gradient values
  optimization->GetDesign()->ReadDesignFromExtern(x.GetPointer());
  StdVector<double> f_grad(n);
  EvalGradObjective(n, x.GetPointer(), true, f_grad);

  StdVector<double> c_grad(n); // worst case

  optimization->GetDesign()->Reset(DesignElement::CONSTRAINT_GRADIENT, DesignElement::DEFAULT);
  for(unsigned int ci = 0; ci < m; ci++)
  {
    double         c    = c_vec[ci];
    double         y    = y_vec[ci];
    double         rho  = rho_vec[ci];
    Approximation* appr = constr[ci];
    Condition*     g    = optimization->constraints.view->Get(ci);
    c_grad.Resize(appr->jac_pattern.GetSize());
    EvalGradConstraint(g, 0, true, true, c_grad); // scale and normalize

    LOG_DBG3(feasPP) << "FP:CGAL ci=" << ci << " g=" << appr->ToString() << " c=" << c << " y=" << y << " rho=" << rho << " case=" << (-y / rho <= c);
    // derivative with respect to x
    for(unsigned int e = 0; e < n; e++)
    {
      grad[e] = f_grad[e];

      if(-y / rho <= c)
      {
        // we need the gradient of c with respect to e but have to consider the sparsity of c
        double grad_c = 0.0;
        if(c_grad.GetSize() == n)
          grad_c = c_grad[e];
        else // we need to traverse
          for(unsigned t = 0, tn = appr->jac_pattern.GetSize(); t < tn; t++)
            if(appr->jac_pattern[t] == e)
              grad_c = c_grad[t]; // sparse patterns are very short!

        grad[e] += y * grad_c + rho*c*grad_c;
        LOG_DBG3(feasPP) << "FP:CGAL ci=" << ci << " g=" << appr->ToString() << " e=" << e << " grad_c=" << grad_c << " -> " << grad[e];
      }
    }

    // derivative with respect to y
    grad[n+ci] = -y / rho <= c ? c : -y / rho;
    LOG_DBG3(feasPP) << "FP:CGAL ci=" << ci << " g=" << appr->ToString() << " dy -> " << grad[n+ci];
  }
  optimization->constraints.view->Done(); // reset slope constraint to global mode

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
  EvalConstraints(n, x_outer.GetPointer(), m, true, g_val.GetPointer(), true);
  // this is copy and paste from EvalGradConstraints()
  // reset values of the constraint gradients before the loop
  // as it also contains a loop over all the design elements
  optimization->GetDesign()->Reset(DesignElement::CONSTRAINT_GRADIENT, DesignElement::DEFAULT);
  for(unsigned int i = 0; i < m; i++)
  {
    constr[i]->outer_val = g_val[i];
    EvalGradConstraint(optimization->constraints.view->Get(i), 0, true, true, constr[i]->outer_grad); // scale and normalize

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

  // as the augmented Lagrangian assumes normalized constraints, everything is 0.0
  lower = 0.0;
  upper = 0.0;

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
      result = common->EvalConstraint(GetCondition(), true, true);
      break;
    case GRAD:
      common->EvalGradConstraint(GetCondition(), 0, true, true, *out); // scale and normalize
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
