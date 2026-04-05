#include "Optimization/Optimizer/MMA.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignElement.hh"
#include "MatVec/Matrix.hh"
#include "Utils/ToolsFull.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Optimization/Optimizer/BFGS.hh"

#include <string>
#include <limits>

DEFINE_LOG(mmaTopOpt, "mma")

using namespace CoupledField;
using std::pow;
using std::max;
using std::min;
using std::abs;
using std::string;

MMA::MMA(Optimization* opt, PtrParamNode pn) : BaseOptimizer(opt, pn, Optimization::MMA_SOLVER)
{

  asymUpdate.SetName("MMA::AsymUpdate");
  asymUpdate.Add(SVANBERG, "svanberg");
  asymUpdate.Add(TOPOPT_ROBUST_SHORT, "topopt_robust_short");
  asymUpdate.Add(TOPOPT_ROBUST_LONG, "topopt_robust_long");
  asymUpdate.Add(FIXED, "fixed");

  subSolverType.SetName("MMA::SubSolverType");
  subSolverType.Add(BFGS_SOLV, "bfgs");
  subSolverType.Add(IP_SOLV, "int_point");
  subSolverType.Add(NEWTON_SOLV, "newton");

  dual = &dual_solver;

  bool do_log = true; // default
  // defaults set in header
  if(this_opt_pn_)
  {
    /** Sub problem solver specification */
    max_sub_iter = this_opt_pn_->Get("max_sub_iter")->As<unsigned int>();
    sub_solve_tol = this_opt_pn_->Get("sub_solve_tol")->As<double>();
    subSolverType_ = subSolverType.Parse(this_opt_pn_->Get("sub_solver_type")->As<string>());
    if(subSolverType_ == IP_SOLV)
      dual = &primal_dual_solver;

    verbose_dual_vars = this_opt_pn_->Get("verbose_dual_vars")->As<bool>();

    penalty_c = this_opt_pn_->Get("constraint_penalty_c")->As<double>();

    kappa = this_opt_pn_->Get("svanbergs_kappa")->As<bool>();

    if(this_opt_pn_->Get("move_limit")->As<string>() == "tuned")
    {
      tuned = new Tuned(this_opt_pn_->Get("tuned", ParamNode::PASS), &move_limit, 5, 5, this);
      move_limit = tuned->max; // set by Tuned
    }
    else
      move_limit = this_opt_pn_->Get("move_limit")->As<double>();

    do_log = this_opt_pn_->Get("log")->As<bool>();

    if(this_opt_pn_->Has("hessian_corr"))
      dual->hess_corr = this_opt_pn_->Get("hessian_corr")->As<bool>();
    else
      dual->hess_corr = subSolverType_ == IP_SOLV; // only here as default (and possible)
    if(dual->hess_corr && subSolverType_ == NEWTON_SOLV)
      throw Exception("Hessian correction not possible or pure Hessian solver - choose IP-Solver");

    /*globallyConvergent = this_opt_pn_->Get("globally_convergent")->As<bool>();*/
    globallyConvergent = this_opt_pn_->Has("globally_convergent") && this_opt_pn_->Get("globally_convergent/enable")->As<bool>();
    if(globallyConvergent)
      rho_init = this_opt_pn_->Get("globally_convergent/rho")->As<double>();

    // asymptotes specification
    if(this_opt_pn_->Has("asymptotes"))
    {
      // is required element in schema
      PtrParamNode asym = this_opt_pn_->Get("asymptotes");

      asymUpdate_ = asymUpdate.Parse(asym->Get("update")->As<string>());

      /** Determines how aggressively the asymptotes are moved
       * asymptotes initialization and increase/decrease
       * these values are used in MMA::GenreteSubProblem() to update low and upp */
      asyminit = asym->Get("initial")->As<double>();
      asymdec = asym->Get("dec")->As<double>();
      asyminc = asym->Get("inc")->As<double>();
      ml_asym = asym->Get("ml_asym")->As<double>();

      if(asymUpdate_ == FIXED) {
        if(!asym->Has("fixed"))
          throw Exception("for MMA asymptotes update mode 'fixed' give subelement 'fixed'");
        asym_fixed_lower = asym->Get("fixed/lower")->As<double>();
        asym_fixed_upper = asym->Get("fixed/upper")->As<double>();
        upperMultiplier = asym->Get("fixed/upper_multiplier")->As<bool>();
        lowerMultiplier = asym->Get("fixed/lower_multiplier")->As<bool>();
        if(asym_fixed_upper <= asym_fixed_lower)
          throw Exception("upper fixed MMA asymptote needs to be larger than lower one");
      }

    } // end asymptotes
    if(subSolverType_ == BFGS_SOLV)
    {
      if(this_opt_pn_->Has("bfgs_setting"))
      {
        PtrParamNode bfgs_set = this_opt_pn_->Get("bfgs_setting");
        {
          bfgs.dual_init = bfgs_set->Get("initial")->As<double>();
          bfgs.dual_low = bfgs_set->Get("low")->As<double>();
          bfgs.dual_upp = bfgs_set->Get("upp")->As<double>();
          bfgs.nsmax= bfgs_set->Get("no_stores")->As<unsigned int>();
        }
      }
    }
  } // end pn reading

  if(do_log)
  {
    log = new std::ofstream(progOpts->GetSimName() + ".mma");
    if(log == nullptr)
       throw Exception("cannot open log file '" + progOpts->GetSimName()  + ".mma' for writing");
  }
  gsp_timer_ = info_->Get(ParamNode::SUMMARY)->Get("mma_generate_sub_prob/timer")->AsTimer(opt_timer).get();

  sps_timer_ = info_->Get(ParamNode::SUMMARY)->Get("mma_solve_sub_prob/timer")->AsTimer(opt_timer).get();

  WriteMMALogHeader();
}

MMA::~MMA()
{
  if(log) { delete log; log = nullptr; }
}

void MMA::ToInfo(PtrParamNode pn)
{
  BaseOptimizer::ToInfo(pn);

  pn->Get("globalConvergent")->SetValue(globallyConvergent);
  pn->Get("kappa")->SetValue(kappa);
  pn->Get("move_limit")->SetValue(tuned ? "tuned" : std::to_string(move_limit));
  if(tuned)
    tuned->ToInfo(pn->Get("tuned_move_limit"));
  PtrParamNode asym = pn->Get("asymptotes");
  asym->Get("update")->SetValue(asymUpdate.ToString(asymUpdate_));
  if(asymUpdate_ == FIXED) {
    asym->Get("lower")->SetValue(asym_fixed_lower);
    asym->Get("lower_is_multiplier")->SetValue(lowerMultiplier);
    asym->Get("upper")->SetValue(asym_fixed_upper);
    asym->Get("upper_is_multiplier")->SetValue(upperMultiplier);

    double low_min = low.Min();
    double low_max = low.Max();
    double upp_min = upp.Min();
    double upp_max = upp.Max();

    if(low_min != low_max)
      asym->Get("L")->SetValue(std::to_string(low_min) + " ... " + std::to_string(low_max));
    else
      asym->Get("L")->SetValue(low_min);

    if(upp_min != upp_max)
      asym->Get("U")->SetValue(std::to_string(upp_min) + " ... " + std::to_string(upp_max));
    else
      asym->Get("U")->SetValue(upp_min);
  } else {
    asym->Get("initial_delta")->SetValue(asyminit);
    asym->Get("inc")->SetValue(asyminc);
    asym->Get("dec")->SetValue(asymdec);
  }

  PtrParamNode ss = pn->Get("subsolver");
  ss->Get("type")->SetValue(subSolverType.ToString(subSolverType_));
  ss->Get("tol")->SetValue(sub_solve_tol);
  ss->Get("maxit")->SetValue(max_sub_iter);
}

void MMA::PostInit()
{
  PostInitScale(1.0);

  BaseOptimizer::PostInit();
  //  assert(optimization->objectives.data.GetSize() == 1); // trivial case only
  ConditionContainer& cc = optimization->constraints;
  // Physics relted initilization
  n = optimization->GetDesign()->GetNumberOfVariables();
  m = cc.view->GetNumberOfActiveConstraints();
  LOG_DBG(mmaTopOpt) << "PI: n=" << n << " m=" << m;

  if(m == 0)
    throw Exception("MMA requires at least one constraint!");

  // set design and bounds
  xval.Resize(n);
  optimization->GetDesign()->WriteDesignToExtern(xval, true);
  optimization->GetDesign()->WriteBoundsToExtern(std_xmin,std_xmax);
  // creates a reference to std_xmin in xmin, the memory still belongs to std_xmin
  xmin.Replace(std_xmin.GetPointer(), std_xmin.GetSize());
  xmax.Replace(std_xmax.GetPointer(), std_xmax.GetSize());

  grad_objective.Resize(n,0.0);
  constraints.Resize(m,0.0);
  grad_constraints.Resize(n*m,0.0);

  //MMA related initilization
  low.Resize(n,0.0); //cdev important to set it to zero
  upp.Resize(n,0.0); //cdev important to set it to zero

  if(asymUpdate_ == FIXED) {
    xold1 = xval;
    xold2 = xval;
  }

  if(globallyConvergent) {
    xold1 = xval;
    xold2 = xval;
    rho_0 = rho_init;
    rho.Resize(m,rho_init);
  }

  alpha.Resize(n);
  beta.Resize(n);
  change.Resize(n);
  p_0j.Resize(n); // for the objective function
  q_0j.Resize(n); // for the objective function
  p_ij.Resize(m, n); // for all the constraints
  q_ij.Resize(m, n); // for all the constraints
  b.Resize(m); // rhs of constraint inequality in subproblem.
  y.Resize(m); // elastic variable
  // a.Resize(m, 0.0); // Aage - min-max handling
  c.Resize(m, penalty_c); // cdev: elastic terms in subproblem

  switch(subSolverType_)
  {
  case BFGS_SOLV:
    bfgs.Init(this);
    break;
  case IP_SOLV:
    primal_dual_solver.Init(this);
    break;
  case NEWTON_SOLV:
    dual_solver.Init(this);
    break;
  }

  opt_timer->Stop();
}


void MMA::DescribeProperties(StdVector<std::pair<std::string, std::string> >& map) const
{
  map.Push_back(std::make_pair("move_limit", std::to_string(move_limit)));
}


void MMA::ComputeObjectiveConstraintsSensitivities()
{
  // To visualize the gradients
  DesignSpace* space = optimization->GetDesign();

  // evaluate properties of initial design
  obj_val = EvalObjective(n, xval.GetPointer(), true);
  EvalGradObjective(n, xval.GetPointer(), true, grad_objective);
  assert(grad_objective.GetSize() == n);
  EvalConstraints(n, xval.GetPointer(), m, true, constraints.GetPointer(), true);

  EvalGradConstraints(n, xval.GetPointer(), m, n*m, true, true, grad_constraints);

  // this is based on topopt implementation, since i have used the same parameters for a,c, asyminit as them, we have to scale similarly
  int res_idx_grad_0 = space->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::MMA_OBJ_GRADIANT);
  for(unsigned int ni=0; ni<n; ++ni)
  {
    if(res_idx_grad_0 >= 0)
    {
      // in a DesignElement we can store special results, but e.g. a slack variable is only a BaseDesignElement
      DesignElement* de = dynamic_cast<DesignElement*>(space->GetDesignElement(ni));
      if(de != NULL && res_idx_grad_0 >= 0)
        de->specialResult[res_idx_grad_0] = grad_objective[ni];
    }
  }

  int res_idx_grad_1 = space->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::MMA_CON_GRADIANT_1);
  int res_idx_grad_2 = space->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::MMA_CON_GRADIANT_2);
  for(unsigned int ni =0; ni<n; ++ni)
  {
    if(res_idx_grad_1 >= 0)
    {
      // in a DesignElement we can store special results, but e.g. a slack variable is only a BaseDesignElement
      DesignElement* de = dynamic_cast<DesignElement*>(space->GetDesignElement(ni));
      if(de != NULL && res_idx_grad_1 >= 0)
        de->specialResult[res_idx_grad_1] = grad_constraints[ni];
    }
    if(res_idx_grad_2 >= 0 && m > 1)
    {
      // in a DesignElement we can store special results, but e.g. a slack variable is only a BaseDesignElement
      DesignElement* de = dynamic_cast<DesignElement*>(space->GetDesignElement(ni));
      if(de != NULL && res_idx_grad_2  >= 0)
        de->specialResult[res_idx_grad_2] = grad_constraints[ni+n];
    }
  }

  LOG_DBG(mmaTopOpt) << "COCS: obj_val=" << obj_val;
  LOG_DBG2(mmaTopOpt) << "COCS: grad_objective=" << grad_objective.ToString();
  LOG_DBG(mmaTopOpt) << "COCS: constraints=" << constraints.ToString();
  LOG_DBG2(mmaTopOpt) << "COCS: grad_constraints=" << grad_constraints.ToString();

}

void MMA::SetPrimalVariables(const Vector<double>& x_in, const Vector<double>& y_in)
{
  this->xval = x_in;
  this->y = y_in;
}

void MMA::SolveProblem()
{
  if(optimization->GetMaxIterations() == 0)
    throw Exception("maximum number of iterations is 0 ");

  ComputeObjectiveConstraintsSensitivities();

  assert(optimization->GetCurrentIteration() == 0);
  CommitIteration(); // write initial configuration
  WriteMMALogMajor();

  int maxit = optimization->GetMaxIterations();

  /** ok is used to check if the sub problem was solvable*/
  bool ok = true;
  while(!optimization->DoStopOptimization() && optimization->GetCurrentIteration() <= maxit && ok)
  {
    ok = SolveSubProblem();

    // new design is stored, also the corresponding function values. Increments iteration
    if(ok)
    {
      ComputeObjectiveConstraintsSensitivities();
      CommitIteration();
      WriteMMALogMajor();
    }
  }

  if(optimization->GetCurrentIteration() >= maxit)
    optimization->DoStopOptimizationHelper(false, "maximum number of iterations reached");
  if(!ok)
    optimization->DoStopOptimizationHelper(false, mma_error);
}

/** Based on TopOpt implementation */

bool MMA::SolveSubProblem()
{
  /** Generate subproblem */
  gsp_timer_->Start();
  SetupSubProblem();// Generate subproblem
  gsp_timer_->Stop();

  /** Copy old design values
   *  will be used to modify asymptotes, so useless when we have fixed asymptotes*/
  if(!(asymUpdate_ == FIXED)) {
    xold2 = xold1;
    xold1 = xval;
  }

  LOG_DBG3(mmaTopOpt) << "S_MMA before: c=" << c.ToString();

  // Solve the MMA subproblem dual with interior point method
  sps_timer_->Start();
  switch(subSolverType_)
  {
  case IP_SOLV:
    sub_prob_iter = primal_dual_solver.IPSubProblemSolver();
    break;
  case NEWTON_SOLV:
    sub_prob_iter = dual_solver.FallbackHessianSolver();
    break;
  case BFGS_SOLV:
    sub_prob_iter = bfgs.BFGSSubProblemSolver();
    break;
  }
  sps_timer_->Stop();

  LOG_DBG3(mmaTopOpt) << "S_MMA before: c=" << c.ToString();

  return sub_prob_iter >= 0; // nothing encoded in a negative iteration number
}

void MMA::UpdateGCAsymptotes()
{
  assert(globallyConvergent);

  // Chaitanya's stuff plus epsilon from Fabian
  if(optimization->GetCurrentIteration() < 1)
  {
    for(unsigned int in =0; in < n; ++in)
    {
      low[in] = min(xmin[in] - 1e-3, low[in]);
      upp[in] = max(xmax[in] + 1e-3, upp[in]);
    }
    LOG_DBG3(mmaTopOpt) << "GSP:GC rho_0= " << rho_0;
    LOG_DBG3(mmaTopOpt) << "GSP:GC rho= " << rho.ToString();
  }
  else
  {
    /** Description provided in K.Svanberg's DCAMM lecture notes section 6*/
    double objective_approx = 0.0; // This is compute the objective approximation @xnew
    /** Deal with the objective*/
    for(unsigned int ni =0; ni< n; ++ni)
    {
      objective_approx += p_0j[ni] /(upp[ni] - xval[ni]) + q_0j[ni] / (xval[ni] - low[ni]);
    }
    objective_approx -= objective_r;

    if(objective_approx < obj_val)
      rho_0 = 2.0*rho_0;

    StdVector<double> constraints_approx(m); // This is compute the constraints approximation @xnew
    /** Deal with the constraint*/
    // #pragma omp parallel for num_threads(CFS_NUM_THREADS)
    for(unsigned int j =0; j<m; ++j)
    {
      constraints_approx[j] = 0.0;
      for(unsigned int ni =0; ni<n; ++ni)
      {
        double asym = p_ij[j][ni] /(upp[ni] - xval[ni]) + q_ij[j][ni] / (xval[ni] - low[ni]);
        constraints_approx[j] += asym;
      }
    }
    for(unsigned int mj=0; mj<m;++mj)
      constraints_approx[mj] += -b[mj];

    for(unsigned int mj =0; mj<m; ++mj)
      if(constraints_approx[mj] < constraints[mj])
        rho[mj] = 2.0*rho[mj];

    /** Decide if the low and upp should be updated. Description provided in K.Svanberg's DCAMM lecture notes section 6
     * We decide to update only if all the function_approximation evaluated at xnew is greater then function evaluated at xnew */
    bool shouldUpdate = true;
    if(!(objective_approx >= obj_val))
      shouldUpdate = false;

    if(shouldUpdate)
    {
      for(unsigned int mj=0; mj<m;++mj)
      {
        if( !(constraints_approx[mj] >= constraints[mj]))
        {
          shouldUpdate = false;
          break;
        }
      }
      for(unsigned int ni =0; ni< n; ++ni)
      {
        low[ni] = xval[ni] - (xold1[ni] - low[ni]);
        upp[ni] = xval[ni] + (xold1[ni] - upp[ni]);
      }
    }

    LOG_DBG2(mmaTopOpt) << "GSP:GC shouldUpdate= " << shouldUpdate;
    LOG_DBG2(mmaTopOpt) << "GSP:GC objective_approx(xnew)= " << objective_approx;
    LOG_DBG2(mmaTopOpt) << "GSP:GC objective(xnew)= " << obj_val;
    LOG_DBG2(mmaTopOpt) << "GSP:GC constraints_approx(xnew)= " << constraints_approx.ToString();
    LOG_DBG2(mmaTopOpt) << "GSP:GC constraints(xnew)= " << constraints.ToString();
    LOG_DBG2(mmaTopOpt) << "GSP:GC rho_0= " << rho_0;
    LOG_DBG2(mmaTopOpt) << "GSP:GC rho= " << rho.ToString();
  }
}

void MMA::UpdateNonGCAsymptotes()
{
  double sign_change = 0.0;
  double gamma = 0.0;

  DesignSpace* space = optimization->GetDesign();

  if(optimization->GetCurrentIteration() <= 2)
  {
    // We have to define the low and upp only once for fixed asymptotes.
    if(asymUpdate_ == FIXED)
    {
      for(unsigned int in =0; in < n; ++in)
      {
        BaseDesignElement* de = space->GetDesignElement(in);
        if(lowerMultiplier)
          low[in] = min(xmin[in]-1e-3, asym_fixed_lower*xmin[in]);
        else
          low[in] = min(de->GetLowerBound()-1e-3, asym_fixed_lower);
        if(upperMultiplier)
          upp[in] = max(xmax[in]+1e-3, asym_fixed_upper*xmax[in]);
        else
          upp[in] = max(de->GetUpperBound()-1e-3, asym_fixed_upper);
      }
    }
    else
    {
      for(unsigned int in =0; in < n; ++in)
      {
        low[in] = xval[in] -  asyminit*(xmax[in] - xmin[in]);
        upp[in] = xval[in] +  asyminit*(xmax[in] - xmin[in]);
      }
    }
  }
  if(asymUpdate_ != FIXED && optimization->GetCurrentIteration() > 2)
  {
    for(unsigned int i =0; i< n; ++i)
    {
      sign_change = (xval[i] - xold1[i]) * (xold1[i] - xold2[i]); // will be negative if there is oscillation of design values
      if(sign_change < 0.0)
        gamma = asymdec; // if there is oscillation decrease the asymptotes, convergence will be slower
      else if (sign_change > 0.0)
        gamma = asyminc; // if there is oscillation increase the asymptotes, for faster convergence.
      else
        gamma = 1.0; // there is no design change
      low[i] = xval[i] - gamma*(xold1[i] - low[i]);
      upp[i] = xval[i] + gamma*(upp[i] - xold1[i]);
      double x_max_min = max(1.0e-5, xmax[i] - xmin[i]);

      /** implementation of robust asymptote based on TopOpt code
       * refer function GenSub(...) in MMA.cc for the case RobustAsymptotesType == 1*/
      if (asymUpdate_ == TOPOPT_ROBUST_SHORT)
      {
        low[i] = max(low[i], xval[i]-10.0*x_max_min);
        low[i] = min(low[i], xval[i]-0.01*x_max_min);
        upp[i] = max(upp[i], xval[i]+0.01*x_max_min);
        upp[i] = min(upp[i], xval[i]+10.0*x_max_min);
      }
      /** implementation of robust asymptote based on TopOpt code
       * refer function GenSub(...) in MMA.c */
      else if(asymUpdate_ == TOPOPT_ROBUST_LONG)
      {
        //Robust Asymptotes
        low[i] = max(low[i], xval[i]-100.0*x_max_min);
        low[i] = min(low[i], xval[i]-1.0e-4*x_max_min);
        upp[i] = max(upp[i], xval[i]+1.0e-4*x_max_min);
        upp[i] = min(upp[i], xval[i]+100.0*x_max_min);

        x_max_min = xmin[i] - 1.0e-5;
        double x_max_max = xmax[i] + 1.0e-5;
        if(xval[i] < x_max_min)
        {
          low[i] = xval[i] - (x_max_max - xval[i])/0.9;
          upp[i] = xval[i] + (x_max_max - xval[i])/0.9;
        }
        if(xval[i] > x_max_max)
        {
          low[i] = xval[i] - (xval[i] - x_max_min)/0.9;
          upp[i] = xval[i] + (xval[i] - x_max_min)/0.9;
        }
      }
    }
  } // end iter >= 2
  LOG_DBG2(mmaTopOpt) << "UNGCA: iter=" << optimization->GetCurrentIteration() << " xmin=" << xmin.ToString();
  LOG_DBG2(mmaTopOpt) << "UNGCA: iter=" << optimization->GetCurrentIteration() << " xmax=" << xmax.ToString();
  LOG_DBG2(mmaTopOpt) << "UNGCA: iter=" << optimization->GetCurrentIteration() << " xval=" << xval.ToString();
  LOG_DBG2(mmaTopOpt) << "UNGCA: iter=" << optimization->GetCurrentIteration() << " L=" << low.ToString();
  LOG_DBG2(mmaTopOpt) << "UNGCA: iter=" << optimization->GetCurrentIteration() << " U=" << upp.ToString();

}


/** According to K.Svanberg's DCAMM lecture notes section 4.
 * can also refer to TopOpt code which is based on K.Svanberg's DCAMM lecture notes */
void MMA::SetupSubProblem()
{
  // Update of low and upp for globallyConvergent
  if(globallyConvergent)
    UpdateGCAsymptotes();
  else
    UpdateNonGCAsymptotes();

  // conditionally visualize the asymptotes
  DesignSpace* space = optimization->GetDesign();
  int res_idx_l = space->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::MMA_LOWER_VAL);
  int res_idx_u = space->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::MMA_UPPER_VAL);
  for(unsigned int i=0; (res_idx_l >= 0 || res_idx_u >=0) && i<n; i++)
  {
    // in a DesignElement we can store special results, but e.g. a slack variable is only a BaseDesignElement
    DesignElement* de = dynamic_cast<DesignElement*>(space->GetDesignElement(i));
    if(de != NULL && res_idx_l >= 0)
      de->specialResult[res_idx_l] = low[i];
    if(de != NULL && res_idx_u >= 0)
      de->specialResult[res_idx_u] = upp[i];
  }
  LOG_DBG2(mmaTopOpt) << "SSP: low=" << low.ToString();
  LOG_DBG2(mmaTopOpt) << "SSP: upp=" << upp.ToString();

  // Formation of pij and qij
  assert(!(globallyConvergent && kappa)); // not both together
  objective_r = -obj_val;  // not only globallyConvergent stuff
  double extra = 0;
  const double feps = 1.0e-6;
  for(unsigned int ni=0; ni < n; ++ni)
  {
    /** MMA and GCMMA – two methods for nonlinear optimization
     * this is chosen to avoid division by zero in subproblem.
     * or MMA and GCMMA – Fortran versions March 2013, Krister Svanberg, (2.8) and (2.9) */
    alpha[ni] = max(xmin[ni], max(low[ni]+ml_asym*(xval[ni]-low[ni]), xval[ni]-move_limit*(xmax[ni]-xmin[ni])));
    beta[ni]  = min(xmax[ni], min(upp[ni]+ml_asym*(upp[ni]-xval[ni]), xval[ni]+move_limit*(xmax[ni]-xmin[ni])));

    double dfdx_pos = max(0.0, grad_objective[ni]);
    double dfdx_neg = max(0.0, -1.0*grad_objective[ni]);

    if(globallyConvergent)
      extra = rho_0*(upp[ni] - low[ni])*0.5;
    else if(kappa)
      extra = 0.001*abs(grad_objective[ni]) + 0.5*feps/(upp[ni] - low[ni]);

    p_0j[ni] = pow(upp[ni] - xval[ni], 2.0) * (dfdx_pos + extra);
    q_0j[ni] = pow(xval[ni] - low[ni], 2.0) * (dfdx_neg + extra);

    objective_r += p_0j[ni] /(upp[ni] - xval[ni]) + q_0j[ni] / (xval[ni] - low[ni]);

    for(unsigned int mj =0; mj<m; ++mj)
    {
      dfdx_pos = max(0.0, grad_constraints[mj*n + ni]);
      dfdx_neg = max(0.0, -1*grad_constraints[mj*n + ni]);

      if(globallyConvergent)
        extra = rho[mj]*(upp[ni] - low[ni])*0.5;
      else if(kappa)
        extra = 0.001*abs(grad_constraints[mj*n + ni]) + 0.5*feps/(upp[ni] - low[ni]);
        // when constraintModification = true p_ij and q_ij are formed according to
        // description provided in K.Svanberg's DCAMM lecture notes section

      p_ij[mj][ni] = pow(upp[ni] - xval[ni], 2.0) * (dfdx_pos + extra);
      q_ij[mj][ni] = pow(xval[ni] - low[ni], 2.0) * (dfdx_neg + extra);
    }
  }

  // if also lambda is zero we cannot construct primal from dual
  assert(p_0j.NormL2() + q_0j.NormL2() > 0.0);

  // Calculation of RHS of the constraints in subproblem
  for(unsigned int j =0; j<m; ++j)
  {
    b[j] = -constraints[j];
    for(unsigned int ni =0; ni<n; ++ni)
      b[j] += p_ij[j][ni] /(upp[ni] - xval[ni]) + q_ij[j][ni] / (xval[ni] - low[ni]);
  }

  LOG_DBG2(mmaTopOpt) << "SSP: b=" << b.ToString();
  LOG_DBG2(mmaTopOpt) << "SSP: objective_r=" << objective_r << " obj_val=" << obj_val;
  LOG_DBG2(mmaTopOpt) << "SSP: alpha=" << alpha.ToString();
  LOG_DBG2(mmaTopOpt) << "SSP: beta=" << beta.ToString();
}

double MMA::EvalDualFunction(Vector<double> &lam)
{
  // To store the primal values
  Vector<double> x_d(n);
  Vector<double> y_d(m);
  // double z_d=0; from Aage min-max handling
  PrimalVarFromDualVar(lam, x_d, y_d); // Compute the new primal values based on the lam
  return FunctionOfDual(lam, x_d, y_d);
}


Vector<double> MMA::EvalDualGrad(Vector<double> &lam)
{
  assert(lam.GetSize() == m && "Size of lam is not equal to number of constraints");

  // To store the primal values
  Vector<double> x_d(n);
  Vector<double> y_d(m);
  // double z_d=0; Aage min-max handling
  PrimalVarFromDualVar(lam, x_d, y_d); // Compute the new primal values based on the lam and skip Aage min-max handling z

  Vector<double> grad(m);

  for(unsigned int jm = 0; jm < m; ++jm)
  {
    grad[jm] = 0.0;
    for(unsigned int in =0; in<n; ++in)
    {
      grad[jm] += p_ij[jm][in] / (upp[in] - x_d[in]) + q_ij[jm][in] / (x_d[in] - low[in]);
    }
    grad[jm] += -b[jm];
  }
  for(unsigned int jm=0; jm<m; ++jm)
    grad[jm] += - y_d[jm]; // - a[jm]*z_d Aage min-max handling

  // since this value will be used to maximize the lagrang multiplier,
  // we send the negative value to minimizer
  for(unsigned int jm=0; jm<m; ++jm)
    grad[jm] = - grad[jm];

  LOG_DBG3(mmaTopOpt) << "BFGS: lam=" << lam.ToString();
  LOG_DBG3(mmaTopOpt) << "BFGS: evalGrad=" << grad.ToString();

  return grad;
}




void MMA::PrimalVarFromDualVar(const Vector<double>& lambda_in, Vector<double>& x_out, Vector<double> &y_out) const
{
  assert(lambda_in.GetSize() == m);
  assert(x_out.GetSize() == n);
  assert(y_out.GetSize() == m);
  LOG_DBG3(mmaTopOpt) << "PVFDV: lam=" << lambda_in.ToString();
  LOG_DBG3(mmaTopOpt) << "PVFDV: p_0j=" << p_0j.ToString();
  LOG_DBG3(mmaTopOpt) << "PVFDV: q_0j=" << q_0j.ToString();

  // double lambda_x_a = 0.0; Aage min-max handling
  for(unsigned int mj=0; mj < m; ++mj)
  {
    assert(lambda_in[mj] >= 0.0);
    y_out[mj] = max(0.0, lambda_in[mj] - c[mj]);
    // lambda_x_a += lambda_in[mj]*a[mj]; Aage min-max handling
  }
  // z_out = max(0.0, 10.0*(lambda_x_a - 1.0)); Aage min-max handling

  // update of xval[] according to K.Svanberg's paper section 4. equation 17-19.
  double pj_x_lambda = 0.0;
  double qj_x_lambda = 0.0;
//#pragma omp parallel for num_threads(CFS_NUM_THREADS) reduction(+:pj_x_lambda,qj_x_lambda)
  for(unsigned int j = 0; j < n; ++j)
  {
    pj_x_lambda = p_0j[j];
    qj_x_lambda = q_0j[j];
    for(unsigned int i=0; i < m; ++i)
    {
      pj_x_lambda += p_ij[i][j]*lambda_in[i];
      qj_x_lambda += q_ij[i][j]*lambda_in[i];
    }

    x_out[j] = (sqrt(pj_x_lambda)*low[j] + sqrt(qj_x_lambda)*upp[j]) / (sqrt(pj_x_lambda) + sqrt(qj_x_lambda));
    assert(!std::isnan(x_out[j]));

    if(x_out[j] < alpha[j])
      x_out[j] = alpha[j];
    if(x_out[j] > beta[j])
      x_out[j] = beta[j];
  }
}

double MMA::FunctionOfDual(const Vector<double> &lam, const Vector<double>& x_in, const Vector<double>& y_in) const
{
  // no_sub_prb_eval++;
  double sum=0;
  for(unsigned int jn=0; jn < n; ++jn)
  {
    double pj_x_lambda = 0.0;
    double qj_x_lambda = 0.0;
    for(unsigned int im=0; im < m; ++im)
    {
      pj_x_lambda += p_ij[im][jn]*lam[im];
      qj_x_lambda += q_ij[im][jn]*lam[im];
    }
    sum += ((p_0j[jn] + pj_x_lambda)/(upp[jn] - x_in[jn])) + ((q_0j[jn] + qj_x_lambda)/(x_in[jn] - low[jn]));
  }

  //objective_r; // r_0 in Svanberg (20)

  assert(b.GetSize() == m && y_in.GetSize() == m && c.GetSize() == m && lam.GetSize() == m);

  double lagrange = -1* b.Inner(lam) + y_in.Inner(c) +.5 * y_in.Inner(y_in) - y_in.Inner(lam); // - lam[im]*a[im]*z_d Aage min-max handling


  double val = -objective_r + sum + lagrange; // + z_d + 0.5*z_d*z_d; Aage min-max handing
  //double val = sum + lagrange; // + z_d + 0.5*z_d*z_d; Aage min-max handing

  LOG_DBG3(mmaTopOpt) << "FOD: b=" << b.ToString() << " y=" << y_in.ToString() << " r0=" << objective_r
                      << " lam=" << lam.ToString() << " sum=" << sum << " l=" << lagrange << " -> val " << val;


  // since this value will be used to maximize the Lagrange multiplier,
  // we send the negative value to minimizer
  return -val;
}

/** according to N.Aage paper equation 15*/
void MMA::GradientOfDual(const Vector<double>& x_in, const Vector<double>& y_in, Vector<double>& dual_grad_out) const
{
  for(unsigned int jm = 0; jm < m; ++jm)
  {
    dual_grad_out[jm] = 0.0;
    for(unsigned int in =0; in<n; ++in)
      dual_grad_out[jm] += p_ij[jm][in] / (upp[in] - x_in[in]) + q_ij[jm][in] / (x_in[in] - low[in]);
  }

  for(unsigned int jm=0; jm<m; ++jm)
    dual_grad_out[jm] += -b[jm] - y_in[jm]; // -a[jm]*z In Aage the min-max handling
}

/** according to N.Aage paper section 3.3 for primal dual solver (mu!)*/
void MMA::HessianOfDual(const Vector<double>& lambda_in, const Vector<double>& mu_in, const Vector<double>& x_in, Matrix<double>& hess_out) const
{
  StdVector<double> pij_qij(n*m); // ( Pij/(Uj - xj)^2 + Qij/(xj - Lj)^2 )
  StdVector<double> grad_funcA(n); // Gradient of MMA approximation of objective
  double pj_x_lambda = 0.0;
  double qj_x_lambda = 0.0;

  LOG_DBG3(mmaTopOpt) << "HOD: l=" << lambda_in.ToString() << " mu=" << mu_in.ToString();


  #pragma omp parallel for num_threads(CFS_NUM_THREADS) reduction(+:pj_x_lambda,qj_x_lambda)
  for(int i=0; i < (int) n; ++i)
  {
    pj_x_lambda = p_0j[i];
    qj_x_lambda = q_0j[i];
    for(unsigned int j = 0; j<m; ++j)
    {
      pj_x_lambda += p_ij[j][i]*lambda_in[j];
      qj_x_lambda += q_ij[j][i]*lambda_in[j];
      pij_qij[i*m + j] = p_ij[j][i] / pow(upp[i] - x_in[i], 2.0) - q_ij[j][i]/pow(x_in[i] - low[i], 2.0);
    }
    grad_funcA[i] = -1.0 / (2.0 * pj_x_lambda / pow(upp[i] - x_in[i], 3.0) + 2.0*qj_x_lambda / pow(x_in[i] - low[i], 3.0));
    double tmp = (sqrt(pj_x_lambda)*low[i] + sqrt(qj_x_lambda)*upp[i]) / (sqrt(pj_x_lambda) + sqrt(qj_x_lambda));
    if(tmp < alpha[i])
      grad_funcA[i] = 0.0;
    if(tmp > beta[i])
      grad_funcA[i] = 0.0;
  }

  StdVector<double> tmp(n*m);

  for(unsigned int j = 0; j<m; ++j)
    for(unsigned int i = 0; i<n; ++i)
      tmp[j*n + i] = pij_qij[i*m + j] *  grad_funcA[i];

  hess_out.InitValue(0);
  for(unsigned int i=0; i< m; ++i)
    for(unsigned int j=0; j <m; ++j)
      for(unsigned int k=0; k<n; ++k)
        hess_out[i][j] += tmp[i * n + k] * pij_qij[k * m + j];

  // double lambda_x_a=0.0; Aage min-max handling
  for(unsigned int j=0; j<m; ++j)
  {
    assert(lambda_in[j] > 0.0);

    //lambda_x_a += lambda[j]*a[j];
    if(lambda_in[j] > c[j])
      hess_out[j][j] += -1.0;
    hess_out[j][j] += -mu_in[j]/lambda_in[j];
  }

  // Aage min-max handling
  // if(lambda_x_a > 0.0)
  //  {
  //  for(unsigned int j=0; j<m; ++j)
  //  {
  //    for(unsigned int k=0; k<m; ++k)
  //      dual_hessian[j*m+k] += -10.0*a[j]*a[k];
  //  }
  //}

  if(dual->hess_corr)
  {
    double trace = hess_out.Trace();
    double corr = 1.0e-4 * trace / m;
    if(-1.0*corr < 1.0e-7)
      corr = -1.0e-7;
    for(unsigned int i=0; i<m; ++i)
      hess_out[i][i] += corr;
  }

  assert(!hess_out.ContainsInf());
  assert(!hess_out.ContainsNaN());
}

/** based on mma.py Dual.hess() */
void MMA::HessianOfDual(const Vector<double>& lambda_in, const Vector<double>& x_in, Matrix<double>& hess_out) const
{
  assert(!dual->hess_corr);

  // Aage paper after (15)
  Matrix<double> dh(m,n); // d_h_ij / d_x
  for(unsigned int i = 0; i < m; i++)
    for(unsigned int j = 0; j < n; j++)
       dh[i][j] = (p_ij[i][j] / Pow2(upp[j] - x_in[j])) - (q_ij[i][j] / Pow2(x_in[j] - low[j]));

  Vector<double> dLi(n); // only diagonal
  for(unsigned int j = 0; j < n; j++)
  {
    if(x_in[j] > alpha[j] && x_in[j] < beta[j])
    {
      double lp = 0;
      double lq = 0;
      for(unsigned int i = 0; i < m; i++) {
        lp += lambda_in[i] * p_ij[i][j];
        lq += lambda_in[i] * q_ij[i][j];
      }
      double tmp = 2*(p_0j[j] + lp) / Pow3(upp[j] - x_in[j]) +  2*(q_0j[j] + lq) / Pow3(x_in[j] - low[j]);
      dLi[j] = 1/tmp;
    }
  }

  // H = -dh (dL)^1 dh^T (13)
  Matrix<double> dLi_dhT(n,m);
  for(unsigned j = 0; j < n; j++)
    for(unsigned i = 0; i < m; i++)
      dLi_dhT[j][i] = dLi[j] * dh[i][j];

  hess_out = dh * dLi_dhT; // skip -1 as we would need to add it at the end again for maximization

  LOG_DBG3(mmaTopOpt) << "HOD: l=" << lambda_in.ToString() << " x_in=" << x_in.ToString();
  LOG_DBG3(mmaTopOpt) << "HOD: dh=" << dh.ToString();
  LOG_DBG3(mmaTopOpt) << "HOD: dLi=" << dLi.ToString();
  LOG_DBG3(mmaTopOpt) << "HOD: dLi_dhT=" << dLi_dhT.ToString();

  LOG_DBG3(mmaTopOpt) << "HOD: hess_out=" << hess_out.ToString();

  assert(!hess_out.ContainsInf());
  assert(!hess_out.ContainsNaN());
}


double MMA::DualResidual(const Vector<double>& lambda_in, const Vector<double>& mu_in, const Vector<double>& x_in, const Vector<double>& y_in, double epsi) const
{
  StdVector<double> res(2*m);
  for(unsigned int j=0; j<m; ++j)
  {
    res[j]=0.0;
    res[j+m]=0.0;
    for(unsigned int i=0; i<n; ++i)
      res[j] += p_ij[j][i]/(upp[i] - x_in[i]) + q_ij[j][i]/(x_in[i] - low[i]);
  }
  for(unsigned int j=0; j<m; ++j)
  {
    res[j] += -b[j] - y_in[j] + mu_in[j]; // - a[j]*z Aage min-max handling
    res[j+m] += mu_in[j]*lambda_in[j] - epsi;
  }

  double result = 0.0;
  for(unsigned int i=0; i<2*m; ++i)
    if(result < abs(res[i]))
      result = abs(res[i]);
  return result;
}


void MMA::LogFileHeader(Optimization::Log& log)
{
  BaseOptimizer::LogFileHeader(log);

  log.AddToHeader("sub_prb_itr");

  if(asymUpdate_ != FIXED) {
    log.AddToHeader("neg_asym_min");
    log.AddToHeader("neg_asym_max");
    log.AddToHeader("pos_asym_min");
    log.AddToHeader("pos_asym_max");
  }

  for(int i = 0; i < optimization->constraints.view->GetNumberOfActiveConstraints(); i++)
  {
    Condition* g = optimization->constraints.view->Get(i);
    // don't add for slope constraints
    if(!g->IsLocalCondition())
      log.AddToHeader("lambda_" + g->ToString());
  }
}


void MMA::LogFileLine(std::ofstream* out, PtrParamNode iteration)
{
  BaseOptimizer::LogFileLine(out, iteration);

  double low_min = low.Min();
  double low_max = low.Max();
  double upp_min = upp.Min();
  double upp_max = upp.Max();

  if(out)
    *out << " \t" << sub_prob_iter;

  iteration->Get("sub_prb_itr")->SetValue(sub_prob_iter);

  if(asymUpdate_ != FIXED) {
    if(out)
      *out << " \t" << low_min << " \t" << low_max << " \t" << upp_min << " \t" << upp_max;
    iteration->Get("L")->SetValue(std::to_string(low_min) + " ... " + std::to_string(low_max));
    iteration->Get("U")->SetValue(std::to_string(upp_min) + " ... " + std::to_string(upp_max));
  }

  const Vector<double> lambda = GetSubSolverLambda();
  for(int i = 0; i < optimization->constraints.view->GetNumberOfActiveConstraints(); i++)
  {
    Condition* g = optimization->constraints.view->Get(i);
    // don't add for slope constraints
    if(!g->IsLocalCondition())
    {
      if(out)
        *out << " \t" << lambda[i];
      iteration->Get("lamba_" + g->ToString())->SetValue(lambda[i]);
    }
  }

  switch(subSolverType_)
  {
  case NEWTON_SOLV:
    break; // TODO
  case IP_SOLV:
    primal_dual_solver.IPLogFileLine(iteration);
    break;
  case BFGS_SOLV:
    bfgs.bfgs->LogFileLine(iteration);
    break;

  }
}

void MMA::WriteMMALogHeader()
{
  if(log)
  {
    *log << "sub solver: " << subSolverType.ToString(subSolverType_) << std::endl;
    *log << "verbose dual vars: " << verbose_dual_vars << std::endl;
  }
}

void MMA::WriteMMALogMajor()
{
  if(!log)
    return;

  Optimization* o = optimization;
  // we write after CommitIteration to have observe values, so decrement!
  assert(o->GetCurrentIteration() > 0);
  *log << "major " << o->GetCurrentIteration()-1;
  *log << " -> cost=" << o->objectives.GetHistoryValue();
  if(o->constraints.active.GetSize() < 5)
    for(auto c : o->constraints.active)
      *log << " " << c->ToString() << "=" << c->GetValue();
  if(o->constraints.observe.GetSize() < 5)
    for(auto c : o->constraints.observe)
      *log << " " << c->ToString() << "=" << c->GetValue();
  *log << std::endl;

}


const Vector<double>& MMA::GetSubSolverLambda() const
{
  switch(subSolverType_)
  {
  case IP_SOLV:
  case NEWTON_SOLV:
    return dual->lambda;
  case BFGS_SOLV:
    return bfgs.lambda;
  }
  // stupid C++ wants a return :(
  assert(false);
  return bfgs.lambda;
}

BFGSHolder::~BFGSHolder()
{
  if(bfgs)
  {
    delete bfgs;
    bfgs = NULL;
  }
}

void BFGSHolder::Init(MMA* mma)
{
  this->mma = mma;
  bfgs = new BFGS(mma->m, mma->sub_solve_tol, mma->max_sub_iter, nsmax, mma);
  lambda.Resize(mma->m); // a copy of bfgs.x to serve as frontent
  lam_v.Resize(mma->m, dual_init);
  up_lam.Resize(mma->m, dual_upp);
  lo_lam.Resize(mma->m, dual_low);
  primal_x.Resize(mma->n);
  primal_y.Resize(mma->m);
}


int BFGSHolder::BFGSSubProblemSolver()
{
  int iter = bfgs->SolveBFGS(lam_v, up_lam, lo_lam);

  // Copy the lambda values back
  lambda = bfgs->x;

  //PrimalVarFromDualVar();
  mma->PrimalVarFromDualVar(lambda, primal_x, primal_y); // skip z from Aage min-max handling
  mma->SetPrimalVariables(primal_x, primal_y);

  LOG_DBG3(mmaTopOpt) << "B_SSP: lambda=" << lambda.ToString();
  LOG_DBG3(mmaTopOpt) << "B_SSP: xval=" << primal_x.ToString();
  LOG_DBG3(mmaTopOpt) << "B_SSP: y=" << primal_y.ToString();
  return iter;
}

void DualSolver::Init(MMA* mma)
{
  this->mma = mma;
  this->n = mma->n;
  this->m = mma->m;

  // dual variables
  lambda.Resize(m,15.0); // Lagrange multiplier
  mu.Resize(m,0.0); // Lagrange multiplier

  // primal variables
  x.Resize(n);
  y.Resize(m);
  y0.Resize(m, 0.0);
}


void DualSolver::Project(const Vector<double>& a, double beta, const Vector<double>& b, Vector<double>& out)
{
  assert(a.GetSize() == b.GetSize());
  assert(a.GetSize() == out.GetSize());

  for(unsigned int i = 0; i < a.GetSize(); i++)
    out[i] = std::max(0.0, std::min(a[i] + beta * b[i], 1e20));

  LOG_DBG3(mmaTopOpt) << "DS:P a=" << a.ToString() << " beta=" << beta << " b=" << b.ToString() << " out=" << out.ToString();
}

Vector<double> DualSolver::Project(const Vector<double>& a, double beta, const Vector<double>& b)
{
  assert(a.GetSize() == b.GetSize());
  Vector<double> out(a.GetSize());

  for(unsigned int i = 0; i < a.GetSize(); i++)
    out[i] = std::max(0.0, std::min(a[i] + beta * b[i], 1e20));

  LOG_DBG3(mmaTopOpt) << "DS:P a=" << a.ToString() << " beta=" << beta << " b=" << b.ToString() << " -> " << out.ToString();

  return out;
}

void DualSolver::WriteHessianMinorHeader()
{
  if(mma->log)
  {
    *(mma->log) << "#minor \t||l-P(l_t)|| \tl_t \t||l_t-l|| \tdir \t hess \tls steps \tstep" << std::endl;
  }
}


void DualSolver::WriteHessianMinorLog(int minor, double pgc_norm, const Vector<double>& lmbda, double delta, const Vector<double>& dir, bool do_hess, int ls_steps, double step)
{
  if(mma->log)
  {
    *(mma->log) << minor << " \t" << pgc_norm << " \t";
    *(mma->log) << (mma->verbose_dual_vars ? lmbda.ToString() : std::to_string(lmbda.NormL2()));
    *(mma->log) << " \t" << delta << " \t";
    *(mma->log) << (mma->verbose_dual_vars ? dir.ToString() : std::to_string(dir.NormL2()));
    *(mma->log) << " \t" << (do_hess ? "true" : "false") << " \t" << ls_steps << " \t" << step << std::endl;
  }
}

int DualSolver::FallbackHessianSolver()
{
  // Simple steepest descent for Hessian steps with projection at bounds
  // algorithm 5.4.1, based on matlab implementation
  // This is not a projected Hessian w.r.t. non-active constraints!

  double alpha = 1e-4;
  Vector<double> lt(m); // temporary lambda (xt in mma.py, xc is lambda)
  Vector<double> grad(m);
  Vector<double> dir(m);  // descent direction from Hessian or gradient
  Vector<double> old_dir;
  Matrix<double> H(m,m); // hesse matrix
  Matrix<double> chol;   // working space for cholesky factorization

  mma->PrimalVarFromDualVar(lambda, x, y);
  double fc = mma->FunctionOfDual(lambda, x,y0);
  LOG_DBG2(mmaTopOpt) << "DS:FHS y=" << y.ToString() << " x="  << x.ToString();
  LOG_DBG2(mmaTopOpt) << "DS:FHS l=" << lambda.ToString() << " fc=" << fc;

  mma->GradientOfDual(x, y0, grad);
  grad *= -1;
  mma->HessianOfDual(lambda, x, H);

  LOG_DBG2(mmaTopOpt) << "DS:FHS l=" << lambda.ToString() << " mu=" << mu.ToString() << " H=" << H.ToString() << " grad=" << grad.ToString();
  // we may not completely trust the cholesky check. Result might still be crap
  bool do_hess = H.CholeskySolveLapack(chol, dir, grad, false); // no excpetion on invalid matrix
  if(!do_hess)
    dir = grad; // when Hessian does not work, use steepest descent
  LOG_DBG2(mmaTopOpt) << "DS:FHS dh=" << do_hess << " dir=" << dir.ToString();

  // set lt
  Project(lambda, -1.0, dir, lt); // lt = self.project(xc - dc, low, up)
  // test if we are with a valid hessian at the lower bound
  if(do_hess && lt.Min() <= 0.0) {
    LOG_DBG2(mmaTopOpt) << "DS:FHS initial fallback do_hess=" << do_hess << " min lt=" << lt.Min();
    do_hess = false;
    dir = grad;
    Project(lambda, -1.0, dir, lt);
  }

  sub_prob_iter = 1;

  WriteHessianMinorHeader();

  double pgc = lambda.NormL2(lt); // || lambda - lt ||
  while(pgc > mma->sub_solve_tol && sub_prob_iter <  mma->max_sub_iter)
  {
    double step = 1.0;
    int ls_steps = 0;
    // we have a valid Project(lambda,-1*step,dir, lt);
    Project(lambda,-1*step,dir, lt); // lt = self.project(lambda - step*dir, low, up)
    LOG_DBG2(mmaTopOpt) << "lt=" << lt.ToString() << " dir=" << dir.ToString() << " step=" << step << " lambda=" << lambda.ToString();
    mma->PrimalVarFromDualVar(lt, x, y);
    double ft = mma->FunctionOfDual(lt, x, y0);
    double fgoal = fc - pow(lambda.NormL2(lt),2) * (alpha/step); // fgoal = fc - np.sum((xc - xt)**2) * (alpha/step)
    double fg_cp = fgoal; // fgoal_checkpoint in case we cancel non-trustworthy Hessian

    LOG_DBG2(mmaTopOpt) << "DS:FHS minor=" << sub_prob_iter << " pgc=" << pgc << " lt=" << lt.ToString()
                        << " s=" << step << " lc=" << lambda.ToString() << " dc=" << dir.ToString() << " ft=" << ft
                        << " fgoal=" << fgoal << " need_ls=" << (ft > fgoal + mma->sub_solve_tol/100);
    while(ft > fgoal + mma->sub_solve_tol/100)
    {
      // two reductions from Hessian are not trustworthy
      if(do_hess && step < .1)
      {
         dir = grad;
         step = 1;
         Project(lambda,-1*step,dir, lt);
         mma->PrimalVarFromDualVar(lt, x, y);
         ft = mma->FunctionOfDual(lt, x, y0);
         fgoal = fg_cp;
         do_hess = false;
         LOG_DBG2(mmaTopOpt) << "DS:FHS cancel Hessian line search";
      }
      else
      {
        step *= .1;
        Project(lambda,-1*step,dir, lt);
        mma->PrimalVarFromDualVar(lt, x, y);
        ft = mma->FunctionOfDual(lt, x, y0);
        fgoal = fc - pow(lambda.NormL2(lt),2) * (alpha/step);
        ls_steps++;
        LOG_DBG2(mmaTopOpt) << "DS:FHS line search step: s=" << step << " dh=" << do_hess << " ft=" << ft;
      }
    }

    double old_do_hess = do_hess;
    old_dir = dir;

    fc = ft;
    mma->PrimalVarFromDualVar(lt, x, y);
    LOG_DBG2(mmaTopOpt) << "DS:FHS after ls: lt=" << lt.ToString() << " primal=" << x.ToString();
    mma->GradientOfDual(x, y0, grad);
    grad *= -1;
    mma->HessianOfDual(lt, x, H);
    LOG_DBG2(mmaTopOpt) << "DS:FHS after ls: H=" << H.ToString() << " grad=" << grad.ToString();

    // In Python we check the Hessian by (ev > 1e-8).all() but EV is too expensive and pos. def. is not sufficient (ev can be very small)
    // So in case the cholesky factorization succeeeds, we still go back, if the projection
    // a) goes too close to zero
    // b) the next pgc overshoots
    bool do_hess = H.CholeskySolveLapack(chol, dir, grad, false);
    if(!do_hess)
      dir = grad;

    LOG_DBG2(mmaTopOpt) << "DS:FHS cholesky solve -> dh" << do_hess << " dir=" << dir.ToString() << " grad=" << grad.ToString() << " H=" << H.ToString();

    double old_pgc = pgc;
    Vector<double> tmp = Project(lt,-1.0,dir);
    pgc = lt.NormL2(tmp); // || xt - self.project(xt - dc, low, up) ||

    if(do_hess && (tmp.NormL2() < 1e-4 || pgc > 10*old_pgc))
    {
      LOG_DBG2(mmaTopOpt) << "DS:FHS abandon Hessian: proj=" << tmp.NormL2() << " pgc=" << pgc << " dir=" << dir.ToString() << " tmp=" << tmp.ToString();
      do_hess = false;
      dir = grad;
      Project(lt,-1.0,dir,tmp);
      pgc = lt.NormL2(tmp);
    }

    WriteHessianMinorLog(sub_prob_iter, pgc, lt, lambda.NormL2(lt), old_dir, old_do_hess, ls_steps, step);

    lambda = lt;
    sub_prob_iter++;
  }

  if(sub_prob_iter <= 1)
    WriteHessianMinorLog(sub_prob_iter, pgc, lt, lambda.NormL2(lt), old_dir, do_hess, 0, 1);

  mma->SetPrimalVariables(x, y);

  return sub_prob_iter;
}

void PrimalDualSolver::Init(MMA* mma)
{
  DualSolver::Init(mma);

  // slack variables for IP
  s.Resize(2*m, 0.0);
}



int PrimalDualSolver::IPSubProblemSolver()
{
  Vector<double> grad(m); // dual
  Matrix<double> hess(m,m); // dual
  Vector<double> dir(m); // descent direction from Hessian

  // I need the below help variable
  bool alreadyRelaxed = false;

  // prepare for logging
  subiters.Resize(0);
  for(unsigned int j=0; j < m; ++j)
  {
    lambda[j] = mma->c[j] * 0.5;
    mu[j] = 1.0;
  }
  double tol = mma->sub_solve_tol * sqrt(m+n);
  double epsi = 1.0;
  double err = 1.0;

  unsigned int loop; //

  while(epsi > tol )
  {
    loop = 0;
    while(err > 0.9*epsi && loop < mma->max_sub_iter)
    {
      loop++;
      // mma->usedSubPrbItr++;

      mma->PrimalVarFromDualVar(lambda, x, y);

      mma->GradientOfDual(x, y, grad);

      for(unsigned int j =0; j<m ; ++j)
        grad[j] = -1.0*grad[j] - epsi/lambda[j];

      mma->HessianOfDual(lambda,mu,x,hess);
      LOG_DBG2(mmaTopOpt) << "DS:IPSPS hess=" << hess.ToString();
      hess.DirectSolve(dir,grad);
      LOG_DBG2(mmaTopOpt) << "DS:IPSPS hess dir=" << dir.ToString();

      // copy dir to the first half of s
      for(unsigned int j=0;j<m;j++)
        s[j] = dir[j];

      for(unsigned int i=0;i<m;i++)
        s[m+i] = -mu[i]+epsi/lambda[i]-s[i]*mu[i]/lambda[i];

      DualLineSearch(); // New value of lambda and mu will be updated

      mma->PrimalVarFromDualVar(lambda, x, y);

      err = mma->DualResidual(lambda, mu, x, y, epsi);

      // keep for verbose output in info xml
      if(progOpts->DoDetailedInfo())
      {
        subiters.Push_back(SubInfo());
        subiters.Last().lambda = lambda;
        subiters.Last().mu = mu;
        subiters.Last().s = s;
        subiters.Last().err = err;
        subiters.Last().iter = loop;
        subiters.Last().epsi = epsi;
      }
    }
    if(loop >= mma->max_sub_iter)
    {
      if(!alreadyRelaxed)
        {epsi=epsi/0.5; alreadyRelaxed = true;}
      else
      {
        std::stringstream ss;
        ss << "MMA subproblem cannot be solved in " << loop << " sub-iterations. err=" << err << " epsilon=" << epsi << " tol=" << tol;
        mma->mma_error = ss.str();
        alreadyRelaxed = false;
        return false;
      }
    }
    epsi = epsi * 0.1;
  }

  mma->SetPrimalVariables(x, y);

  return sub_prob_iter;
}

/** according to N.Aage paper section 3.2*/
void PrimalDualSolver::DualLineSearch()
{
  double theta=1.005;

  for (unsigned int i=0;i<m;i++)
  {
    if (theta < -1.01*s[i]/lambda[i])
      theta = -1.01*s[i]/lambda[i];
    if (theta < -1.01*s[i+m]/mu[i])
      theta = -1.01*s[i+m]/mu[i];
  }
  theta = 1.0/theta;
  for(unsigned int i=0;i<m;i++)
  {
    lambda[i] = lambda[i] + theta*s[i];
    mu[i] = mu[i] + theta*s[i+m];
    assert(lambda[i] >= 0 && mu[i] >= 0);
  }
  LOG_DBG3(mmaTopOpt) << "SSP: theta=" << theta;
}


void PrimalDualSolver::IPLogFileLine(PtrParamNode iteration)
{
  if(progOpts->DoDetailedInfo())
  {
    PtrParamNode sub =iteration->Get("subproblem");
    for(unsigned int is =0; is < subiters.GetSize(); is++)
    {
      SubInfo& si = subiters[is];
      PtrParamNode ipn = sub->Get("iter", ParamNode::APPEND);
      ipn->Get("number")->SetValue(is);
      ipn->Get("iter")->SetValue(si.iter);
      ipn->Get("epsi")->SetValue(si.epsi);
      ipn->Get("err")->SetValue(si.err);
      assert(si.lambda.GetSize() == m);
      assert(si.s.GetSize() == 2 * m);
      for(unsigned int ci = 0; ci < m; ci++)
      {
        string gn = mma->optimization->constraints.view->Get(ci)->ToString();
        ipn->Get("lambda_" + gn)->SetValue(si.lambda[ci]);
        ipn->Get("mu_" + gn)->SetValue(si.mu[ci]);
        ipn->Get("s_grad_" + gn)->SetValue(si.s[ci]);
        ipn->Get("s_slack_" + gn)->SetValue(si.s[m + ci]);
      }
    }
  }
}

