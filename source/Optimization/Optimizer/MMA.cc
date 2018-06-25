#include "Optimization/Optimizer/MMA.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignElement.hh"
#include "MatVec/Matrix.hh"
#include "Utils/tools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ProgramOptions.hh"

#include <string>
#include <limits>

DECLARE_LOG(mmaTopOpt)
DEFINE_LOG(mmaTopOpt, "mma")

using namespace CoupledField;
using std::pow;
using std::max;
using std::min;
using std::abs;
using std::string;

MMA::MMA(Optimization* opt, PtrParamNode pn) : BaseOptimizer(opt, pn, Optimization::MMA_SOLVER)
{
  n = 0; // number of design variables
  m = 0; // number of oncstrains
  PostInitScale(1.0);
  if(this_opt_pn_ != NULL)
  {
    max_sub_iter = this_opt_pn_->Get("max_sub_iter")->As<unsigned int>();
    asyminit = this_opt_pn_->Get("asym_init")->As<double>();
    asymdec = this_opt_pn_->Get("asym_dec")->As<double>();
    asyminc = this_opt_pn_->Get("asym_inc")->As<double>();
    robustAsymptotes = this_opt_pn_->Get("robust_asymptote")->As<bool>();
    fixedAsymptotes = this_opt_pn_->Has("fixed_asymptotes") && this_opt_pn_->Get("fixed_asymptotes/enable")->As<bool>();
    if(fixedAsymptotes) {
      asym_fixed_lower = this_opt_pn_->Get("fixed_asymptotes/lower")->As<double>();
      asym_fixed_upper = this_opt_pn_->Get("fixed_asymptotes/upper")->As<double>();
      if(asym_fixed_upper <= asym_fixed_lower)
        throw Exception("upper fixed MMA asymptote needs to be larger than lower one");
    }
    if(fixedAsymptotes && robustAsymptotes)
      throw Exception("cannot have robust and fixed asymptotes concurrently for MMA");
  }

  assert(!(robustAsymptotes && fixedAsymptotes)); // cannot have both the asymptotes


  gsp_timer_ = info_->Get(ParamNode::SUMMARY)->Get("GenreteSubProblem/timer")->AsTimer().get();
  gsp_timer_->SetSub();

  sps_timer_ = info_->Get(ParamNode::SUMMARY)->Get("SolveSubProblem/timer")->AsTimer().get();
  sps_timer_->SetSub();

}

MMA::~MMA()
{

}

void MMA::ToInfo(PtrParamNode pn)
{
  PtrParamNode m = pn->Get("mma")->Get(ParamNode::HEADER);
  PtrParamNode asym = m->Get("asymptotes");
  asym->Get("fixed")->SetValue(fixedAsymptotes);
  if(fixedAsymptotes) {
    asym->Get("fixed_lower")->SetValue(asym_fixed_lower);
    asym->Get("fixed_upper")->SetValue(asym_fixed_upper);
  }
  asym->Get("robust")->SetValue(robustAsymptotes);
  asym->Get("initial_asymptote")->SetValue(asyminit);
  asym->Get("increase_asymptote")->SetValue(asyminc);
  asym->Get("decrease_asymptote")->SetValue(asymdec);
}

void MMA::PostInit()
{
  assert(optimization->objectives.data.GetSize() == 1); // trivial case only
  ConditionContainer& cc = optimization->constraints;
  // Physics relted initilization
  n = optimization->GetDesign()->GetNumberOfVariables();
  m = cc.view->GetNumberOfActiveConstraints();

  // set design and bounds
  xval.Resize(n);
  optimization->GetDesign()->WriteDesignToExtern(xval, true);
  optimization->GetDesign()->WriteBoundsToExtern(std_xmin,std_xmax);
  // creates a reference to std_xmin in xmin, the meomorz still belongs to std_xmin
  xmin.Replace(std_xmin.GetSize(), std_xmin.GetPointer(), false);
  xmax.Replace(std_xmax.GetSize(), std_xmax.GetPointer(), false);

  grad_compliance.Resize(n,0.0);
  constrains.Resize(m,0.0);
  grad_constrains.Resize(n*m,0.0);
  //MMA related initilization
  low.Resize(n,0.0); //cdev important to set it to zero
  upp.Resize(n,0.0); //cdev important to set it to zero
  xold1 = xval;
  xold2 = xval;
  alpha.Resize(n);
  beta.Resize(n);
  change.Resize(n);
  p_0j.Resize(n); // for the objective function
  q_0j.Resize(n); // for the objective function
  p_ij.Resize(m, n); // for all the constrain function
  q_ij.Resize(m, n); // for all the constrain function
  b.Resize(m); // rhs of constrain inequality in subproblem.
  lamda.Resize(m); // Lagrange multiplier
  mu.Resize(m); // Lagrange multiplier
  y.Resize(m); // elastic variable
  z = 0.0; // elastic variable
  a.Resize(m, 0.0); // cdev: terms in subproblem which i do not understand
  c.Resize(m, 1000.0); // cdev: terms in subproblem which i do not understand
  d.Resize(m, 0.0); // cdev: terms in subproblem which i do not understand
  s.Resize(2*m, 0.0); // slack variables
  dual_gradient.Resize(m, 0.0);
  dual_hessian.Resize(m*m, 0.0);

}

void MMA::ComputeObjectiveConstraintsSensitivities()
{
  // To viaualize the gradients
  DesignSpace* space = optimization->GetDesign();

  // evaluate properties of initial design
  compliance = EvalObjective(n, xval.GetPointer(), false);
  EvalGradObjective(n, xval.GetPointer(), false, grad_compliance);
  EvalConstraints(n, xval.GetPointer(), m, false, constrains.GetPointer(), true);
//  con[0] = con[0] - 0.5;

  EvalGradConstraints(n, xval.GetPointer(), m, n*m, false, false, grad_constrains);

  /*
  if(optimization->GetCurrentIteration() == 0)
  {
    obj_scale = 10.0/compliance;
  }
  */
  compliance = compliance*obj_scale;
  // this is based on topopt implementation, since i have used the same parameters for a,c, asyminit as them, we have to scale similarly
  int res_idx_grad_0 = space->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::MMA_OBJ_GRADIANT);
  for(unsigned int ni=0; ni<n; ++ni){
    //grad_compliance[ni] = grad_compliance[ni]*obj_scale;

    if(res_idx_grad_0 >= 0 )
    {
      // in a DesignElement we can store special results, but e.g. a slack variable is only a BaseDesignElement
      DesignElement* de = dynamic_cast<DesignElement*>(space->GetDesignElement(ni));
      if(de != NULL && res_idx_grad_0 >= 0)
        de->specialResult[res_idx_grad_0 ] = grad_compliance[ni];
    }
  }



  int res_idx_grad_1 = space->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::MMA_CON_GRADIANT_1);
  int res_idx_grad_2 = space->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::MMA_CON_GRADIANT_2);
  for(unsigned int ni =0; ni<n; ++ni)
  {
    if(res_idx_grad_1 >= 0 )
      {
        // in a DesignElement we can store special results, but e.g. a slack variable is only a BaseDesignElement
        DesignElement* de = dynamic_cast<DesignElement*>(space->GetDesignElement(ni));
        if(de != NULL && res_idx_grad_1 >= 0)
          de->specialResult[res_idx_grad_1 ] = grad_constrains[ni];
      }
    if(res_idx_grad_2 >= 0 && m > 1)
      {
        // in a DesignElement we can store special results, but e.g. a slack variable is only a BaseDesignElement
        DesignElement* de = dynamic_cast<DesignElement*>(space->GetDesignElement(ni));
        if(de != NULL && res_idx_grad_2  >= 0)
          de->specialResult[res_idx_grad_2 ] = grad_constrains[ni+n];
      }
  }


  LOG_DBG3(mmaTopOpt) << "COCS: compliance=" << compliance;
  LOG_DBG3(mmaTopOpt) << "COCS: obj_scale=" << obj_scale;
  LOG_DBG3(mmaTopOpt) << "COCS: grad_compliance=" << grad_compliance.ToString(0);
  LOG_DBG3(mmaTopOpt) << "COCS: constrains=" << constrains.ToString(0);
  LOG_DBG3(mmaTopOpt) << "COCS: grad_constrains=" << grad_constrains.ToString(0);

  }

void MMA::SolveProblem()
{
  if(optimization->GetMaxIterations() == 0)
  {
    std::cout << std::endl << "maximum number of iterations is 0 " << std::endl;
    return;
  }
  ComputeObjectiveConstraintsSensitivities();

  assert(optimization->GetCurrentIteration() == 0);

  //while(!optimization->DoStopOptimization() && optimization->GetCurrentIteration() <= optimization->GetMaxIterations())
  int maxit = optimization->GetMaxIterations();

  bool ok = true;
  while(!optimization->DoStopOptimization() && optimization->GetCurrentIteration() <= maxit && ok)
  {
    //  StupidTest(); return ;
    AdjustMoveLimits();
    ok = SolveMMA();

    // new design is stored, also the correspoding function values. Increments iteration
    optimization->CommitIteration();
    if(ok)
      ComputeObjectiveConstraintsSensitivities();
  }

  PtrParamNode summary = optimization->optInfoNode->Get(ParamNode::SUMMARY);
  summary->Get("break/converged")->SetValue(ok);

  if(!ok)
   summary->SetWarning(mma_error_);

}

void MMA::AdjustMoveLimits() // cdev: Tested no errors
{

  assert(xmin.GetSize() >= n);
  assert(xval.GetSize() >= n);
  assert(xmax.GetSize() >= n);

  DesignSpace* space = optimization->GetDesign();
  int res_idx_l = space->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::MMA_LOWER_VAL);
  int res_idx_u = space->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::MMA_UPPER_VAL);

  for(unsigned int i =0; i<n; ++i)
  {
    BaseDesignElement* de = space->GetDesignElement(i);

    if(!fixedAsymptotes){
      xmin[i] = max(de->GetLowerBound(), xval[i]-move_limits);
      xmax[i] = min(de->GetUpperBound(), xval[i]+move_limits);
    }

    if(res_idx_l >= 0 || res_idx_u >=0)
    {
      // in a DesignElement we can store special results, but e.g. a slack variable is only a BaseDesignElement
      DesignElement* de = dynamic_cast<DesignElement*>(space->GetDesignElement(i));
      if(de != NULL && res_idx_l >= 0)
        de->specialResult[res_idx_l] = xmin[i];
      if(de != NULL && res_idx_u >= 0)
        de->specialResult[res_idx_u] = xmax[i];
    }
  }
  LOG_DBG3(mmaTopOpt) << "AML: xmin=" << xmin.ToString(0);
  LOG_DBG3(mmaTopOpt) << "AML: xmax=" << xmax.ToString(0);
}

bool MMA::SolveMMA()
{
  gsp_timer_->Start();
  // Generate subproblem
  GenreteSubProblem();

  gsp_timer_->Stop();

  // Copy old design values
  xold2 = xold1;
  xold1 = xval;

  // Solve the MMA subproblem dual with interior point method
  sps_timer_->Start();
  bool ok = SolveSubProblem();
  sps_timer_->Stop();
  return ok;
}

void MMA::GenreteSubProblem() // cdev: Tested decimal errors in low and upp and subsequently in alpha and beta, may be dude to reding files
{
  double sign_change = 0.0;
  double gamma = 0.0;

  // For visualization
//  DesignSpace* space = optimization->GetDesign();
//  int res_idx = space->GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::MMA_ASYMPTOTE);

  if(optimization->GetCurrentIteration() < 3)
  {
    if(fixedAsymptotes){
      for(unsigned int i =0; i < n; ++i)
      {
        low[i] = -0.1;
        upp[i] = 10.0*xmax[i];
      }
    }
    else {
        for(unsigned int i =0; i < n; ++i)
        {
          low[i] += xval[i] -  asyminit*(xmax[i] - xmin[i]);
          upp[i] += xval[i] +  asyminit*(xmax[i] - xmin[i]);
        }
    }
  }
  if(optimization->GetCurrentIteration() > 2)
  {
    for(unsigned int i =0; i< n; ++i)
    {
      sign_change = (xval[i] - xold1[i]) * (xold1[i] - xold2[i]); // will be negative if there is oscillation of design values
      if( sign_change < 0.0) gamma = asymdec; // if there is oscillation decrese the asymptotes, convergence will be slower
      else if (sign_change > 0.0) gamma = asyminc; // if there is oscillation increase the asymptotes, for faster convergence.
      else gamma = 1.0; // there is no design change
      low[i] = xval[i] - gamma*(xold1[i] - low[i]);
      upp[i] = xval[i] + gamma*(upp[i] - xold1[i]);
      double x_min_tmp = max(1.0e-5, xmax[i] - xmin[i]);
      double x_max_tmp = 0.0;

      if(fixedAsymptotes) {
          low[i] = 0.0;
          upp[i] = 10.0*xmax[i];
      }
      else if(robustAsymptotes) {
        //Robust Asymptotes
        low[i]=max(low[i],xval[i]-100.0*x_min_tmp);
        low[i]=min(low[i],xval[i]-0.01*x_min_tmp);
        upp[i]=max(upp[i],xval[i]+0.01*x_min_tmp);
        upp[i]=min(upp[i],xval[i]+100.0*x_min_tmp);
        low[i] = min(max(low[i], xval[i] - 100.0*x_min_tmp ),xval[i] - 1.0e-4*x_min_tmp);
        upp[i] = min(max(upp[i], xval[i] + 1.0e-4*x_min_tmp ),xval[i] + 100.0*x_min_tmp);
        x_min_tmp = xmin[i] - 1.0e-5;
        x_max_tmp = xmax[i] + 1.0e-5;
        if(xval[i] < x_min_tmp)
        {
          low[i] = xval[i] - (x_max_tmp - xval[i])/0.9;
          upp[i] = xval[i] + (x_max_tmp - xval[i])/0.9;
        }
        if(xval[i] > x_max_tmp)
        {
          low[i] = xval[i] - (xval[i] - x_min_tmp)/0.9;
          upp[i] = xval[i] + (xval[i] - x_min_tmp)/0.9;
        }
      }
      else {
        low[i]=max(low[i],xval[i]-10.0*x_min_tmp);
        low[i]=min(low[i],xval[i]-1.0e-4*x_min_tmp);
        upp[i]=max(upp[i],xval[i]+1.0e-4*x_min_tmp);
        upp[i]=min(upp[i],xval[i]+10.0*x_min_tmp);
      }
    }
  }

  // Formation of pij and qij

  double dfdx_pos, dfdx_neg, extra;
  double feps = 1.0e-6;
  for(unsigned int ni=0; ni < n; ++ni)
  {
    alpha[ni] = max(xmin[ni], 0.9*low[ni]+0.1*xval[ni]);
    beta[ni] = min(xmax[ni], 0.9*upp[ni]+0.1*xval[ni]);
    dfdx_pos = max(0.0, grad_compliance[ni]);
    dfdx_neg = max(0.0, -1.0*grad_compliance[ni]);
    extra = 0.001*Abs(grad_compliance[ni]) + 0.5*feps/(upp[ni] - low[ni]);
    p_0j[ni] = pow(upp[ni] - xval[ni], 2.0) * (dfdx_pos + extra);
    q_0j[ni] = pow(xval[ni] - low[ni], 2.0) * (dfdx_neg + extra);

    for(unsigned int mj =0; mj<m ; ++mj)
    {
      if(mj == 0) // constrain modification version from topopt code
      {
        dfdx_pos = max(0.0, grad_constrains[mj*n + ni]);
        dfdx_neg = max(0.0, -1*grad_constrains[mj*n + ni]);

        extra = 0.001*Abs(grad_constrains[mj*n + ni]) + 0.5*feps/(upp[ni] - low[ni]);
        p_ij[mj][ni] = pow(upp[ni] - xval[ni], 2.0) * (dfdx_pos + extra) ;
        q_ij[mj][ni] = pow(xval[ni] - low[ni], 2.0) * (dfdx_neg + extra) ;
      }
      else
      {
        dfdx_pos = max(0.0, grad_constrains[mj*n + ni]);
        dfdx_neg = max(0.0, -1*grad_constrains[mj*n + ni]);
        p_ij[mj][ni] = pow(upp[ni] - xval[ni], 2.0) * dfdx_pos ;
        q_ij[mj][ni] = pow(xval[ni] - low[ni], 2.0) * dfdx_neg ;
      }
    }
  }
  // Calculation of the MMA approximation function
  for(unsigned int j =0; j<m; ++j)
  {
    b[j] = 0.0;
    for(unsigned int i =0; i<n; ++i)
    {
      double asym = p_ij[j][i] /(upp[i] - xval[i]) + q_ij[j][i] / (xval[i] - low[i]);
      b[j] += asym;
    }
  }
  for(unsigned int j=0; j<m;++j)
  {
    b[j] += -constrains[j];
  }

  LOG_DBG3(mmaTopOpt) << "GSP:Sub low=" << low.ToString(0);
  LOG_DBG3(mmaTopOpt) << "GSP:Sub upp=" << upp.ToString(0);
  LOG_DBG3(mmaTopOpt) << "GSP:Sub alpha=" << alpha.ToString(0);
  LOG_DBG3(mmaTopOpt) << "GSP:Sub beta=" << beta.ToString(0);
}

bool MMA::SolveSubProblem()
{
  // prepare for logging
  subiters.Resize(0);
  for(unsigned int j=0; j < m; ++j)
  {
    lamda[j] = c[j] * 0.5;
    mu[j] = 1.0;
  }
  double tol = 1.0e-9 * sqrt(m+n);
  double epsi = 1.0;
  double err = 1.0;

  unsigned int loop; //

  LOG_DBG3(mmaTopOpt) << "SSP: Problem Iteration = " << optimization->GetCurrentIteration() << " @@@@@@@@@@@@@@@@@@@@@@@@@@";
  LOG_DBG3(mmaTopOpt) << "SSP: xval=" << xval.ToString(0);
  LOG_DBG3(mmaTopOpt) << "SSP: y=" << y.ToString(0);
  LOG_DBG3(mmaTopOpt) << "SSP: z=" << z;
  LOG_DBG3(mmaTopOpt) << "SSP: low=" << low.ToString(0);
  LOG_DBG3(mmaTopOpt) << "SSP: alpha=" << alpha.ToString(0);
  LOG_DBG3(mmaTopOpt) << "SSP: upp=" << upp.ToString(0);
  LOG_DBG3(mmaTopOpt) << "SSP: beta=" << beta.ToString(0);
  LOG_DBG3(mmaTopOpt) << "SSP: p_0j=" << p_0j.ToString(0);
  LOG_DBG3(mmaTopOpt) << "SSP: q_0j=" << q_0j.ToString(0);
  LOG_DBG3(mmaTopOpt) << "SSP: p_ij=" << p_ij.ToString(0);
  LOG_DBG3(mmaTopOpt) << "SSP: q_ij=" << q_ij.ToString(0);
  LOG_DBG3(mmaTopOpt) << "SSP: conA=" << b.ToString(0);
  LOG_DBG3(mmaTopOpt) << "SSP: lamda=" << lamda.ToString(0);
  LOG_DBG3(mmaTopOpt) << "SSP: mu=" << mu.ToString(0);
  LOG_DBG3(mmaTopOpt) << "SSP: a=" << a.ToString(0);
  LOG_DBG3(mmaTopOpt) << "SSP: c=" << c.ToString(0);
  LOG_DBG3(mmaTopOpt) << "SSP: d=" << d.ToString(0);
  LOG_DBG3(mmaTopOpt) << "SSP: epsi=" << epsi;
  LOG_DBG3(mmaTopOpt) << "SSP: err=" << err;

  while(epsi > tol )
  {
    loop = 0;
    while(err > 0.9*epsi && loop < max_sub_iter)
    {
      ++loop;
      LOG_DBG3(mmaTopOpt) << "SSP: Sub Prob Loop = " << loop <<" :::::::::::::::::::::::::::::::::::::::::::::: ";
      LOG_DBG3(mmaTopOpt) << "SSP: loop=" << loop << " epsi=" << epsi << " tol=" << tol;

      PrimalVarFromDualVar();

      LOG_DBG3(mmaTopOpt) << "SSP: P_Fr_D xval=" << xval.ToString(0);
      LOG_DBG3(mmaTopOpt) << "SSP: P_Fr_D y=" << y.ToString(0);
      LOG_DBG3(mmaTopOpt) << "SSP: P_Fr_D z=" << z;

      GradientOfDual();

      LOG_DBG3(mmaTopOpt) << "SSP: Grad_D grad=" << dual_gradient.ToString(0);

      for(unsigned int j =0; j<m ; ++j)
      {
        dual_gradient[j] = -1.0*dual_gradient[j] - epsi/lamda[j];
      }

      LOG_DBG3(mmaTopOpt) << "SSP: Grad_Mod grad=" << dual_gradient.ToString(0);

      HessianOfDual();

      LOG_DBG3(mmaTopOpt) << "SSP: Hess_D Hess=" << dual_hessian.ToString(0);

      Factorize(dual_hessian, m);

      LOG_DBG3(mmaTopOpt) << "SSP: Factorize Hess=" << dual_hessian.ToString(0);

      Solve(dual_hessian, dual_gradient,m);

      LOG_DBG3(mmaTopOpt) << "SSP: Solve grad=" << dual_gradient.ToString(0);

      for (unsigned int j=0;j<m;j++){
        s[j]=dual_gradient[j];
      }
      for (unsigned int i=0;i<m;i++){
        s[m+i]= -mu[i]+epsi/lamda[i]-s[i]*mu[i]/lamda[i];
      }

      LOG_DBG3(mmaTopOpt) << "SSP: S=" << s.ToString(0);

      DualLineSearch(); // New value of lamda and mu will be updated

      LOG_DBG3(mmaTopOpt) << "SSP: lamda=" << lamda.ToString(0);
      LOG_DBG3(mmaTopOpt) << "SSP: mu=" << mu.ToString(0);

      PrimalVarFromDualVar();

      LOG_DBG3(mmaTopOpt) << "SSP: P_Fr_D xval=" << xval.ToString(0);
      LOG_DBG3(mmaTopOpt) << "SSP: P_Fr_D y=" << y.ToString(0);
      LOG_DBG3(mmaTopOpt) << "SSP: P_Fr_D z=" << z;

      err = DualResidual( epsi);

      LOG_DBG3(mmaTopOpt) << "SSP: err=" << err;
      LOG_DBG3(mmaTopOpt) << "SSP: --------------------------------------------------------------------------- \n \n" ;

      // keep for verbose output in info xml
      if(progOpts->DoDetailedInfo())
      {
        subiters.Push_back(SubInfo());
        subiters.Last().lambda = lamda;
        subiters.Last().mu = mu;
        subiters.Last().s = s;
        subiters.Last().err = err;
        subiters.Last().iter = loop;
        subiters.Last().epsi = epsi;
      }
    }

    if(loop >= max_sub_iter)
    {
      std::stringstream ss;
      ss << "MMA subproblem cannot be solved in " << loop << " sub-iters. err=" << err << " epsilon=" << epsi << " tol=" << tol;
      mma_error_ = ss.str();
      return false;
    }

    epsi=epsi*0.1;

    LOG_DBG3(mmaTopOpt) << "SSP: .......................................................................... \n \n" ;
  }

  return true;
}

void MMA::PrimalVarFromDualVar()
{
  double lamda_x_a = 0.0;
  for(unsigned int j=0; j < m; ++j)
  {
    if(lamda[j] < 0.0) lamda[j] = 0.0;
    y[j] = max(0.0, lamda[j] - c[j]);
    lamda_x_a += lamda[j]*a[j];
  }
  z = max(0.0, 10.0*(lamda_x_a - 1.0));
  double pj_x_lamda = 0.0;
  double qj_x_lamda = 0.0;
  for(unsigned int i = 0; i < n; ++i)
  {
    pj_x_lamda = p_0j[i];
    qj_x_lamda = q_0j[i];
    for(unsigned int j=0; j < m; ++j)
    {
      pj_x_lamda += p_ij[j][i]*lamda[j];
      qj_x_lamda += q_ij[j][i]*lamda[j];
    }
    xval[i] = (sqrt(pj_x_lamda)*low[i] + sqrt(qj_x_lamda)*upp[i]) / (sqrt(pj_x_lamda) + sqrt(qj_x_lamda));
    if(xval[i] < alpha[i]) xval[i] = alpha[i];
    if(xval[i] > beta[i]) xval[i] = beta[i];
  }
}

void MMA::GradientOfDual()
{
  for(unsigned int j = 0; j < m; ++j)
  {
    dual_gradient[j] = 0.0;
    for(unsigned int i =0; i<n; ++i)
    {
      dual_gradient[j] += p_ij[j][i] / (upp[i] - xval[i]) + q_ij[j][i] / (xval[i] - low[i]);
    }
  }
for(unsigned int j=0; j<m; ++j)
  {
    dual_gradient[j] += -b[j] - a[j]*z - y[j];
  }
}

void MMA::HessianOfDual()
{
  StdVector<double> pij_qij(n*m); // ( Pij/(Uj - xj)^2 + Qij/(xj - Lj)^2 )
  StdVector<double> grad_funcA(n); // Gradinet of MMA approximation of objective
  double pj_x_lamda = 0.0;
  double qj_x_lamda = 0.0;

  for(unsigned int i=0; i<n; ++i)
  {
    pj_x_lamda = p_0j[i];
    qj_x_lamda = q_0j[i];
    for(unsigned int j = 0; j<m; ++j)
    {
      pj_x_lamda += p_ij[j][i]*lamda[j];
      qj_x_lamda += q_ij[j][i]*lamda[j];
      pij_qij[i*m + j] = p_ij[j][i] / pow(upp[i] - xval[i], 2.0) - q_ij[j][i]/pow(xval[i] - low[i], 2.0);
    }
    grad_funcA[i] = -1.0 / (2.0 * pj_x_lamda / pow(upp[i] - xval[i], 3.0) + 2.0*qj_x_lamda / pow(xval[i] - low[i], 3.0));
    double tmp = (sqrt(pj_x_lamda)*low[i] + sqrt(qj_x_lamda)*upp[i]) / (sqrt(pj_x_lamda) + sqrt(qj_x_lamda));
    if(tmp < alpha[i]) {grad_funcA[i] = 0.0;}
    if(tmp > beta[i]) {grad_funcA[i] = 0.0;}
  }

  StdVector<double> tmp(n*m);

  for(unsigned int j = 0; j<m; ++j)
  {
    for(unsigned int i = 0; i<n; ++i)
    {
      tmp[j*n + i] = 0.0;
      tmp[j*n + i] += pij_qij[i*m + j] *  grad_funcA[i];
    }
  }

  for(unsigned int i=0; i< m; ++i)
  {
    for(unsigned int j=0; j <m; ++j)
    {
      dual_hessian[i*m + j] = 0.0;
      for(unsigned int k=0; k<n; ++k)
      {
        dual_hessian[i * m + j] += tmp[i * n + k] * pij_qij[k * m + j];
      }
    }
  }

  double lamda_x_a=0.0;
  for(unsigned int j=0; j<m; ++j)
  {
    if(lamda[j] < 0.0) { lamda[j] = 0.0; }
    lamda_x_a += lamda[j]*a[j];
    if(lamda[j] > c[j]) { dual_hessian[j*m+j] += -1.0; }
    dual_hessian[j*m+j] += -mu[j]/lamda[j];
  }

  if(lamda_x_a > 0.0)
  {
    for(unsigned int j=0; j<m; ++j)
    {
      for(unsigned int k=0; k<m; ++k)
        dual_hessian[j*m+k] += -10.0*a[j]*a[k];
    }
  }

  double dual_hessian_trace = 0.0;
  for(unsigned int i=0; i<m; ++i)
    dual_hessian_trace += dual_hessian[i*m + i];

  double dual_hessian_correction = 1.0e-4*dual_hessian_trace/(double)m;
  if(-1.0*dual_hessian_correction < 1.0e-7) {dual_hessian_correction = -1.0e-7;}
  for(unsigned int i=0; i<m; ++i)
    dual_hessian[i*m + i] += dual_hessian_correction;
}

void MMA::Factorize(StdVector<double> &K, const unsigned int nn)
{
  for (unsigned int  ss=0;ss<nn-1;ss++){
    for (unsigned int  i=ss+1;i<nn;i++){
      K[i*nn+ss] = K[i*nn+ss] / K[ss*nn+ss];
      for (unsigned int  j=ss+1;j<nn;j++){
        K[i*nn+j] = K[i*nn+j] - K[i*nn+ss]*K[ss*nn+j];
      }
    }
  }
}

void MMA::Solve(StdVector<double> &K, StdVector<double> &x, const int nn)
{
  for (int i=1;i<nn;i++){
    double a = 0.0;
    for (int j=0;j<i;j++){
      a = a - K[i*nn+j]*x[j];
    }
    x[i] = x[i] + a;
  }
  x[nn-1] = x[nn-1]/K[(nn-1)*nn+(nn-1)];
  for (int i=nn-2;i>=0;i--){
    double a=x[i];
    for (int j=i+1;j<nn;j++){
      a = a - K[i*nn+j]*x[j];
    }
    x[i] = a/K[i*nn+i];
  }
}

void MMA::DualLineSearch()
{
  double theta=1.005;
  for (unsigned int i=0;i<m;i++){
    if (theta < -1.01*s[i]/lamda[i]){
      theta = -1.01*s[i]/lamda[i];
    }
    if (theta < -1.01*s[i+m]/mu[i]){
      theta = -1.01*s[i+m]/mu[i];
    }
  }
  theta = 1.0/theta;
  for (unsigned int i=0;i<m;i++){
    lamda[i] = lamda[i] + theta*s[i];
    mu[i] = mu[i] + theta*s[i+m];
  }
  LOG_DBG3(mmaTopOpt) << "SSP: theta=" << theta;
}

double MMA::DualResidual(double epsi)
{
  StdVector<double> res(2*m);
  for(unsigned int j=0; j<m; ++j)
  {
    res[j]=0.0;
    res[j+m]=0.0;
    for(unsigned int i=0; i<n; ++i)
    {
      res[j] += p_ij[j][i]/(upp[i] - xval[i]) + q_ij[j][i]/(xval[i] - low[i]);
    }
  }
  for(unsigned int j=0; j<m; ++j)
  {
    res[j] += -b[j] - a[j]*z - y[j] + mu[j];
    res[j+m] += mu[j]*lamda[j] - epsi;
  }

  double result = 0.0;
  for(unsigned int i=0; i<2*m; ++i)
    if(result < Abs(res[i])) {result = Abs(res[i]);}
  return result;
}


void MMA::LogFileHeader(Optimization::Log& log)
{
  log.AddToHeader("neg_asym_min");
  log.AddToHeader("neg_asym_max");
  log.AddToHeader("pos_asym_min");
  log.AddToHeader("pos_asym_max");

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
  double xmin_min = xmin.Min();
  double xmin_max = xmin.Max();
  double xmax_min = xmax.Min();
  double xmax_max = xmax.Max();

  if(out)
    *out << " \t" << xmin_min << " \t" << xmin_max << " \t" << xmax_min << " \t" << xmax_max;

  iteration->Get("neg_asym_min")->SetValue(xmin_min);
  iteration->Get("neg_asym_max")->SetValue(xmin_max);
  iteration->Get("pos_asym_min")->SetValue(xmax_min);
  iteration->Get("pos_asym_max")->SetValue(xmax_max);

  for(int i = 0; i < optimization->constraints.view->GetNumberOfActiveConstraints(); i++)
  {
    Condition* g = optimization->constraints.view->Get(i);
    // don't add for slope constraints
    if(!g->IsLocalCondition())
    {
      if(out)
        *out << " \t" << lamda[i];
      iteration->Get("lamba_" + g->ToString())->SetValue(lamda[i]);
    }
  }

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
      unsigned int m = optimization->constraints.view->GetNumberOfActiveConstraints();
      assert(si.lambda.GetSize() == m);
      assert(si.s.GetSize() == 2 * m);
      for(unsigned int ci = 0; ci < m; ci++)
      {
        string gn = optimization->constraints.view->Get(ci)->ToString();
        ipn->Get("lambda_" + gn)->SetValue(si.lambda[ci]);
        ipn->Get("mu_" + gn)->SetValue(si.mu[ci]);
        ipn->Get("s_grad_" + gn)->SetValue(si.s[ci]);
        ipn->Get("s_slack_" + gn)->SetValue(si.s[m + ci]);
      }
    }
  }
}

// Used for debugging, Helper function to read files used in InitilizeFromFile()
bool MMA::is_number(const std::string& s)
{
    try
    {
        std::stod(s);
    }
    catch(...)
    {
        return false;
    }
    return true;
}

/**
 * Used for debugging, will read the data from PETSC output file produced by topopt and inilitize.
 * */
void MMA::InitilizeFromFile(std::string filename, double *dp){
  std::ifstream input(filename);
    if (input) {
        std::string eachline;
        int i=0;
        while (std::getline(input, eachline)) {
            if(is_number(eachline))
            {
                std::istringstream iss(eachline);
                dp[i]=  std::stod(eachline);
                ++i;
            }
        }
  }
}

/** Used for debugging, content of the function modified based on the function being chceked
 * 1. We initilize all the data used by the function
 * 2. Call the function
 * 3. Output the data you want to check
 */
void MMA::StupidTest()
{

  //Set movie limits

  n=512;
  m=1;
  Matrix<double> temp(n,m);
  xval.Resize(n);  InitilizeFromFile("xval_DH.txt", xval.GetPointer());
  low.Resize(n);  InitilizeFromFile("L_DH.txt", low.GetPointer());
  upp.Resize(n);  InitilizeFromFile("U_DH.txt", upp.GetPointer());
  p_ij.Resize(temp); InitilizeFromFile("pij_DH.txt", p_ij[0]);
  q_ij.Resize(temp); InitilizeFromFile("qij_D.txt", q_ij[0]);
  alpha.Resize(n); InitilizeFromFile("alpha_GS.txt", alpha.GetPointer());
  beta.Resize(n); InitilizeFromFile("beta_GS.txt", beta.GetPointer());
  p_0j.Resize(n); InitilizeFromFile("p0_GS.txt", p_0j.GetPointer());
  q_0j.Resize(n); InitilizeFromFile("q0_GS.txt", q_0j.GetPointer());
  c[0] = 1000.0;
  lamda[0] = 500.000000;
  mu[0] = 1.0;
  a[0] = 0.0;

  HessianOfDual();

  LOG_DBG3(mmaTopOpt) << "StupidTest: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ " ;
  LOG_DBG3(mmaTopOpt) << "PD: dual_hessian[j]=" << dual_hessian[0];


  LOG_DBG3(mmaTopOpt) << "StupidTest: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ " ;
  assert (true && "Stupid test done");
}

//
//void MMA::StupidTest()
//{
//
//  //Set movie limits
//
//  n=512;
//  m=1;
//  k=4;
//  Matrix<double> temp(n,m);
//  asyminit=0.500000;
//  asymdec = 0.700000;
//  asyminc = 1.200000;
//  Xmin = 0.001000 ; Xmax = 1.000000 ; move_limits = 0.200000;
//  xval.Resize(n);  InitilizeFromFile("xval_GS.txt", xval.GetPointer());
//  xmin.Resize(n); InitilizeFromFile("xmin_GS.txt", xmin.GetPointer());
//  xmax.Resize(n); InitilizeFromFile("xmax_GS.txt", xmax.GetPointer());
//  low.Resize(n); InitilizeFromFile("L_GS.txt", low.GetPointer());
//  upp.Resize(n); InitilizeFromFile("U_GS.txt", upp.GetPointer());
//  xold1.Resize(n); InitilizeFromFile("xo1_GS.txt", xold1.GetPointer());
//  xold2.Resize(n); InitilizeFromFile("xo2_GS.txt", xold2.GetPointer());
//  alpha.Resize(n); InitilizeFromFile("alpha_GS.txt", alpha.GetPointer());
//  beta.Resize(n); InitilizeFromFile("beta_GS.txt", beta.GetPointer());
//  grad_func.Resize(n); InitilizeFromFile("dfdx_GS.txt", grad_func.GetPointer());
//  p_0j.Resize(n); InitilizeFromFile("p0_GS.txt", p_0j.GetPointer());
//  q_0j.Resize(n); InitilizeFromFile("q0_GS.txt", q_0j.GetPointer());
//  grad_con.Resize(n); InitilizeFromFile("dgdx_GS.txt", grad_con.GetPointer());
//  p_ij.Resize(temp); InitilizeFromFile("pij_GS.txt", p_ij[0]);
//  q_ij.Resize(temp); InitilizeFromFile("pij_GS.txt", q_ij[0]);
//
//
//
//
//  GenreteSubProblem();
//
//  LOG_DBG3(mmaTopOpt) << "StupidTest: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ " ;
//  LOG_DBG3(mmaTopOpt) << "GSP: xval=" << xval.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "GSP: low=" << low.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "GSP: upp=" << upp.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "GSP: xmin=" << xmin.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "GSP: xmax=" << xmax.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "GSP: alpha=" << alpha.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "GSP: beta=" << beta.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "GSP: grad_func=" << grad_func.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "GSP: po=" << p_0j.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "GSP: qo=" << q_0j.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "GSP: grad_con=" << grad_con.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "GSP: pij=" << p_ij.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "GSP: qij=" << q_ij.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "GSP: conA=" << conA[0];
//
//  LOG_DBG3(mmaTopOpt) << "StupidTest: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ " ;
//
//}

//
//void MMA::StupidTest()
//{
//  n=512;
//  m=1;
//
//  xval.Resize(n); InitilizeFromFile("xval.txt", xval.GetPointer());
//  xmin.Resize(n); InitilizeFromFile("xmin.txt", xmin.GetPointer());
//  xmax.Resize(n); InitilizeFromFile("xmax.txt", xmax.GetPointer());
//
//  Xmin = 0.001;
//  Xmax = 1.0;
//  grad_func.Resize(n,0.0); InitilizeFromFile("dfdx.txt", grad_func.GetPointer());
//
//  con.Resize(m,0.0); con[0] = -0.000000;
//
//  grad_con.Resize(n*m,0.0); InitilizeFromFile("dgdx.txt", grad_con.GetPointer());
//
//  low.Resize(n,0.0); //cdev important to set it to zero
//  upp.Resize(n,0.0); //cdev important to set it to zero
//  xold1 = xval;
//  xold2 = xval;
//  alpha.Resize(n);
//  beta.Resize(n);
//  change.Resize(n);
//  p_0j.Resize(n); // for the objective function
//  q_0j.Resize(n); // for the objective function
//  p_ij.Resize(m, n); // for all the constrain function
//  q_ij.Resize(m, n); // for all the constrain function
//  conA.Resize(m); // MMA approximation of constrain functions
//  lamda.Resize(m); // Lagrange multiplier
//  mu.Resize(m); // Lagrange multiplier
//  y.Resize(m); // elastic variable
//  z = 0.0; // elastic variable
//  a.Resize(m, 0.0); // cdev: terms in subproblem which i do not understand
//  c.Resize(m, 1000.0); // cdev: terms in subproblem which i do not understand
//  d.Resize(m, 0.0); // cdev: terms in subproblem which i do not understand
//  s.Resize(2*m, 0.0); // slack variables
//  dual_gradient.Resize(m, 0.0);
//  dual_hessian.Resize(m*m, 0.0);
//
//
//
//  LOG_DBG3(mmaTopOpt) << "StupidTest: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ " ;
//  LOG_DBG3(mmaTopOpt) << "SSP: xval=" << xval.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "SSP: xmin=" << xmin.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "SSP: xmax=" << xmax.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "SSP: grad_func=" << grad_func.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "SSP: grad_con=" << grad_con.ToString(0);
//  LOG_DBG3(mmaTopOpt) << "StupidTest: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ " ;
//
//
//  SolveMMA();
//}

