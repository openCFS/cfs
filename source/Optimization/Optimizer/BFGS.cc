#include "Optimization/Optimizer/BFGS.hh"
#include "Optimization/Optimizer/MMA.hh"
#include "Utils/ToolsFull.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ProgramOptions.hh"

DEFINE_LOG(bfgs_LOG, "bfgs_log")

using namespace CoupledField;
using std::max;
using std::min;

BFGS::BFGS(unsigned int n, double tol, unsigned int maxit, unsigned int nsmax, MMA *problem)
{
  this->prob= problem;
  this->n=n;
  this->tol=tol;
  this->maxit=maxit;
  this->nsmax=nsmax;

  x0.Resize(n);
  low.Resize(n);
  upp.Resize(n);
}


int BFGS::SolveBFGS(Vector<double>& xc, Vector<double>& up, Vector<double>& lo)
{
  bfgs_details.Resize(0); // prepare for logging

  LOG_DBG3(bfgs_LOG) << "SolB: initial xc=" << xc.ToString();
  LOG_DBG3(bfgs_LOG) << "SolB: initial up=" << up.ToString();
  LOG_DBG3(bfgs_LOG) << "SolB: initial lo=" << lo.ToString();

  /** Sanity check for the dimension */
  assert(n== xc.GetSize() && n== up.GetSize() && n== lo.GetSize());
  for(unsigned int in =0; in<n; ++in)
    assert(lo[in] < up[in]);

  // Project x inside the bounds and make it feasible
  xc = kk_proj(xc, up, lo);

  double alp=1.0e-4;
  /** For limited memory hessian computation we store a pairs {sk, yk}
   * nsmax determines how many points we store*/

  Matrix<double> ystore(n,nsmax+1);
  Matrix<double> sstore=ystore;
  int ns= -1; //TODO: I have not understood what this is? This value is taken from TimKelly code

  // Start of iteration
  unsigned int itc = 0;
  double fc;
  Vector<double> gc(n);
  fc = prob->EvalDualFunction(xc);
  gc = prob->EvalDualGrad(xc);

  unsigned int iarm=0;

  /** To store iteration history
   * [norm(projgrad), f, number of step length reductions,
   * iteration count,relative size of active set]
   */
  // Matrix<double> ithist(maxit,5);

  Vector<double> x_gc(n);
  x_gc = xc - gc;

  Vector<double> pgc(n);
  pgc = xc - kk_proj(x_gc, up, lo);

  // active index set
  Vector<int> alist(n);

  Vector<double> tst(n);
  tst = up-lo;

  double lim1 = 0.5*tst.Min();
  double epsilon = min(lim1, pgc.NormL2());

  /** interestingly ia is never used - might indicate a bug! TOOD
  // check if the x is on the bounds, if it is on the bounds then increase the number of active index
  unsigned int ia=0;
  for(unsigned int in =0; in < n; ++in)
  {
    if (xc[in] == up[in] || xc[in] == lo[in])
      ++ia;
  } */

  for(unsigned int in =0; in < n; ++in)
  {
    if (min(up[in] - xc[in], xc[in] - lo[in]) < epsilon)
      alist[in] = 1;
  }

  // Optimization loop begins
  while(pgc.NormL2() > tol && itc <= maxit)
  {
    double lam=1;

    Vector<double> dsd(n);
    dsd = -gc;

    dsd = bfgsrp(ystore, sstore, ns, dsd, alist);

    Vector<double> neg_gc = gc;
    neg_gc.ScalarMult(-1.0);

    dsd += proja(neg_gc, alist);

    Vector<double> xc_lam_dsd = dsd ;
    xc_lam_dsd.ScalarMult(lam);
    xc_lam_dsd += xc;
    Vector<double> xt = kk_proj(xc_lam_dsd, up , lo);

    double ft=0.0;
    ft = prob->EvalDualFunction(xt);

    ++itc;
    iarm =0;

    // Vector<double> pl = xc - xt;
    Vector<double> pl(n);
    pl = xc - xt;

    /** fgoal=fc-(gc'*pl)*alp;
     * First I compute the inner product (gc'*pl)
     * and then compute the final value.
     */
    double fgoal = gc.Inner(pl);
    fgoal = fgoal*alp;
    fgoal = fc - fgoal;

    /** Simple line search */
    while (ft > fgoal)
    {
      lam = lam * 0.1;
      ++iarm;
      Vector<double> xta = xt;
      Vector<double> xc_lam_dsd = dsd ;
      xc_lam_dsd.ScalarMult(lam);
      xc_lam_dsd += xc;
      xt = kk_proj(xc_lam_dsd, up, lo);
      pl = xc -xt;
      ft = prob->EvalDualFunction(xt);

      if(iarm > 10)
      {
        x =xc;
        (void) std::setprecision(9); // cast to get rid of ignoring return value of function declared with 'nodiscard' attribute
        if(progOpts->DoDetailedInfo())
        {
          bfgs_details.Push_back(BFGSInfo());
          bfgs_details.Last().fc = fc;
          bfgs_details.Last().iter = itc;
          bfgs_details.Last().norm_pgc = pgc.NormL2();
          bfgs_details.Last().pgc = pgc;
          bfgs_details.Last().xc = xc;
        }
        return itc;
      }
      fgoal = gc.Inner(pl);
      fgoal = fgoal*alp;
      fgoal = fc - fgoal;

    }
    fc = prob->EvalDualFunction(xt);

    Vector<double> gp;

    gp = prob->EvalDualGrad(xt);

    Vector<double> y(n);
    y = gp - gc;

    Vector<double> s(n);
    s= xt - xc;

    gc = gp;
    xc = xt;
    Vector<double> xc_gc(n);
    xc_gc = xc - gc;

    pgc = xc - kk_proj(xc_gc, up, lo);
    epsilon = min(lim1, pgc.NormL2());

    // Resetting alist to zero.
    for (unsigned int in =0; in <n ; ++in)
      alist[in] = 0;

    for(unsigned int in =0; in<n; ++in)
    {
      if(min(up[in] - xc[in], xc[in] - lo[in]) < epsilon)
        alist[in] = 1;
    }
    y = proji(y,alist);
    s = proji(s,alist);

    /** restart if y'*s is not positive or we're out of room */
    double yts = y.Inner(s);
    if (yts <= 0)
      ns= -1;

    if(ns == (int)nsmax)
      ns = -1;
    else if(yts > 0)
    {
      ++ns;
      double alpha = sqrt(yts);
      for(unsigned int in =0; in < n; ++in)
      {
        sstore(in, ns) = s[in]/alpha;
        ystore(in, ns) = y[in]/alpha;
      }
    }

    int ia = 0;
    for(unsigned int in =0; in<n; ++in)
    {
      if(xc[in] == up[in] || xc[in] == lo[in])
      {
        ++ ia;
      }
    }

    if(progOpts->DoDetailedInfo())
    {
      bfgs_details.Push_back(BFGSInfo());
      bfgs_details.Last().nactive = ia;
      bfgs_details.Last().fc = fc;
      bfgs_details.Last().iter = itc;
      bfgs_details.Last().norm_pgc = pgc.NormL2();
      bfgs_details.Last().pgc = pgc;
      bfgs_details.Last().xc = xc;
    }

  } // end of iteration loop
  LOG_DBG3(bfgs_LOG) <<"SolB: Sub prob iterations: " <<itc <<" Sol= "<<xc.ToString();

  x = xc;
  return itc;
}


/** Projection onto the active set
 * If x < lo; then x = low
 * if x > up; then x = upp*/
Vector<double> BFGS::kk_proj(const Vector<double> &x, const Vector<double> &up, const Vector<double> &lo)
{
  Vector<double> proj(n);
  for(unsigned int in=0; in < n; ++in)
  {
    proj[in] = min(up[in], x[in]);
    proj[in] = max(x[in], lo[in]);
  }
  return proj;
}

/** projection onto epsilon-active set */
Vector<double> BFGS::proja(const Vector<double>& x, const Vector<int>& alist)
{
  Vector<double> px = x;
  for(unsigned int in=0; in<n; ++in)
    if(alist[in]==0)
      px[in] = 0;

  return px;
}

/** projection onto epsilon-inactive set */
Vector<double> BFGS::proji(const Vector<double>& x, const Vector<int>& alist)
{
  Vector<double> px = x;
  for(unsigned int in=0; in<n; ++in)
    if(alist[in]==1)
      px[in] = 0;

  return px;
}

Vector<double> BFGS::bfgsrp( Matrix<double> &ystore,  Matrix<double> &sstore, const int ns,
    const Vector<double> &dsd, const Vector<int> &alist)
{
  Vector<double> dnewt = proji(dsd, alist);
  if(ns == -1)
    return dnewt;

  Vector<double> sstore_ns(n);
  sstore.GetCol(sstore_ns, ns);
  Vector<double> sstore_ns_proji = proji(sstore_ns, alist);

  Vector<double> ystore_ns(n);
  ystore.GetCol(ystore_ns, ns);
  Vector<double> ystore_ns_proji = proji(ystore_ns, alist);

  for(unsigned int in =0; in < n; ++in)
  {
    assert(ns >= 0);
    sstore(in, ns) = sstore_ns_proji[in];
    ystore(in, ns) = ystore_ns_proji[in];
  }

  sstore.GetCol(sstore_ns, ns);
  double beta = sstore_ns.Inner(dsd);

  ystore.GetCol(ystore_ns, ns);
  ystore_ns.ScalarMult(beta);
  dnewt = dsd - ystore_ns;

  Vector<int>xlist(n);
  dnewt=bfgsrp(ystore,sstore,ns-1,dnewt,xlist);

  ystore.GetCol(ystore_ns, ns);

  sstore.GetCol(sstore_ns, ns);
  double b_y_d = beta-ystore_ns.Inner(dnewt);
  sstore_ns.ScalarMult(b_y_d);
  dnewt=dnewt+sstore_ns;
  dnewt = proji(dnewt, alist);

  return dnewt;
}

void BFGS::LogFileLine(PtrParamNode iteration)
{

  if(progOpts->DoDetailedInfo())
  {
    PtrParamNode sub =iteration->Get("bfgs_subproblem");
    for(unsigned int ib=0; ib < bfgs_details.GetSize(); ++ib)
    {
      BFGSInfo& bi = bfgs_details[ib];
      PtrParamNode ipn = sub->Get("iter", ParamNode::APPEND);
      ipn->Get("number")->SetValue(ib);
      ipn->Get("iter")->SetValue(bi.iter);
      ipn->Get("current_x")->SetValue(bi.xc);
      ipn->Get("current_func")->SetValue(bi.fc);
      ipn->Get("norm_grad")->SetValue(bi.norm_pgc);
      ipn->Get("active_con")->SetValue(bi.nactive);
      ipn->Get("grad")->SetValue(bi.pgc);
    }
  }
}



