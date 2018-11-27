#include "Optimization/Optimizer/BFGS.hh"


#include "Utils/tools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ProgramOptions.hh"

DECLARE_LOG(bfgs_LOG)
DEFINE_LOG(bfgs_LOG, "bfgs_log")


using namespace CoupledField;
using std::pow;
using std::max;
using std::min;
using std::abs;
using std::string;

//TestFunction::TestFunction(unsigned int n_, unsigned int m_):n(n_),m(m_){}
//unsigned int GetNoConstraints(){return m;}
/** Function is defined to e^(x^2)
 * So the gradint is going to be 2x*e^(x^2)
 */
double TestFunction::EvalDualFucntion(Vector<double> xin){
  assert(xin.GetSize() == n);
  double fout =0;
//  for (unsigned int in=0; in <n; ++in )
//    fout += exp(pow(xin[in], 2));
  for (unsigned int in=0; in <n; ++in )
  {
    fout += pow(xin[in], 2);
    if(in%2)
      fout += std::cos(xin[in]);
    else
      fout += std::sin(xin[in]);

  }
  return fout;
}

Vector<double> TestFunction::EvalDualGrads(Vector<double> xin){
  assert(xin.GetSize() == n);
  Vector<double> gout(n);
  for (unsigned int in=0; in <n; ++in)
  {
    gout[in] = 2.0*xin[in];
    if(in%2)
      gout[in] -= std::sin(xin[in]);
    else
      gout[in] += std::cos(xin[in]);

  }
  return gout;
}

BFGS::BFGS(unsigned int n_, double tol_, unsigned int maxit_, unsigned int nsmax_, MMA *problem)
          :prob(problem), n(n_), tol(tol_), maxit(maxit_), nsmax(nsmax_)
{
  x0.Resize(n);
  low.Resize(n);
  upp.Resize(n);
  prob_t = NULL;

}

BFGS::BFGS(unsigned int n_, double tol_, unsigned int maxit_, unsigned int nsmax_,  TestFunction *problem)
          :prob_t(problem), n(n_), tol(tol_), maxit(maxit_), nsmax(nsmax_)
{
  x0.Resize(n);
  low.Resize(n);
  upp.Resize(n);
  prob = NULL;

}

BFGS::~BFGS() {}

void BFGS::SolveBFGS(Vector<double>& xc, Vector<double>& up, Vector<double>& lo)
{

  LOG_DBG3(bfgs_LOG) << "SolB: initial xc=" << xc.ToString(0);
  LOG_DBG3(bfgs_LOG) << "SolB: initial up=" << up.ToString(0);
  LOG_DBG3(bfgs_LOG) << "SolB: initial lo=" << lo.ToString(0);

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
  Matrix<double>sstore=ystore;
  int ns= -1; //TODO: I have not understood what this is? This value is taken from TimKelly code

  // Start of iteration
  unsigned int itc = 0;
  double fc;
  Vector<double> gc;
  if(prob_t == NULL)
  {
    fc = prob->EvalDualFucntion(xc);
    gc = prob->EvalDualGrads(xc);
  }
  else
  {
    fc = prob_t->EvalDualFucntion(xc);
    gc = prob_t->EvalDualGrads(xc);
  }

  unsigned int iarm=0;

  /** To store iteration history
   * [norm(projgrad), f, number of step length reductions,
   * iteration count,relative size of active set]
   */
  // Matrix<double> ithist(maxit,5);

  Vector<double> x_gc = xc - gc; // TODO: this is comnputed only because the kk_proj did not accept x-gc as a argument
  Vector<double> pgc = xc - kk_proj(x_gc, up, lo);

  // active index set
  unsigned int ia=0;
  Vector<int> alist(n);
  Vector<double> tst = up-lo;
  double lim1 = 0.5*tst.Min();
  double epsilon = min(lim1, pgc.NormL2());

  // check if the x is on the bounds, if it is on the bounds then increase the number of active index
  for(unsigned int in =0; in < n; ++in)
  {
    if (xc[in] == up[in] || xc[in] == lo[in])
      ++ia;
  }

  for(unsigned int in =0; in < n; ++in)
  {
    if (min(up[in] - xc[in], xc[in] - lo[in]) < epsilon)
      alist[in] = 1;
  }

  // Optimization loop begins
  while(pgc.NormL2() > tol && itc <= maxit)
  {
    double lam=1;
    Vector<double> dsd = -gc;
    dsd = bfgsrp(ystore, sstore, ns, dsd, alist);

    Vector<double> neg_gc = gc;
    neg_gc.ScalarMult(-1.0);

    dsd += proja(neg_gc, alist);

    Vector<double> xc_lam_dsd = dsd ;
    xc_lam_dsd.ScalarMult(lam);
    xc_lam_dsd += xc;
    Vector<double> xt = kk_proj(xc_lam_dsd, up , lo);

    double ft=0.0;
    if(prob_t == NULL)
      ft = prob->EvalDualFucntion(xt);
    else
      ft = prob_t->EvalDualFucntion(xt);

    ++itc;
    iarm =0;

    Vector<double> pl = xc - xt;

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
      if(prob_t == NULL)
        ft = prob->EvalDualFucntion(xt);
      else
        ft = prob_t->EvalDualFucntion(xt);

      if(iarm > 10)
      {
        x =xc;
        return;
      }
      fgoal = gc.Inner(pl);
      fgoal = fgoal*alp;
      fgoal = fc - fgoal;

    }
    if(prob_t == NULL)
      fc = prob->EvalDualFucntion(xt);
    else
      fc = prob_t->EvalDualFucntion(xt);

    Vector<double> gp;

    if(prob_t == NULL)
      gp = prob->EvalDualGrads(xt);
    else
      gp = prob_t->EvalDualGrads(xt);

    Vector<double> y = gp - gc;
    Vector<double> s = xt - xc;
    gc = gp;
    xc = xt;
    Vector<double> xc_gc = xc - gc;
    pgc = xc - kk_proj(xc_gc, up, lo);
    epsilon = min(lim1, pgc.NormL2());

    // Resetting alist to zero.
    for (unsigned int in =0; in <n ; ++in)
      alist[in] = 0;

    int ial = 0;
    for(unsigned int in =0; in<n; ++in)
    {
      if(min(up[in] - xc[in], xc[in] - lo[in]) < epsilon)
      {
        alist[in] = 1;
        ++ ial;
      }
    }
    y = proji(y,alist);
    s = proji(s,alist);

    /** restart if y'*s is not positive or we're out of room */
    double yts = y.Inner(s);
    if (yts <= 0)
      ns= -1;

    if(ns == nsmax)
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
    LOG_DBG3(bfgs_LOG) << "SolB: at iteration : "<<itc<<" xc= " << xc.ToString(0)<<" fc= " << fc<<" fgoal= " << fgoal;
  } // end of iteration loop
  LOG_DBG3(bfgs_LOG) <<"SolB: Sub prob iterations: " <<itc <<" Sol= "<<xc.ToString();
  x = xc;
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


