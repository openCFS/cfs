#ifndef BFGS_HH_
#define BFGS_HH_

#include "General/Exception.hh"

#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"
#include "Utils/StdVector.hh"
#include "Optimization/Optimizer/MMA.hh"


using namespace CoupledField;
using std::pow;
using std::max;
using std::min;
using std::abs;
using std::string;


class TestFunction {
  public:
    TestFunction(unsigned int n_):n(n_){}
    unsigned int GetNoDesign(){return n;}
    double EvalDualFucntion(Vector<double> xin);
    Vector<double> EvalDualGrads(Vector<double> xin);

  private:
    unsigned int n;

};



class BFGS {

  public:
  // Constructor
  BFGS(unsigned int n_, double tol_, unsigned int maxit_, unsigned int nsmax_, MMA *problem);
  BFGS(unsigned int n_, double tol_, unsigned int maxit_, unsigned int nsmax_, TestFunction *problem);

  // Destructor
  ~BFGS();

  // BFGS Problem Solver
  void SolveBFGS(Vector<double>& x0_, Vector<double>& upp_, Vector<double>& low_);

  Vector<double> x; // solution


//  // we have more BFGS_Info than subproblem
//  struct SubInfo {
//    StdVector<double> xval; // Lagrange multipler for all constraints
//    StdVector<double> mu;
//    StdVector<double> s; // slack variables for all contraints + 1 for IP slack
//    double err;
//    double epsi; // for the fixed number of subproblems to be solved
//    int iter; // for the sub probelm IP interation count
//  };

  private:

  /** @see BaseOptimier */
  void LogFileHeader(Optimization::Log& log);
  void LogFileLine(std::ofstream* out, PtrParamNode iteration);

  /** Projection onto the active set
   * If x < lo; then x = low
   * if x > up; then x = upp*/
  Vector<double> kk_proj(const Vector<double> &x, const Vector<double> &up, const Vector<double> &lo);


  Vector<double> bfgsrp( Matrix<double> &ystore,  Matrix<double> &sstore,
          const int ns, const Vector<double> &dsd, const Vector<int> &alist);

  /** projection onto epsilon-active set */
  Vector<double> proja(const Vector<double>& x, const Vector<int>& alist);

  /** projection onto epsilon-inactive set */
  Vector<double> proji(const Vector<double>& x, const Vector<int>& alist);

  MMA *prob;
  TestFunction *prob_t;
  unsigned int n = 0; // dimension of the design variable

  Vector<double> x0; // initial design


  StdVector<double> upp, low; // upper and lower bound of the design variable
  double tol = 1.0e-6; // = termination criterion norm(grad) < tol optional, default = 1.d-6
  unsigned int maxit = 1000; // maximum iterations
  unsigned int nsmax=maxit;


};

#endif /* BFGS_HH_ */
