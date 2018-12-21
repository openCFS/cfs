#ifndef BFGS_HH_
#define BFGS_HH_

#include "General/Exception.hh"
#include "Optimization/Optimization.hh"
#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"
#include "Utils/StdVector.hh"

using std::pow;
using std::max;
using std::min;
using std::abs;
using std::string;

namespace CoupledField
{

class MMA;

class BFGS {

  public:
  // Constructor
  BFGS();
  BFGS(unsigned int n_, double tol_, unsigned int maxit_, unsigned int nsmax_, MMA* problem);

  void Initilize(unsigned int n_, double tol_, unsigned int maxit_, unsigned int nsmax_, MMA* problem);

  // Destructor
  ~BFGS();

  // BFGS Problem Solver
  void SolveBFGS(Vector<double>& x0_, Vector<double>& upp_, Vector<double>& low_);

  Vector<double> x; // solution


  // we have more BFGS_Info than subproblem
  struct BFGSInfo {
    Vector<double> xc; // arg values
    double fc; // function values
    double norm_pgc; // norm of gradients
    Vector<double> pgc; // arg values
    unsigned int iter; // for the sub probelm IP interation count
    unsigned int nactive; // number of active constraints
  };

  StdVector<BFGSInfo> bfgs_details;

  /** @see BaseOptimier */
  virtual void LogFileLine(std::ofstream* out, PtrParamNode iteration);
  //  void LogFileHeader(Optimization::Log& log);

  private:



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
  unsigned int n = 0; // dimension of the design variable

  Vector<double> x0; // initial design


  StdVector<double> upp, low; // upper and lower bound of the design variable
  double tol = 1.0e-6; // = termination criterion norm(grad) < tol optional, default = 1.d-6
  unsigned int maxit = 1000; // maximum iterations
  unsigned int nsmax=maxit;
}; // end of class
} // end of namespace

#endif /* BFGS_HH_ */
