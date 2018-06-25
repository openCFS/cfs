#ifndef MMA_HH_
#define MMA_HH_


#include "Optimization/Optimizer/BaseOptimizer.hh"
//#include "Optimization/Optimizer/MMASubProblem.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "MatVec/Vector.hh"

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace boost::numeric::ublas;

namespace CoupledField
{
class MMA;
template <class TYPE> class SCRS_Matrix;

class MMA : public BaseOptimizer
{
public:
  MMA(Optimization* optimization, PtrParamNode pn);

  virtual ~MMA();

protected:
  void PostInit();

  /** @see BaseOptimizer::ToInfo() */
  void ToInfo(PtrParamNode pn);

  /** Optimization Problem
   *      min func(xval)
   *   s.t:   con_i <= 0 ; i=[1 , m]
   *          xmin_i <= xval_i <= xmax_i ; i =[1, n]
   */

  // physics evaluation
  unsigned int n =0; // n number of design variables
  unsigned int m = 0; // m n number of  contrains

  // Design variable, min bound and max bound
  StdVector<double> xval;
  /** this is what is called the asymptotes in MMA and low/upp in Svanberg's mmasub.m.\
   * The std_xmin case holds the data for WriteBoundsToExtern. */
  StdVector<double> std_xmin;
  StdVector<double> std_xmax;

  /** is a math object to std_* with no own memory
   * @see Vector::Replace */
  Vector<double> xmin;
  Vector<double> xmax;

  double compliance = -1.0; // compliance (objective)
  StdVector<double> grad_compliance; // gradient of compliance function, dimension n
  StdVector<double> constrains; // vector with contrain values, dimension m
  StdVector<double> grad_constrains; // gradient of contrains, dimension m x n.
  double obj_scale = 0; //scaling factor for objective ToDO: Idea behind choosing value for scaling is not clear.

  /** MMA approximation
   *      p_ij = (upp_j - xval_j)^2 * (d.f_i/d.xval_j)      if (d.f_i/d.xval_j) > 0
   *          = 0                                           if (d.f_i/d.xval_j) <= 0
   *
   *      q_ij = 0                                          if (d.f_i/d.xval_j) >= 0
   *          = - (xval_j - low_j)^2 * (d.f_i/d.xval_j)     if  (d.f_i/d.xval_j) < 0
   *
   *      r_i = f_i(xavl) - summation(p_ij/(upp_j - xval_j) + q_ij/(xval_j - low_j))
   *
   *      funcA_i = r_i + summation(p_ij/(upp_j - xval_j) + q_ij/(xval_j - low_j))
   */

  StdVector<double> p_0j, q_0j; // objective approx
  Matrix<double> p_ij, q_ij; // constraint approx
  StdVector<double> b; // rhs of constrain inequality in subproblem.

  /** MMA Sub Problem
   *      min funcA(xval) + z + 1/2 z^2 + summation(y_i*c_i + 1/2 y_i^2 ; i = [1, m])
   *   s.t:   consA_i(xval) - a_i*z - y_i <=0 ; i = [1, m]
   *          alf_i <= xval <= bet_i ; i =[1, n]
   *          y_i >= 0
   *          z >= 0
   */

  // MMA mathematics parameters
  double  asyminit = 0.9;
  double asymdec = 0.5;
  double asyminc = 1.2;
  double move_limits = 0.2; // asymptotes speed control

  /** for fixed asymptotes only */
  double asym_fixed_lower = -1;
  double asym_fixed_upper = -1;

  unsigned int max_sub_iter = 10; // maximum iteration for subproblem solver
  StdVector<double> a, c, d; // penalty parameter for the subproblem
  StdVector<double> alpha, beta; // sub problem move limits, dimension n
  StdVector<double> y; //elastic variable to make the problem always fesible
  double z=0.0; // elastic variable
  bool robustAsymptotes = false; // based on the implementation in topopt petsc code, check in MMA.cc in function GenSub
  bool fixedAsymptotes = false;

  /** Globally convergent version refer svanberg DCMMA section 6 */
  bool globallyConvergent = false;
  StdVector<double> rho; // will be useful only when globallyConvergent is required


  // Other help variables
  StdVector<double> lamda, mu, s; // Lagrange multipliers and slack variables
  StdVector<double> low, upp; // Asymptotes bound
  StdVector<double> xold1, xold2; // Save old design vaiables to tune asymptotes
  StdVector<double> change; // Design change
  double tol = 0.0001;

  StdVector<double> dual_gradient, dual_hessian;

  // we have more SubInfo than subproblem iteration as the iterations of the interior point method are also output
  struct SubInfo {
    StdVector<double> lambda; // Lagrange multipler for all constraints
    StdVector<double> mu;
    StdVector<double> s; // slack variables for all contraints + 1 for IP slack
    double err;
    double epsi; // for the fixed number of subproblems to be solved
    int iter; // for the sub probelm IP interation count
  };

  StdVector<SubInfo> subiters;

protected:

  void SolveProblem();

  void AdjustMoveLimits();

  bool SolveMMA();

  void GenreteSubProblem();

  void ComputeObjectiveConstraintsSensitivities();

  /** @see BaseOptimizer::ToInfo() */
//  void ToInfo(PtrParamNode pn);

  /** @see BaseOptimier */
  void LogFileHeader(Optimization::Log& log);
  void LogFileLine(std::ofstream* out, PtrParamNode iteration);

  void StupidTest();
  void InitilizeFromFile(std::string filename, double * vec);
  bool is_number(const std::string& s);

private:

  bool SolveSubProblem();

  void PrimalVarFromDualVar();

  void GradientOfDual();

  void HessianOfDual();

  void Factorize(StdVector<double> & , const unsigned int );

  void Solve(StdVector<double> & , StdVector<double> &, const int );

  void DualLineSearch();

  double DualResidual(double);

  /** generate sub problem */
  Timer* gsp_timer_ = NULL;

  /** generate sub problem */
  Timer* sps_timer_ = NULL;

  /** when the subproblem failes, the error message set */
  std::string mma_error_;
};


} // end of namespace


#endif /* MMA_HH_ */
