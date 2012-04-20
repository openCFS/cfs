#ifndef FEASPP_HH_
#define FEASPP_HH_


#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "Optimization/Optimizer/FeasSubProblem.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
using namespace boost::numeric::ublas;

namespace CoupledField
{

class Approximation;
class FeasPP;
template <class TYPE> class SCRS_Matrix;

/** this is a C++ implementation of MMA using IPOPT as subproblem solver. A feature is
 * to forward feasibility constraints to the subproblem, see thesis of Sonja Lehmann.
 * The only dependency is IPOPT and we use the dependencies of the IPOPT CFS solver for
 * the FeasPP subsolver IPOPT. Hence USE_IPOPT needs to be enabled to build FeasPP.
 */

class FeasPP : public BaseOptimizer
{
public:
  FeasPP(Optimization* optimization, PtrParamNode pn);

  virtual ~FeasPP();

  /** initialize problem
   * @see BaseOptimizer::PostInit() */
  void PostInit();

  /** Construct the global Hessian from the functions */
  void SetupHessian();

  /** number of total design variables */
  unsigned int n;

  /** number of constraint functions */
  unsigned int m;

  Approximation* obj; // created in PostInit()

  StdVector<Approximation*> constr;

  /** the cfs design variable */
  StdVector<double> x_outer;

  /** lower asymptote */
  StdVector<double> L;

  /** upper asymptote */
  StdVector<double> U;

  /** are there constraints which are not to be approximated?
      The objective can be queried by its own flag.  */
  bool non_approx_constraints;

  bool approx_linear;

  bool approx_feasibility;

  /** In the feasibility paper in (6) is a convexication term for the objective only.
   * The parameter is a scalar, fixed tau. Setting to 0.0 disables the term! */
  double convex_tau;

  compressed_matrix<double>* hessian;

protected:

  void SolveProblem();

private:

  /** setup and solve the subproblem by ipopt */
  void SolveSubProblem();

  /** updates the design and the outer function values and gradients */
  void UpdateToCurrentStep();

  typedef struct
  {
    bool old_point_is_optimal;
    int steps;
    double org_dx;
    double curr_dx;
  } BTI;

  /** performs backtracking linesearch. obj(x_new) >= obj(x_old)= ob->outer_value
   * Note that the objective is evaluated, the design_id is set but the constraints
   * might not be evaluated -> UpdateSystem().
   * The final design is stored in the CFS-Design!
   * @param the solution of the subproblem.
   * @return number of function evaluations and the norms of the designs*/
  BTI Backtracking(const Vector<double>& x_old, const Vector<double>& d, const StdVector<double>& x_new);

  typedef enum { NONE, BACKTRACKING } Globalization;

  static Enum<Globalization> global;

  Globalization global_;

  SmartPtr<FeasSubProblem> ipopt;
};

/** this is either the approximation of a function or the function itself */
class Approximation
{
public:
  /** @param constraint_idx -1 for objective */
  Approximation(FeasPP* feas_pp, int constraint_idx, bool approximate);

  typedef enum { FUNC, GRAD, HESSIAN } Eval;

  /** evaluate function value or first derivative according to the MMA approximation */
  double Evaluate(const double* x_inner, Eval eval, StdVector<double>* out = NULL);

  void AddHessianPattern(compressed_matrix<double>& hessian);

  /** helper for logging */
  std::string ToString();

  /** the gradient of the real (outer) function) */
  StdVector<double> outer_grad;

  /** this stores the sparsity pattern of the Jacobian, might be dense  */
  StdVector<unsigned int> jac_pattern;

  /** this stores the sparsity pattern of the Hessian. Is diagonal for the approximations.
   * The it has 2 columns and * rows with first column i and second j. For diagonal i = j.
   * Stores only the upper triagonal as the Hessian needs to be symmetric */
  Matrix<unsigned int> hess_pattern;

  /** bounds */
  double lower;
  double upper;

  /** function value of the outer function */
  double outer_val;

  /** shall this function be approximated */
  bool approximate;
private:

  /** Remember do reset the condition container via Done()!
   * @see constraint_idx. */
  Function* GetFunction();

  Condition* GetCondition();

  double EvalApproximation(const double* x_inner, Eval eval, StdVector<double>* out);

  double EvalDirect(const double* x_inner, Eval eval, StdVector<double>* out);

  /** We cannot store a function pointer because of the local constraints where we need an index.
   * -1 is for objective */
  int constraint_idx;

  FeasPP* common;
};



} // end of namespace


#endif /* FEASPP_HH_ */
