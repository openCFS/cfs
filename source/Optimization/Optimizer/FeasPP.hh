#ifndef FEASPP_HH_
#define FEASPP_HH_


#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "Optimization/Optimizer/FeasSubProblem.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "MatVec/vector.hh"

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

  /** number of constraint functions
   * @see m_c and m_e */
  unsigned int m;

  Approximation* obj; // created in PostInit()

  StdVector<Approximation*> constr;

  /** the cfs design variable */
  Vector<double> x_outer;

  /** lower asymptote */
  StdVector<double> L;

  /** upper asymptote */
  StdVector<double> U;

  /** As we might use design bounds, this are the lower design bounbs. Either from design or the design bound constraints
   * Indices by design type */
  StdVector<double> lower_bound;
  StdVector<double> upper_bound;

  /** are there constraints which are not to be approximated?
      The objective can be queried by its own flag.  */
  bool non_approx_constraints;

  bool approx_linear;

  bool approx_feasibility;

  /** In the feasibility paper in (6) is a convexication term for the objective only.
   * The parameter is a scalar, fixed tau. Setting to 0.0 disables the term! */
  double convex_tau;

  compressed_matrix<double>* hessian;

  /** are the design bounds for the sub-problem the ones from the original problem or dynamic, depending on L/U and x */
  bool dynamic_design_bounds;

  /** do we approximate benson vanderbei constraints in the subproblem by determinant constraints? Then also the lagrange multiplyer
   * from the subproblem will be transformed. See the feasibility paper! */
  bool approx_vanderbei_by_determinants;

protected:

  void SolveProblem();

  /** @see BaseOptimizer::ToInfo() */
  void ToInfo(PtrParamNode pn);

private:

  /** setup and solve the subproblem by ipopt */
  void SolveSubProblem();

  /** updates the design and the outer function values and gradients */
  void UpdateToCurrentStep();

  /** calculates the augmented Lagrangian for globalization
   * Feasibility paper (16)
   * @param x the functions will be evaluated at x
   * @param rho vector of penalty parameters, first the 'normal' constraints, then the feasibility constraints. */
  double CalcAugmentedLagrangian(const Vector<double>& x, const Vector<double>& y, const Vector<double>& rho);

  /** the gradient with respect to x and y written all to grad
   * @param grad size n + m */
  void CalcGradAugmentedLagrangian(const Vector<double>& x, const Vector<double>& y, const Vector<double>& rho, Vector<double>& grad);

  /** strictly feasibility papert (18) and (19). necessary to ensure convergence */
  double CalcEta(double tau, const Vector<double>& x_vec, const Vector<double>& z_vec);

  /** calculates the penality parameter sfor the augmented Lagrangian.
   * Feasibility paper (20) and (21)
   * @param diff = norm(z-x)^2
   * @param rho old and to be overwritten */
  void CalcPenaltyRho(double eta, double diff, const Vector<double>& y_vec,  const Vector<double>& v_vec, Vector<double>& rho) const;

  typedef struct
  {
    bool old_point_is_optimal;
    int steps;
    double stepwidth;
    double org_dx;
    double curr_dx;
  } LSR; // line search result

  /** performs backtracking linesearch. obj(x_new) >= obj(x_old)= ob->outer_value
   * Note that the objective is evaluated, the design_id is set but the constraints
   * might not be evaluated -> UpdateSystem().
   * The final design is stored in the CFS-Design!
   * @param the solution of the subproblem.
   * @return number of function evaluations and the norms of the designs*/
  LSR Backtracking(const Vector<double>& x_old, const Vector<double>& x_new);

  /**
   * @param k the iteration. 0 for the very first call after leaving initial design
   * @param v the new lambda which gets updated by linesearch! */
  LSR AugmentedLagrangianLineSearch(int k, const Vector<double>& x, const Vector<double>& z, const Vector<double>& y, Vector<double>& v, PtrParamNode in);

  /** update the asymptotes as long as they are not set to fixed!
   * @param force_reduction to react on subproblem problems */
  void UpdateAsymptotes(const Vector<double>&x_outer, int iter, bool force_reduction = false);

  /** calc a not KKT like stopping criteria following (7.24) in Sonja's thesis */
  // double CalcStopingCriteria(const Vector<double>& x, const Vector<double>& x_old, const Vector<double>& y);

  /** calc the three KKT values and outputs them to the param node
   * @approx true for inner problem (ipopt-result), false for outer-problem */
  void CalcKKT(const Vector<double>& x, PtrParamNode in,  bool sub, bool det, bool max_norm = true);

  double CalcAngle();

  /** assume the current design to be FMO tensors and output them */
  void DumpFMPTensors();

  /** translates between determinant constraints and benson vanderbei constraints */
  Function::Type TranslateFeasibilityConstraint(Function::Type type) const;

  typedef enum { NONE, BACKTRACKING, AUG_LAGRANGIAN } Globalization;

  Enum<Globalization> global;

  Globalization global_;

  typedef enum { FIXED, MMA } Asymptotes;

  Enum<Asymptotes> asymptotes;

  Asymptotes asymptotes_;

  SmartPtr<FeasSubProblem> ipopt;

  /** number of non-feasibility constraints. m = m_c + me */
  unsigned int m_c;

  /** number of feasibility constraints */
  unsigned int m_e;

  /** This is the set of penalty parameters rho for the augmented Lagrangian */
  Vector<double> rho;

  /** the initial value for rho, optionally to be set in xml */
  double rho_init_;

  /** constant eta in (22) - different from eta_k (18) ... (21) */
  double rho_eta_;

  /** decrease factor for augmented Lagrangian Armijo rule: my in (23) in (0,1) */
  double decrease_;

  /** Augmented Lagrangian: Stepwidth for Armijo rule: beta in Step 5 in Algorithm 1 in (0,1) */
  double stepwidth_;

  /** the minimal step width factor form augmented Lagrangian and backtracking */
  double min_step_;

  /** MMA parameters for asymptotes update */
  double mma_shrink_;
  double mma_grow_;
  double mma_dist_;

  /** this design history is used and controlled by UpdateAsymptotes() */
  std::pair<int, Vector<double> > prev_x_;
  std::pair<int, Vector<double> > prev_prev_x_;

  /** bounds for the asymptotes */
  double u_max_;
  double l_min_;

  /** optional refinements in case backtracking returned to old solution. Then it might makes sense to
   * solve the subproblem again with tighter tolerances */
  double max_refine_;
  double refine_steps_;

  /** maximal number of asymptotes reduction in case the subproblem cannot be solved. Enlarging is automatically
   * by UpdateAsymptotes when solutions don't oscillate */
  int max_reductions_;

  /** the KKT criterion for CalcKKT() */
  double kkt_;

  /** do we have a complete set of design bound constraints? */
  bool full_design_bound_constraints_;

};

/** this is either the approximation of a function or the function itself */
class Approximation
{
public:
  /** @param constraint_idx -1 for objective */
  Approximation(FeasPP* feas_pp, int constraint_idx, bool approximate);

  /** needs to be called to finish constructor. But determinant_shift is not known at constructor time! */
  void PostInit();

  typedef enum { FUNC, GRAD, HESSIAN } Eval;

  /** evaluate function value or first derivative according to the MMA approximation */
  double Evaluate(const double* x_inner, Eval eval, StdVector<double>* out = NULL);

  void AddHessianPattern(compressed_matrix<double>& hessian);

  /** gives the position within the gradient for a special design.
   * @return for dense gradients this is design otherwise a search is performed */
  unsigned int FindGradIndex(unsigned int design) const;

  /** When we have benson vanderbei constraints in the outer problem (to calc KKT and augmented Lagrangian), we need
   * solve the subproblem with determinant constraints. Then the Lagrange multipliers obtained from IPOPT need to be
   * transformed according to the feasibility paper. If this dies not apply for this function the input value is
   * returned.
   * You need to make sure, that the current CFS design is the final design of ipopt!! */
  double TransformMultiplyer(double lambda_ipopt);

  /** helper for logging
   * @param determinant @see GetCondition() */
  std::string ToString(bool determinant = false);

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

  /** Remember do reset the condition container via Done()!
   * @see constraint_idx.
   * @param determinant see GetCondition() */
  Function* GetFunction(bool determinant = false);

  /** It is important to use only this method to get a condition and not optimization->conditions.view->Get()!
   * The reason is that we might have to differentiate between determinant constraints and benson vanderbei constraints.
   * Nevertheless, one needs to call optimization->conditions.view->Done() when all constraints are traversed!
   * @param determiant gives corresponding determinant constraint instead of benson vanderbei if it applyies. Otherwise ignored */
  Condition* GetCondition(bool determinant = false);

  /** When this constraint_idx points to a benson vanderbei constraint, the det_shift added to constraint_idx gives
   * the corresponding determinant constraint which is required for the inner problem */
  int determinant_shift;

private:

  double EvalApproximation(const double* x_inner, Eval eval, StdVector<double>* out);

  double EvalDirect(const double* x_inner, Eval eval, StdVector<double>* out);

  /** We cannot store a function pointer because of the local constraints where we need an index.
   * -1 is for objective */
  int constraint_idx;

  FeasPP* common;
};



} // end of namespace


#endif /* FEASPP_HH_ */
