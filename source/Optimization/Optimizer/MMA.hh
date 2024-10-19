#ifndef MMA_HH_
#define MMA_HH_


#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

namespace CoupledField
{
class MMA;
class BFGS;

/** Here we collect the data for the BFGS subsolver. It holds
 * the rather self contained BFGS and keeps other dual data out of
 * MMA. Due to the self contained stuff, it is not the most efficient approach */
class BFGSHolder
{
public:
  BFGSHolder() {};

  ~BFGSHolder();

  /** activate the BFGS engine */
  void Init(MMA* mma);

  /** to be repeatedly called */
  int BFGSSubProblemSolver();

  // To store the primal values for BFGS

  unsigned int nsmax=10;
  double dual_low = 0.0;
  double dual_upp = 1.0e9;
  double dual_init = 0.1;

  /** copy of BFGS::x to serve as interface */
  Vector<double> lambda;

  BFGS* bfgs = NULL;

private:
  Vector<double> lam_v;
  Vector<double> up_lam;
  Vector<double> lo_lam;

  Vector<double> primal_x;
  Vector<double> primal_y;

  // set by Init() in cas
  MMA*  mma = NULL;

};

/** In contrast to BFGSHolder which keeps a general purpose self BFGS solver,
 * we have here Hessian based dual solvers specific for MMA.
 * This is a little bit of copy and paste from BFGSHolder/BFGS but this is
 * the way Chaitanya did the stuff originally.
 * This serves for the interior point solver (from Chaitanya) and a Hessian
 * solver with steepest descent based on mma.py (Fabian, based on Kelley) */
class DualSolver
{
public:
  virtual ~DualSolver() {};

  virtual void Init(MMA* mma);

  /** projected [0,inf] solver with Hessian and steepest decent fallback.
   * The reference implementation is in mma.py.
   * Based on the projected steepest descent solver from Kelley */
  int FallbackHessianSolver();

  void WriteHessianMinorHeader();

  // Lagrange multipliers
  Vector<double> lambda, mu; // vector size m

  /** shall we correct the Hessian trace? Only for IP-Solver */
  bool hess_corr = false;

protected:

  /** projects 0 >= a + beta * b <= 1e20 */
  static void Project(const Vector<double>& a, double beta, const Vector<double>& b, Vector<double>& out);
  static Vector<double> Project(const Vector<double>& a, double beta, const Vector<double>& b);

  /** size primal variables */
  unsigned int n = 0;

  /** primal inequality constraints */
  unsigned int m = 0;

  unsigned int sub_prob_iter = 0;

  // primal variables
  Vector<double> x;
  Vector<double> y; // we don't use it here

  MMA* mma = NULL;

private:


  void WriteHessianMinorLog(int minor, double pgc, const Vector<double>& lmbda, double delta_lmbda, const Vector<double>& old_dir, bool do_hess, int ls_steps, double step);

  Vector<double> y0; // to ignore the y-part

};

/** the IP-Solver from Aage is a primal dual solver */
class PrimalDualSolver : public DualSolver
{
public:

  void Init(MMA* mma) override;

  /** interior point solver implemented by Chaitanya
   * @return number if sub-problem-iterations */
  int IPSubProblemSolver();

  void IPLogFileLine(PtrParamNode iteration);

private:

  void DualLineSearch();

  Vector<double> s; // vector size 2*m

  // we have more SubInfo than subproblem iteration as the iterations of the interior point method are also output
  struct SubInfo {
    Vector<double> lambda; // Lagrange multipler for all constraints
    Vector<double> mu;
    Vector<double> s; // slack variables for all constraints + 1 for IP slack
    double err;
    double epsi; // for the fixed number of subproblems to be solved
    int iter; // for the sub problem IP iteration count
  };

  StdVector<SubInfo> subiters;

};


/** the MMA (method of moving asymptotes) implementation originates from Chaitanya Dev as a
 * hiwi job around 2018.
 *
 * Used papers are:
 * Svanberg, "The method of moving asymptotes—a new method for structural optimization." International journal for numerical methods in engineering 24.2 (1987): 359-373.
 * Svanberg, DCAMM Lecture Notes 1998
 * Aage, Lazarov. "Parallel framework for topology optimization using the method of moving asymptotes." Structural and multidisciplinary optimization 47.4 (2013): 493-505
 */
class MMA : public BaseOptimizer
{
public:
  MMA(Optimization* optimization, PtrParamNode pn);

  ~MMA();

  void PostInit() override;

  /** @see BaseOptimizer::ToInfo() */
  void ToInfo(PtrParamNode pn) override;

  /** @see BaseOptimizer */
  void SolveProblem() override;

  /** see BaseOptimizer::LogFileLine() */
  void LogFileLine(std::ofstream* out, PtrParamNode iteration) override;

  void DescribeProperties(StdVector<std::pair<std::string, std::string> >& map) const override;

  void PythonSetProperty(pyObject* args) override;

  /** inefficiently evaluations dual function, as the primal variables are generated and thrown away */
  double EvalDualFunction(Vector<double> &xin);

  /** inefficiently evaluates gradient of dual function. Better use Gradient of Dual */
  Vector<double> EvalDualGrad(Vector<double> &xin);

  /** set primal variables with the solution from the subproblem */
  void SetPrimalVariables(const Vector<double>& x_in, const Vector<double>& y_in);

    // We skip here final double &z_out from Aage min-max handling
  void PrimalVarFromDualVar(const Vector<double>& lambda_in, Vector<double>& x_out, Vector<double> &y_out) const;

  double FunctionOfDual(const Vector<double> &lambda_in, const Vector<double>& x_in, const Vector<double>& y_in) const;

  void GradientOfDual(const Vector<double>& x_in, const Vector<double>& y_in, Vector<double>& dual_grad_out) const;

  /** from Aage, Lazarov, 2013 for the primal dual solver in section 3.2 and 3.3 */
  void HessianOfDual(const Vector<double>& lambda_in, const Vector<double>& mu_in, const Vector<double>& x_in, Matrix<double>& hess_out) const;

  /** pure dual based on mma.py.
   * does -1 for maximization */
  void HessianOfDual(const Vector<double>& lambda_in, const Vector<double>& x_in, Matrix<double>& hess_out) const;

  /** could actually move to PrimalDualSolver ?! */
  double DualResidual(const Vector<double>& lambda_in, const Vector<double>& mu_in, const Vector<double>& x_in, const Vector<double>& y_in, double epsi) const;

  typedef enum { SVANBERG, TOPOPT_ROBUST_SHORT, TOPOPT_ROBUST_LONG , FIXED } AsymUpdate;

  typedef enum { NEWTON_SOLV, IP_SOLV, BFGS_SOLV} SubSolverType;

  Enum<AsymUpdate> asymUpdate;
  Enum<SubSolverType> subSolverType;

  /** n is the number of primal design variables - value assigned in MMA::PostInit() */
  unsigned int n = 0;
  /** m is the number of  constraints - value assigned in MMA::PostInit() */
  unsigned int m = 0;
  unsigned int max_sub_iter = 20; // maximum iteration for subproblem solver
  double sub_solve_tol = 1.0e-4; // tolerance for subproblem solver

  /** output dual variables per iteration */
  bool verbose_dual_vars = false;

  // In Aage the parameter for the min-max variable z
  // Vector<double> a;
  Vector<double> c; // penalty parameter for the subproblem

  /** when the subproblem fails, the error message set */
  std::string mma_error;

  /** optional write textual .mma log file */
  std::ofstream* log = nullptr;

private:

  /** @see BaseOptimier - different from our own .mma log */
  void LogFileHeader(Optimization::Log& log) override;

  /** header for textual .mma log */
  void WriteMMALogHeader();

  /** we have different data than in the .plot.dat and we cannot be plotted.
   * Call immediately after CommitIteration()! */
  void WriteMMALogMajor();

  bool SolveSubProblem();

  void SetupSubProblem();

  void UpdateGCAsymptotes();
  void UpdateNonGCAsymptotes();

  void ComputeObjectiveConstraintsSensitivities();

  /** allows access to the subsolvers stuff */
  const Vector<double>& GetSubSolverLambda() const;

  /** General Optimization Problem
   *      min compliance(xval)
   *   s.t:   constraints_i <= 0 ; i=[1 , m]
   *          xmin_i <= xval_i <= xmax_i ; i =[1, n]
   */



  /** this is what is called the asymptotes in MMA and low/upp in Svanberg's mmasub.m.\
   * The std_xmin case holds the data for WriteBoundsToExtern. */
  StdVector<double> std_xmin;
  StdVector<double> std_xmax;
  /** is a math object to std_* with no own memory
   * @see Vector::Replace */
  Vector<double> xmin;
  Vector<double> xmax;

  /** objective value at the current iteration point
   * values updated in MMA::ComputeObjectiveConstraintsSensitivities() */
  double obj_val = -1.0;

  /** gradient/sensitivity of objective at the current iteration point
   * values updated in MMA::ComputeObjectiveConstraintsSensitivities() */
  StdVector<double> grad_objective;

  /** constraints evaluated at the current iteration point
  * values updated in MMA::ComputeObjectiveConstraintsSensitivities() */
  StdVector<double> constraints;

  /** gradient/sensitivity of constraints evaluated at the current iteration point
  * values updated in MMA::ComputeObjectiveConstraintsSensitivities() */
  StdVector<double> grad_constraints; // gradient of contraints, dimension m x n.

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

  Vector<double> p_0j, q_0j; // objective approximation
  Matrix<double> p_ij, q_ij; // constraint approximation
  Vector<double> b; // rhs of inequality constraints  in subproblem.

  /** MMA Sub Problem
   *      min funcA(xval) + z + 1/2 z^2 + summation(y_i*c_i + 1/2 y_i^2 ; i = [1, m])
   *   s.t:   consA_i(xval) - a_i*z - y_i <=0 ; i = [1, m]
   *          alf_i <= xval <= bet_i ; i =[1, n]
   *          y_i >= 0
   *          z >= 0
   */

  AsymUpdate asymUpdate_ = SVANBERG; // TOPOPT_ROBUST_SHORT, TOPOPT_ROBUST_LONG, SVANBERG, FIXED
  SubSolverType subSolverType_ = NEWTON_SOLV; //IP , BFGS

  /** Determines how aggresively the asymptotes are moved
   * asymptotes initialization and increase/decrease
   * these values are used in MMA::GenreteSubProblem() to update low and upp
   * can be initiated from xml. See the constructor */
  double asyminit = 0.9;
  double asymdec = 0.7;
  double asyminc = 1.43;
  double ml_asym = 0.1; // for alpha/beta update, called ALBEFA in MMA and GCMMA – Fortran versions March 2013
  double move_limit = 0.5; // for alpha/beta update, called XXMOVE in MMA and GCMMA – Fortran versions March 2013, not given in the 87 paper

  /** for fixed asymptotes only
  * according to K.Svanberg's paper section 3. equation 10.
  * according to the paper the reasonable choice for asym_fixed_lower is 0 and for asym_fixed_upper is 10
  * updated in
  * used in MMA::GenreteSubProblem()*/
  double asym_fixed_lower = 0;
  double asym_fixed_upper = 10;

  double penalty_c = 1000.0;

  /** sub problem move limits
   * according to K.Svanberg's DCAMM lecture notes section 4.
   * this is chosen to avoid division by zero in sub problem, refer K.Svanberg's paper section 3. equation 8
   * values updated in MMA::GenreteSubProblem()
   * values used in MMA::PrimalVarFromDualVar()*/
  Vector<double> alpha, beta;


  // elastic variable to solve non smooth problems like min-max in Agge - not used
  // double z=0.0;

  /** For fixed asymptotes in MMA::GenerateSubProblem()
   * the upper asymptotes can be set by multiplying a constant(asym_fixed_upper) to the current design value
   * or can just be a constant number, this is controlled by this bool, if true we are multiplying by a constant  */
  bool upperMultiplier = true;
  bool lowerMultiplier = true;

  /** When constraintModification = false p_ij and q_ij are formed according to
   * K.Svanberg's DCAMM lecture notes section 4
   * when constraintModification = true p_ij and q_ij are formed according to
   * K.Svanberg's orginal paper*/
  bool kappa = false;

  /** Globally convergent version refer K.Svanberg's DCAMM lecture notes section 6 */
  bool globallyConvergent = false;
  double rho_init = 0.0001; // set in xml in globalyCongergent/rho
  StdVector<double> rho;
  double rho_0=0;
  double objective_r = 0.0; // The r_i in the function approx

  /** Primal Design variable, min bound and max bound */
  Vector<double> xval;

  /** The idea for y and z is provided in K.Svanberg's DCAMM lecture notes section 2.
   * Further, the values are used based on TopOpt implementation,
   * the description for which is N.Aage's paper in section 3.1
   * values updated in MMA::PrimalVarFromDualVar()
   * values used in MMA::GradientOfDual() ; MMA::DualResidual(); */
  Vector<double> y; //elastic variable to make the problem always feasible

  /** Defined in K.Svanberg's paper section 3.
   * values updated in MMA::GenreteSubProblem()
   * and then used in same function to compute p_0j, q_0j, p_ij and q_ij
   * then to compute primal values in MMA::PrimalVarFromDualVar() */
  Vector<double> low, upp; // Asymptotes bound

  /** Used for heuristic to update the asymptotes
   * values update in MMA::SolveMMA()
   * values used in MMA::GenreteSubProblem()*/
  Vector<double> xold1, xold2;
  Vector<double> change; // Design change

  /** becomes active only when we call Init() */
  BFGSHolder bfgs;

  /** robust fallback Newton solver */
  DualSolver dual_solver;
  /** IP-Solver from Aage */
  PrimalDualSolver primal_dual_solver;

  int sub_prob_iter = 0;

  /** shortcut to either &dual_solver (default for bfgs case) or &primal_dual_solver */
  DualSolver* dual = nullptr;

  /** generate sub problem */
  Timer* gsp_timer_ = nullptr;

  /** generate sub problem */
  Timer* sps_timer_ = nullptr;
};


} // end of namespace


#endif /* MMA_HH_ */
