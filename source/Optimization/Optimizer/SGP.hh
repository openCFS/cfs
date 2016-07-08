#ifndef SGP_HH_
#define SGP_HH_


#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "MatVec/Vector.hh"

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

namespace CoupledField
{

class SGPApproximation;
class SGP;
template <class TYPE> class SCRS_Matrix;

/** this is a C++ implementation of SGP using a brute force algorithm as subproblem solver.
 */

class SGP : public BaseOptimizer
{
public:
  SGP(Optimization* optimization, PtrParamNode pn);

  virtual ~SGP();

  /** initialize problem
   * @see BaseOptimizer::PostInit() */
  void PostInit();

  /** number of total design variables */
  unsigned int n;

  /** number of constraint functions
   * @see m_c and m_e */
  unsigned int m;

  /** number of all finite elements */
  unsigned int n_elem;

  SGPApproximation* obj; // created in PostInit()

  StdVector<SGPApproximation*> constr;

  /** the cfs design variable */
  Vector<double> x_outer;

  /** density rho and rotation angle theta */
  StdVector<double> rho_outer, theta_outer;

  StdVector<Matrix<double> > E_outer;

  /** lower asymptote */
  StdVector<Matrix<double> > L;

  /** upper asymptote */
  StdVector<Matrix<double> > U;

  /** values for lower and upper bound */
  double lbound, ubound;

  /** base material matrix */
  Matrix<double> E_0;

  // penalty parameter
  double pmin,pmax;

  // simp penalty parameter
  int penal;

  // volume fraction
  double vvol;

protected:

  void SolveProblem();

  /** @see BaseOptimizer::ToInfo() */
  void ToInfo(PtrParamNode pn);

private:

  /** setup and solve the subproblem by ipopt */
  void SolveSubProblem();

  /** updates the design and the outer function values and gradients */
  void UpdateToCurrentStep();

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

  /** update the asymptotes as long as they are not set to fixed!
   * @param force_reduction to react on subproblem problems */
  void UpdateAsymptotes(const Vector<double>&x_outer, int iter, bool force_reduction = false);

  /** writes design to rho_outer, theta_outer and E_outer */
  void DesignToOuter();

  /** writes rho_outer, theta_outer and E_outer back to design */
  void OuterToDesign();

  /** create 2D rotation matrix for angle alpha */
  void GetRotationMatrix(double alpha, Matrix<Double> & rot);

  /** create outer derivative from full design gradient */
  void GetOuterDerivative(StdVector<Matrix<double> > & out, StdVector<Double>  obj_grad);

  /** assume the current design to be FMO tensors and output them */
  void DumpFMPTensors();

  void TestRotation();

  typedef enum { NONE, BACKTRACKING, AUG_LAGRANGIAN } Globalization;

  typedef enum { FIXED, MMA } Asymptotes;

  /** this design history is used and controlled by UpdateAsymptotes() */
  std::pair<int, Vector<double> > prev_x_;
  std::pair<int, Vector<double> > prev_prev_x_;

};

/** this is either the approximation of a function or the function itself */
class SGPApproximation
{
public:

  SGPApproximation(SGP* sgp);

  /** needs to be called to finish constructor.*/
  void PostInit();

  typedef enum { FUNC, GRAD} Eval;

  /** evaluate function value or first derivative according to the MMA approximation */
  double Evaluate(const double* x_inner, Eval eval, StdVector<double>* out = NULL);

  /** evaluate function according to the SGP approximation */
  double SubSolve(Eval eval, StdVector<Matrix<double> > df, double ppen, StdVector<double> & rho_outer, StdVector<double> & theta_outer);

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

  /**

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

private:

  double EvalApproximation(double rho_inner, double theta_inner, Eval eval, Matrix<double> E_0, Matrix<double> BB, Matrix<double> E_tmptmp, Matrix<double> L0, double ppen);
  double EvalDirect(const double* x_inner, Eval eval, StdVector<double>* out);

  void CalcE_inner(Matrix<double> & E_inner, Matrix<double> E_0, double theta_inner);

  /** We cannot store a function pointer because of the local constraints where we need an index.
   * -1 is for objective */
  int constraint_idx;

  SGP* common;
};



} // end of namespace


#endif /* SGP_HH_ */
