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

  StdVector<Condition*> constr;

  /** volume of the structure */
  Double volume, volume_observe;

  /** volume constraint bound */
  Double volume_bound;

  /** value of the compliance */
  Double compliance;

  /** constraint value and gradient of the outer function filtering_gaps*/
  Double filtering_gaps, filtering_gaps_observe;
  StdVector<Double> filtering_gap_grad;
  StdVector<Double> volume_grad;
  Double filtering_gaps_bound;

  /** the cfs design variable */
  Vector<Double> x_outer;

  /** density rho and rotation angle theta, thicknesses s1, s2 */
  StdVector<Double> rho_outer, theta_outer, s1_outer, s2_outer;

  StdVector<Matrix<Double> > E_outer;

  StdVector<Matrix<Double> > E_inner;

  /** lower asymptote */
  StdVector<Matrix<Double> > L;

  /** upper asymptote */
  StdVector<Matrix<Double> > U;

  /** values for lower and upper bound */
  Double lbound, ubound;

  /** convergence tolerance for outer iterations */
  Double tolerance;

  /** convergence tolerance for volume bisection */
  Double volume_tolerance;

  /** tau of model function in subsolve */
  Double tau;

  /** turn of finite differences derivative check*/
  bool derivative_check;

  /** base material matrix */
  Matrix<Double> E_0;

  // penalty parameter
  Double pmin_vol,pmax_vol,ppen_vol,pmini,pmaxi,ppeni,pmin_filt,pmax_filt,ppen_filt;

  // counts outer iterations without any progress
  int worseCounter = 0;

  // bisection iteration counter
  int bisect;

  struct InnerVariable {
    void Read(PtrParamNode pn, DesignSpace* space);
    DesignElement::Type type = DesignElement::NO_TYPE;
    int steps = -1;
    /** As we might use design bounds, this are the lower design bounds. Either from design or the design bound constraints
     * Indices by design type */
    Double lower_bound = -1;
    Double upper_bound = -1;
    Double inc = -1;
  };



  StdVector<InnerVariable> inner_variables;

  //solver type, e.g. brute_force
  std::string solver_type;

  bool HasInnerVar(BaseDesignElement::Type type) const;

  InnerVariable& GetInnerVar(DesignElement::Type type);

  // transfer function for penalization of density
  TransferFunction* tf;

  // merit function value
  Double merit;

  /** designMaterial helper object for two-scale optimization */
  DesignMaterial* helper_dm = NULL;

protected:

  void SolveProblem();

  void ToInfoIter(PtrParamNode pm_iter);

  void ToInfo();


  /** helper for Postinit */
  void SetConfiguration();

  /** kind of inner variables configuration */
  typedef enum {NO_CONF= -1, DENSITY_ROTANGLE = 0, STIFF1_STIFF2 = 1, FOMO = 2, FMO = 3, FOMO_TOP = 4} Configuration;

  Configuration configuration = NO_CONF;
private:

  /** setup and solve the subproblem by ipopt */
  void SolveSubProblem();

  /** updates the design and the outer function values and gradients */
  void UpdateToCurrentStep(bool inner = false);

  /** Performs a gradient check with central difference quotient for necessary derivatives */
  StdVector<Double> GradientCheck(Double & max_grad_error);

  typedef struct
  {
    bool old_point_is_optimal;
    int steps;
    Double stepwidth;
    Double org_dx;
    Double curr_dx;
  } LSR; // line search result

  /** performs backtracking linesearch. obj(x_new) >= obj(x_old)= ob->outer_value
   * Note that the objective is evaluated, the design_id is set but the constraints
   * might not be evaluated -> UpdateSystem().
   * The final design is stored in the CFS-Design!
   * @param the solution of the subproblem.
   * @return number of function evaluations and the norms of the designs*/
  LSR Backtracking(const Vector<Double>& x_old, const Vector<Double>& x_new);

  /** update the asymptotes as long as they are not set to fixed!
   * @param force_reduction to react on subproblem problems */
  void UpdateAsymptotes(const Vector<Double>&x_outer, int iter, bool force_reduction = false);

  /** writes design to rho_outer, theta_outer and E_outer, , E_outer is not updated for inner = true */
  void DesignToOuter(bool inner = false, bool only_update_outer = false);

  /** writes rho_outer, theta_outer and E_outer back to design. If bool filter true only tensor entry designs are updated. */
  void OuterToDesign(bool filter = false);

  /** create 2D rotation matrix for angle alpha */
  void GetRotationMatrix(Double alpha, Matrix<Double> & rot);

  /** create outer derivative from full design gradient */
  void GetOuterDerivative(StdVector<Matrix<Double> > & out, StdVector<Double>  obj_grad);

  /** helper function for derivative check, return only necessary gradient entries from design gradient */
  void GetOuterDerivativeVector(StdVector<Double> & out, StdVector<Double> obj_grad) ;

  /** assume the current design to be FMO tensors and output them */
  void DumpFMPTensors();

  void TestRotation();

  typedef enum { NONE, BACKTRACKING, AUG_LAGRANGIAN } Globalization;

  typedef enum { FIXED, MMA } Asymptotes;

  /** this design history is used and controlled by UpdateAsymptotes() */
  std::pair<int, Vector<Double> > prev_x_;
  std::pair<int, Vector<Double> > prev_prev_x_;

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
  Double Evaluate(const Double* x_inner, Eval eval, StdVector<Double>* out = NULL);

  /** evaluate function according to the SGP approximation */
  Double SubSolve(Eval eval, StdVector<Matrix<Double> > df, Double ppen, StdVector<Double> & s1, StdVector<Double> & s2, StdVector<Double> & theta_outer);

  /** evaluate function according to the SGP approximation for density_rotangle configuration*/
  Double SubSolve_Density_Rotangle(Eval eval, StdVector<Matrix<Double> > df, Double ppen, StdVector<Double> & rho_outer, StdVector<Double> & theta_outer);

  /** evaluate function according to the SGP approximation, parameterization FOMO_Top*/
  Double SubSolve_FOMO_Top(Eval eval, StdVector<Matrix<Double> > df, Double ppen, StdVector<Double> & rho_outer, StdVector<Double> & theta_outer);

  /** evaluate function according to the SGP approximation, parameterization FOMO*/
  Double SubSolve_FOMO(Eval eval, StdVector<Matrix<Double> > df, Double ppen, StdVector<Double> & theta_outer, Double Vloc);

  /** evaluate function according to the SGP approximation, parameterization FMO*/
  Double SubSolve_FMO(Eval eval, StdVector<Matrix<Double> > df, Double ppen, Double Vloc);

  /** helper for logging
   * @param determinant @see GetCondition() */
  std::string ToString();

  /** the gradient of the real (outer) function) */
  StdVector<Double> outer_grad;

  /** this stores the sparsity pattern of the Jacobian, might be dense  */
  StdVector<unsigned int> jac_pattern;

  /** this stores the sparsity pattern of the Hessian. Is diagonal for the approximations.
   * The it has 2 columns and * rows with first column i and second j. For diagonal i = j.
   * Stores only the upper triagonal as the Hessian needs to be symmetric */
  Matrix<unsigned int> hess_pattern;

  /** bounds */
  Double lower;
  Double upper;

  /** shall this function be approximated */
  bool approximate;

  /** Remember do reset the condition container via Done()!
   * @see constraint_idx.
   * @param determinant see GetCondition() */
  Function* GetFunction();


private:

  Double EvalApproximation(Double sum_inner_vars, Eval eval, Matrix<Double> BB, Matrix<Double> E_tmptmp, Double ppen,int index);
  Double CalcAnalyticSol_FOMO_Top(Double &rho1, Double &rho2, Double & rho, Vector<Double> & ev,  Matrix<Double> & ev_vector,  Eval eval, Matrix<Double> BB, Double theta_inner, Double ppen, int index);
  Double CalcAnalyticSol_FOMO(Double &rho1, Double &rho2, Vector<Double> & ev,  Matrix<Double> & ev_vector,  Eval eval, Matrix<Double> BB, Matrix<Double> L, Double theta_inner, Double ppen, Double Vloc);
  Double CalcAnalyticSol_FMO(Double &rho1, Double &rho2, Double & rho3, Vector<Double> & ev,  Matrix<Double> & ev_vector,  Eval eval, Matrix<Double> BB, Matrix<Double> L, Double ppen, Double Vloc);
  Double EvalDirect(const Double* x_inner, Eval eval, StdVector<Double>* out);

  void CalcE_inner(Matrix<Double> & E_inner, Double s1, Double s2,Double theta,Matrix<Double> & tmp);

  void CalcE_inner_Density_Rotangle(Matrix<Double> & E_inner, Matrix<Double> E_0, Double theta_inner);


  /** We cannot store a function pointer because of the local constraints where we need an index.
   * -1 is for objective */
  int constraint_idx;

  SGP* common;
};



} // end of namespace


#endif /* SGP_HH_ */
