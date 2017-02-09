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
  double volume, volume_observe;

  /** volume constraint bound */
  double volume_bound;

  /** value of the compliance */
  double compliance;

  /** constraint value and gradient of the outer function filtering_gaps*/
  double filtering_gaps, filtering_gaps_observe;
  StdVector<double> filtering_gap_grad;
  StdVector<double> volume_grad;
  double filtering_gaps_bound;

  /** the cfs design variable */
  Vector<double> x_outer;

  /** density rho and rotation angle theta, thicknesses s1, s2 */
  StdVector<double> rho_outer, theta_outer, s1_outer, s2_outer;

  StdVector<Matrix<double> > E_outer;

  StdVector<Matrix<double> > E_inner;

  /** lower asymptote */
  StdVector<Matrix<double> > L;

  /** upper asymptote */
  StdVector<Matrix<double> > U;

  /** values for lower and upper bound */
  double lbound, ubound;

  /** convergence tolerance for outer iterations */
  double tolerance;

  /** convergence tolerance for volume bisection */
  double volume_tolerance;

  /** tau of model function in subsolve */
  double tau;

  /** turn of finite differences derivative check*/
  bool derivative_check;

  /** base material matrix */
  Matrix<double> E_0;

  // penalty parameter
  double pmin_vol,pmax_vol,ppen_vol,pmini,pmaxi,ppeni,pmin_filt,pmax_filt,ppen_filt;

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
    double lower_bound = -1;
    double upper_bound = -1;
    double inc = -1;
  };



  StdVector<InnerVariable> inner_variables;

  //solver type, e.g. brute_force
  std::string solver_type;

  bool HasInnerVar(BaseDesignElement::Type type) const;

  InnerVariable& GetInnerVar(DesignElement::Type type);

  // transfer function for penalization of density
  TransferFunction* tf;

  // merit function value
  double merit;

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
  StdVector<double> GradientCheck(double & max_grad_error);

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

  /** writes design to rho_outer, theta_outer and E_outer, , E_outer is not updated for inner = true */
  void DesignToOuter(bool inner = false, bool only_update_outer = false);

  /** writes rho_outer, theta_outer and E_outer back to design. If bool filter true only tensor entry designs are updated. */
  void OuterToDesign(bool filter = false);

  /** create 2D rotation matrix for angle alpha */
  void GetRotationMatrix(double alpha, Matrix<Double> & rot);

  /** create outer derivative from full design gradient */
  void GetOuterDerivative(StdVector<Matrix<double> > & out, StdVector<Double>  obj_grad);

  /** helper function for derivative check, return only necessary gradient entries from design gradient */
  void GetOuterDerivativeVector(StdVector<double> & out, StdVector<Double> obj_grad) ;

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
  double SubSolve(Eval eval, StdVector<Matrix<double> > df, double ppen, StdVector<double> & s1, StdVector<double> & s2, StdVector<double> & theta_outer);

  /** evaluate function according to the SGP approximation for density_rotangle configuration*/
  double SubSolve_Density_Rotangle(Eval eval, StdVector<Matrix<double> > df, double ppen, StdVector<double> & rho_outer, StdVector<double> & theta_outer);

  /** evaluate function according to the SGP approximation, parameterization FOMO_Top*/
  double SubSolve_FOMO_Top(Eval eval, StdVector<Matrix<double> > df, double ppen, StdVector<double> & rho_outer, StdVector<double> & theta_outer);

  /** evaluate function according to the SGP approximation, parameterization FOMO*/
  double SubSolve_FOMO(Eval eval, StdVector<Matrix<double> > df, double ppen, StdVector<double> & theta_outer, double Vloc);

  /** evaluate function according to the SGP approximation, parameterization FMO*/
  double SubSolve_FMO(Eval eval, StdVector<Matrix<double> > df, double ppen, double Vloc);

  /** helper for logging
   * @param determinant @see GetCondition() */
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

  /** shall this function be approximated */
  bool approximate;

  /** Remember do reset the condition container via Done()!
   * @see constraint_idx.
   * @param determinant see GetCondition() */
  Function* GetFunction();


private:

  double EvalApproximation(double sum_inner_vars, Eval eval, Matrix<double> BB, Matrix<double> E_tmptmp, double ppen,int index);
  double CalcAnalyticSol_FOMO_Top(double &rho1, double &rho2, double & rho, Vector<double> & ev,  Matrix<double> & ev_vector,  Eval eval, Matrix<double> BB, double theta_inner, double ppen, int index);
  double CalcAnalyticSol_FOMO(double &rho1, double &rho2, Vector<double> & ev,  Matrix<double> & ev_vector,  Eval eval, Matrix<double> BB, Matrix<double> L, double theta_inner, double ppen, double Vloc);
  double CalcAnalyticSol_FMO(double &rho1, double &rho2, double & rho3, Vector<double> & ev,  Matrix<double> & ev_vector,  Eval eval, Matrix<double> BB, Matrix<double> L, double ppen, double Vloc);
  double EvalDirect(const double* x_inner, Eval eval, StdVector<double>* out);

  void CalcE_inner(Matrix<double> & E_inner, double s1, double s2,double theta,Matrix<double> & tmp);

  void CalcE_inner_Density_Rotangle(Matrix<double> & E_inner, Matrix<double> E_0, double theta_inner);


  /** We cannot store a function pointer because of the local constraints where we need an index.
   * -1 is for objective */
  int constraint_idx;

  SGP* common;
};



} // end of namespace


#endif /* SGP_HH_ */
