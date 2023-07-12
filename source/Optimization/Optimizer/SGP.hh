#ifndef SGP_HH_
#define SGP_HH_

#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "MatVec/Vector.hh"

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

// check if Intel MKL is available
#include <def_use_blas.hh>
#include "sgp/sgp.hpp"

namespace CoupledField
{
class BaseOptimizer;
/** this is an interface to the external SGPInterface solver for simultaneous topology and local orientation optimization
 */

class SGP : public sgp::Algorithm
{
public:
  SGP(Optimization* optimization, BaseOptimizer* base, PtrParamNode pn, sgp::InitParams initargs);

  void Init();

  virtual ~SGP();

  /** initialize problem
   * @see BaseOptimizer::PostInit() */
  void PostInit();

  void eval_Gamma_epsilon(const sgp::Material_Vector& x, std::vector<sgp::Matrix>& Gamma, std::vector<sgp::Vector>& epsilon) const override;
  void eval_dJ(const sgp::Material_Vector& x, std::vector<sgp::Matrix>& dJ) const override;

  /** diff_filt: for visualization of x - F*x, where F is filter matrix */
  double eval_J(const sgp::Material_Vector& x) const override;

  /** helper function to pass detailed filtering info to cfs, such that it can be requested in xml-file via optResult*/
  void store_filt_results(const sgp::Material_Vector& x_filt_sin, const sgp::Material_Vector& x_filt_cos, const sgp::Material_Vector& diff_filt_sin, const sgp::Material_Vector& diff_filt_cos) const override;
  bool abort(const sgp::Info& info) const override;

  void SolveProblem();

  static void SetEnums();

  typedef enum { ZERO_ASYMPTOTES, EXACT } Approximation;
  typedef enum { LINEAR, NON_LINEAR } Filtering;

  /** to convert string/enum for this type */
  static Enum<Approximation> approximation;
  static Enum<Filtering> filtering;

protected:

  void ToInfo(PtrParamNode summary, sgp::Info info);

private:

  /** Initialize external sgp's data structure */
  void InitSGP();

  /** setup and solve the subproblem by SGP */
  void SolveSubProblem();

  /** Performs a gradient check with foward finite difference approximation of tensor derivative for each element*/
  void GradientCheck();

  /** pass additional info, e.g. filter difference x - F*x, to cfs
   * helper func - called in eval_J() */
  void WriteAdditionalResults();

  /** Sets bound values for design parameters and constraints */
  void SetBounds();

  /** Transfer initial design from cfs to external sgp lib */
  void SetInitialDesign();

  /** Transfer filter matrices to external lib and disable filtering in cfs*/
  void SetFilterMatrices();

  /** Invertes global stiffness FE matrix by solving K*inv = I
   *  this is expensive as we have to perform this step for ndof right-hand sides
   */
  Matrix<double> GetInverseStiffMat(const Function* f) const;

  /**
   * In dependence on value of approximation_ = {ZERO_ASYMPTOTES, EXACT}:
   * - Computes derivative dJ/dD of const function w.r.t to material tensor
   * - Computes \Gamma (stored in 'deriv') and \epsilon terms for exact approximation
   * Helper function for inherited eval_dJ() and eval_Gamma_epsilon()
   */
  void CalcTensDeriv(StdVector<Matrix<double>>& deriv, StdVector<Vector<double>>& eps) const;

  /**
   * Computes derivate of const function w.r.t to material tensor for element 'elem'
   * Helper function for CalcTensDeriv
   */
  Matrix<double> CalcDJdDe(BaseFE* ptFe, BaseBDBInt* bdb, const Elem* elem, const StdVector<LocPoint>& intPoints, const StdVector<double>& weights, Vector<double>& u_e);

  /**
   * Computes Gamma and epsilon terms for exact approximation:
   * This involves inversion of global FE stiffness matrix
   */
  void ComputeElementGammaEps(Matrix<Double>& gamma, Vector<Double>& eps, const Matrix<double>& invK, const Elem* elem, const StdVector<Matrix<double>> bMat, Vector<double> u) const;

  /**
   * Helper function for ComputeElementGammaEps: computes Gamma on one FE element
   */
  void ComputeElementGamma(Matrix<Double>& gamma, const StdVector<Integer> eqnVec, const Matrix<double>& invK, Matrix<double> Be) const;

  /** Helper function for setting a string valued option
   *  currently known options are:
   *  approximation: can be 'exact' or 'asymptotes'
   *  regularization: can be 'linear_filter or 'nonlinear_filter'
   */
  void SetStringValue(const std::string& key, const std::string& value);

  /** Helper function for setting an integer option
   *  currently known options are:
   *  levels: number of levels of hierarichal grid for subproblem
   *  samples: number of samples per level on hierarichal grid
   *  max_bisections: number of max. bisection steps for volume multiplier
   */
  void SetIntegerValue(const std::string& key, int value);

  /** Helper function for setting a real valued option
   *  currently known options are:
   *  tau_init: initial value for globalization penalty
   *  tau_factor: update factor for globalization penalty tau, e.g. tau *= tau_factor
   *  penalty_filter: penalty value for regularization term
   *  pmin_vol: lower and upper bounds for volume multiplier
   *  lmax_vol: upper bounds for volume multiplier
   */
  void SetNumericValue(const std::string& key, double value);

  /** Disables filtering in cfs by setting type of each filter of each design element to NO_FILTERING
   * We need the filter matrix only for external sgp lib*/
  void DisableFilters();

  /** helper function: returns a list with integration points for FE elem */
  void GetIntPoints(const BaseBDBInt* bdb, const Elem* elem, BaseFE** ptFe, StdVector<LocPoint>& intPoints, StdVector<double>& weights) const;

  /** helper function: calculates tensor derivative */
  Matrix<double> CalcDJdD(BaseBDBInt* bdb, BaseFE* ptFe, const Elem* elem, shared_ptr<ElemShapeMap> esm, const StdVector<LocPoint>& intPoints, const StdVector<double>& weights, const Vector<double>& ue, const Matrix<double>& umat);

  Approximation approximation_ = ZERO_ASYMPTOTES;
  Filtering filtering_ = LINEAR;

  /** perform gradient check via finite differences?*/
  bool gradientCheck_ = false;

  /** number of total design variables */
  unsigned int n_designs_;

  /** number of design variables per FE, e.g. n_vars = 4 if we have density + 3 angles */
  unsigned int n_vars_ = 0;

  /** number of constraint functions
   * @see m_c and m_e */
  unsigned int m;

  /** number of all finite elements */
  unsigned int n_elems_;

  /** number of rows of material (stiffness) tensor (3 for 2D and 6 for 3D) */
  unsigned int dim_tens_;

  /** volume of the structure */
  double volume, volume_observe, volume_cfs;

  /** volume constraint bound */
  double volume_bound_;

  /** lower asymptote */
  double L_ = 0;

  /** values for lower and upper bound */
  double lbound_, ubound_;

  /** number of levels of hierarichal grid for subproblem */
  unsigned int levels_;

  /** number of samples per level on hierarichal grid */
  unsigned int samples_per_level_density_ = 0;
  unsigned int samples_per_level_angle_ = 0;

  /** tolerance for stopping criterion of optimization problem */
  double tolerance_ = 1e-6;
  /** tolerance for stopping criterion of volume bisection */
  double volume_tolerance_ = 1e-6;

  /** number of max. volume bisections */
  unsigned int max_bisections_;

  /** number of max. globalization steps */
  unsigned int max_globs_ = 20;

  /** turn of finite differences derivative check*/
  bool derivative_check_;

  /** core material matrix */
  Matrix<double> E_0_;

  /** lower and upper bound for volume multiplier lambda */
  double pmin_vol_ = 0;
  double pmax_vol_ = 500;

  /** penalty for regularization terms */
  double p_filt_dens_ = 0;
  double p_filt_angle_ = 0;

  /** initial value for globalization penalty tau */
  double tau_init_ = 1e-9;
  double tau_factor_ = 10;

  /** Number of intergration points per element.
   *  This variable is only accessed when approximation_ == EXACT,
   *  then we assume a regular mesh */
  double nip_ = 0;

  /**
   * flags indicating which variables are subject to optimization
   */
  bool hasDensity_;
  bool hasAngle_; // it is possible to only optimize the density
  bool hasRotAngleFirst_; // the angle in 2D or first angle in 3D
  bool hasRotAngleSecond_;
  bool hasRotAngleThird_;


  // stores indices of
  StdVector<int> var_idx_;

  /** should subsolver generate data for 1D plot exact vs. model?*/
  bool generatePlotData_ = false;


  /** algorithmic parameters for external SGP
   * Need to reset this by calling Init() when doing a restart*/

  sgp::SolveParams solve_params_;

  // counts outer iterations without any progress
  int worseCounter = 0;

  /** Pointer to the problem. We could get it globally but this way it is more explicit */
  Optimization* optimization_;

  /** our BaseOptimizer from which we would also but try to avoid double inheritance this way */
  BaseOptimizer* base_;

  /** ParamNode from base */
  PtrParamNode base_pn_;

  /** The optimizer ParamNode - away from the constructor to support restart */
  PtrParamNode optimizer_pn_;

  /** Pointer to optimization node in .info.xml  */
  PtrParamNode sgp_node_;

  /** design space */
  DesignSpace* space_;

  /** Pointer to grid object*/
  Grid* grid_;

}; // end SGP


} // end of namespace


#endif /* SGP_HH_ */
