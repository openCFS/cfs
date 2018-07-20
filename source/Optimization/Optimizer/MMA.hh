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

  /** General Optimization Problem
   *      min compliance(xval)
   *   s.t:   constrains_i <= 0 ; i=[1 , m]
   *          xmin_i <= xval_i <= xmax_i ; i =[1, n]
   */

  /** n is the number of design variables - value assigned in MMA::PostInit() */
  unsigned int n = 0;
  /** m is the number of  contrains - value assigned in MMA::PostInit() */
  unsigned int m = 0;


  /** Design variable, min bound and max bound */
  StdVector<double> xval;

  /** this is what is called the asymptotes in MMA and low/upp in Svanberg's mmasub.m.\
   * The std_xmin case holds the data for WriteBoundsToExtern. */
  StdVector<double> std_xmin;
  StdVector<double> std_xmax;
  /** is a math object to std_* with no own memory
   * @see Vector::Replace */
  Vector<double> xmin;
  Vector<double> xmax;

  /** compliance (objective) value at the current iteration point
   * values updated in MMA::ComputeObjectiveConstraintsSensitivities() */
  double compliance = -1.0;

  /** gradient/sensitivity of compliance(objective) at the current iteration point
   * values updated in MMA::ComputeObjectiveConstraintsSensitivities() */
  StdVector<double> grad_compliance;

  /** constrains evaluated at the current iteration point
  * values updated in MMA::ComputeObjectiveConstraintsSensitivities() */
  StdVector<double> constrains;

  /** gradient/sensitivity of constrains evaluated at the current iteration point
  * values updated in MMA::ComputeObjectiveConstraintsSensitivities() */
  StdVector<double> grad_constrains; // gradient of contrains, dimension m x n.

  /** ToDO: Idea behind choosing value for scaling is not clear.
   * was copied from TopOpt implementation now NOT USED.
   * was used in MMA::ComputeObjectiveConstraintsSensitivities()
   * I removed it and it still worked */
  double obj_scale = 0;


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

  /** Determines how aggresively the asymptotes are moved
   * asymptotes initialization and increase/decrease
   * these values are used in MMA::GenreteSubProblem() to update low and upp
   * can be initiated from xml. See the constructor */
  double  asyminit = 0.9;
  double asymdec = 0.9;
  double asyminc = 1.1;

  /** controls the seed at which outer limits are update
   * used in MMA::AdjustMoveLimits() */
  double move = 0.2;
  bool moveLimits = true;

  /** for fixed asymptotes only
  * according to K.Svanberg's paper section 3. equation 10.
  * according to paer the resonable choice for asym_fixed_lower is 0 and for asym_fixed_upper is 10
  * updated in
  * used in MMA::GenreteSubProblem()*/
  double asym_fixed_lower = -1;
  double asym_fixed_upper = -1;

  unsigned int max_sub_iter = 10; // maximum iteration for subproblem solver
  StdVector<double> a, c, d; // penalty parameter for the subproblem

  /** sub problem move limits
   * according to K.Svanberg's DCAMM lecture notes section 4.
   * this is chosen to avoid division by zero in sub problem, refer K.Svanberg's paper section 3. equation 8
   * values updated in MMA::GenreteSubProblem()
   * values used in MMA::PrimalVarFromDualVar()*/
  StdVector<double> alpha, beta;

  /** The idea for y and z is provided in K.Svanberg's DCAMM lecture notes section 2.
   * Further, the values are used based on TopOpt implementation,
   * the decription for which is N.Aage's paper in section 3.1
   * values updated in MMA::PrimalVarFromDualVar()
   * values used in MMA::GradientOfDual() ; MMA::DualResidual(); */
  StdVector<double> y; //elastic variable to make the problem always fesible
  double z=0.0; // elastic variable to solve non smooth problems like min-max

  /** Differnet Heuristic to update low and upp in subproblem
   * used in MMA::GenreteSubProblem() */
  bool robustAsymptotes = false;
  bool fixedAsymptotes = false;
  /** For fixed asymptotes in MMA::GenreteSubProblem()
   * the upper asymptotes can be set by multiplying a constant(asym_fixed_upper) to the current design value
   * or can just be a costant number, this is controlled by this bool, if true we are multplying by a constant  */
  bool upperMultiplier = true;
  bool lowerMultiplier = true;

  /** When constraintModification = false p_ij and q_ij are formed according to
   * K.Svanberg's DCAMM lecture notes section 4
   * when constraintModification = true p_ij and q_ij are formed according to
   * K.Svanberg's orginal paper*/
  bool constraintModification = false;

  /** Globally convergent version refer K.Svanberg's DCAMM lecture notes section 6
   * ToDo: NOT IMPLEMENTED*/
  // bool globallyConvergent = false;
  // StdVector<double> rho;


  // Lagrange multipliers
  StdVector<double> lamda, mu; // vector size m
  StdVector<double> s; // vector size 2*m

  /** Defined in K.Svanberg's paper section 3.
   * values updated in MMA::GenreteSubProblem()
   * and then used in same function to compute p_0j, q_0j, p_ij and q_ij
   * then to compute primal values in MMA::PrimalVarFromDualVar() */
  StdVector<double> low, upp; // Asymptotes bound

  /** Used for heuristic to update the asymptotes
   * values update in MMA::SolveMMA()
   * values used in MMA::GenreteSubProblem()*/
  StdVector<double> xold1, xold2;
  StdVector<double> change; // Design change

  /** used in MMA::SolveSubProblem() */
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

  /*
  void StupidTest();
  void InitilizeFromFile(std::string filename, double * vec);
  bool is_number(const std::string& s);
  */
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
