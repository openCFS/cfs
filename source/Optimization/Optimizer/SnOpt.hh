#ifndef SNOPT_HH_
#define SNOPT_HH_

#include <stdint.h>
#include <string>

typedef int32_t integer;
typedef double doublereal;

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{

class Optimization;
class Timer;

using std::string;

/** Interface class to the SnOpt-Optimizer
 *  @see http://www.sbsi-sol-optimize.com/asp/sol_product_snopt.htm for details */ 
class SnOpt : public BaseOptimizer
{
public:
  /** @param optimization the problem we optimize
   * @param pn here we can have options - might be NULL! */
  SnOpt(Optimization* optimization, PtrParamNode pn);
  
  virtual ~SnOpt();
  
  /** translates the snopt call to a call of cfs-functions */
  int Callback(integer* Status, const integer n,
      doublereal* x, integer* needF, integer* nF, doublereal* F,
      integer* needG, integer* nG, doublereal* G,
      char* cu, integer* lencu,
      integer* iu, integer* leniu,
      doublereal* ru, integer* lenru);

private:
  /** Return the infinty value. According to snopt-manual this has to be at least 1e20 */
  virtual double GetInfBound() const
  {
    return 1e20;
  }
  
  /** mainly calls sninit_ (which must be called prior to anything else from snopt)
   *  and resizes the data StdVectors according to the problem dimensions
   *  Also set snopt options to values defined in xml */
  inline void Init();
  
  /** define the file names where snopt writes the output and open the files
   *  using the snopen_ function defined in SnOptInterface.hh - does not work! problem
   *  comes from snopt which currently (version 7.2) does not respect the settings properly */
  inline void setSnoptOutputFiles();

  /** Does as the name suggests: calls snopta_ to solve the optimization problem */ 
  void SolveProblem();
  
  /** put some interesting information into the info.xml file */
  void InfoXMLOutput();
  
  /** set problem parameters */
  inline bool get_nlp_info();
  
  /** evaluates the cost function */
  inline bool eval_f(int n, const double* x, double& obj_value);
  
  /** evaluates the gradient of the cost function */
  inline bool eval_grad_f(int n, const double* x, double* grad_f);
  
  /** evaluate the constraints */
  inline bool eval_g(int n, const double* x, int m, double* g);
  
  /** evaluate the gradient of the constraints */
  inline bool eval_jac_g(int n, const double* x, int m, int nele_jac, double* values);
  
  /** Use SnOpt-function snmema_ to determine the amount of memory needed.
   *  snmema_ just gives the absolute minimal values, so we overestimate by a 
   *  certain factor. This is atm hardcoded.
   *  If it should not be sufficient for a certain problem, one could put it in xml */
  inline void AdjustWorkArrayMemory();
  
  /** convenience function to collect all option modifications; 
   *  must be called after get_nlp_info! */
  inline void SetSnOptOptions();
  
  /** mainly fills the index StdVectors iGfun, jGvar, iAfun, jAvar */
  inline void initJacobians();
  
  /** if we have linear constraints we only have to set them up once;
   *  this is done in SolveProblem before the call to snopta_ */
  inline void setupLinearConstraints();
  
  /** Helper function for setting an integer option
   *  currently known options are:
   *  minor_iterations_limit: maximum (total!) number of iterations of the simplex method or the QP algo
   *  superbasics_limit: default min{500, #nonlin-vars + 1}, limit of storage allocated for suberbasics variables
   *  timing_level: default = 3, 0 suppresses output
   *  verify_level: \in [-1,.., 3], finite difference checks for gradients
   *  major_print_level: \in {0, 1, 11, 111, 1111, 11111}, controls amount of output
   *  minor_print_level: \in {0, 1, 11}, controls amount of output for QP subproblem
   */
  void SetIntegerValue(const std::string& key, integer value);
  
  /** Helper function for setting a real valued option
   *  currently known options are:
   *  major_optimality_tolerance: final accuracy of dual variables
   *  major_feasibility_tolerance: how accurately should the nonlin constraints be fulfilled
   *  minor_feasibility_tolerance: tolerance with which all variables satisfy upper/lower bounds
   *  linesearch_tolerance: \in [0, 1], default 0.9, higher is less accurate
   */
  void SetNumericValue(const std::string& key, double value);
  
  /** Helper function for setting a string valued option 
   *  currently known options are: 
   *  qpsolver: can be one of Cholesky, CG, QN 
   *    from the snopt-manual: 
   *      - Cholesky: most robust, but expensive if number of superbasics is large
   *      - CG: good for large number of degrees of freedom (superbasics > 2000)
   *      - QN: good if number of major iterations is small and superbasics is large
   */
  void SetStringValue(const std::string& key, const std::string& value);

  /** number of function evaluations */
  integer f_evals;

  /** number of function evaluations */
  integer g_evals;
  
  /** indicates the usage of names for F rows and cols, value of 1 indicates default names generated by snopt */
  integer nxname;
  integer nFname;
    
  /** snopt variables for output */
  integer iSumm; // choose 6 for output to stdout
  integer iPrint; // choose 6 for output to stdout
  integer iSpecs;
  integer DerOpt;
  
  /** some variables that snopt needs for info passing
   *  they must all have size at least 500 */
  integer minrw;
  integer miniw;
  integer mincw;
  integer lenrw;
  integer leniw;
  integer lencw;
  
  integer n;    // number of variables
  integer nF;   // is the number of problem functions (objective and constraints, linear and nonlinear)
  integer nG;   // number of nonzero nonlinear gradient values
  integer lenG; // dimension of the coordinate arrays iGfun, jGvar
  integer nA;   // number of nonzero linear gradient values
  integer lenA; // dimension of the coordinate arrays iAfun, jAvar
  
  integer nS; // number of superbasic variables
  integer nInf;
  doublereal sInf;
  
  // 0 = cold start
  // 1 = same as 0 but more meaningful when a basis file is given
  // 2 = warm start: xstate and Fstate define a valid starting point
  integer Start;
  
  /** SnOpt is hard to kill, so if we get the request to stop the optimization, we set this trigger
   *  to really end the optimization! */
  bool stop;
 
  /** reports the result of the call to snOptA */
  integer INFO;
  
  /** return value of snOptA */
  integer EXIT;
  
  const integer ObjRow; // which row of F do we use as objective? -> by definition, it will be 1! NOTE: Me must add one to mesh with fortran
  doublereal ObjAdd; // add this value to objective when outputting
  
  /** number of linear constraints */
  int lin_constraints;
  
  /** number of nonlinear constraints */
  int nonlin_constraints;
  
  StdVector<doublereal> rw;
  StdVector<integer> iw;
  StdVector<char> cw; // must be 8 times larger than lencw!!!
  
  StdVector<doublereal> x;     // start values
  StdVector<doublereal> xlow;  // holds the lower bounds on x
  StdVector<doublereal> xupp;   // holds the upper bounds on x
  StdVector<doublereal> xmul;
  StdVector<integer> xstate; // is a set of initial states for each x  (0,1,2,3,4,5)
  
  StdVector<doublereal> F;     // nonlinear objective and constraint terms
  StdVector<doublereal> Flow;  // holds the lower bounds on F
  StdVector<doublereal> Fupp;   // holds the upper bounds on F
  StdVector<doublereal> Fmul;   // set of initial values for the dual variables
  StdVector<integer> Fstate;
  
  // G and A always have n colums!
  StdVector<doublereal> G;  // nonlinear objective and constraint gradient terms
  StdVector<integer> iGfun; // (iGfun(k),jGvar(k)), k = 1,2,...,nG, define the
  StdVector<integer> jGvar; // coordinates of the nonzero nonlinear problem derivatives
  
  StdVector<doublereal> A;     // linear objective and constraint gradient terms
  StdVector<integer> iAfun; // (iAfun(k),jAvar(k)), k = 1,2,...,nA, define the
  StdVector<integer> jAvar; // coordinates of the nonzero linear problem derivatives
 
  /** interface for eval_grad_f and eval_jac_g requires a StdVector */
  StdVector<doublereal> gradhelper;
  
  /** this is the size of the objective gradients */
  unsigned int n_obj_grad;

  /** We do not know when the major iterations and define it as the last function evaluation
   * before a new gradient step is evaluated. Such we can do the CommitIteration() only "post mortem" :(
   * This is therefore the "static" helper to be used in Callback(). */
  bool perform_commit_iteration_;
  
  /** filename of snopt output file */
  std::string outfilename;

  /** Timer for SnOpt */ 
  boost::shared_ptr<Timer> timer_;
  /** Timer for callback function */ 
  boost::shared_ptr<Timer> timer_callback_;
};

} // end of namespace
#endif /*SNOPT_HH_*/
