#ifndef SNOPT_HH_
#define SNOPT_HH_

#include <stdint.h>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{

class Optimization;
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
  
  void LogFileHeader(Optimization::Log& log) override;

  void LogFileLine(std::ofstream* out, PtrParamNode iteration) override;


  /** translates the snopt call to a call of cfs-functions */
  int Callback(int32_t* Status, const int32_t n,
      double* x, int32_t* needF, int32_t* nF, double* F,
      int32_t* needG, int32_t* nG, double* G,
      char* cu, int32_t* lencu,
      int32_t* iu, int32_t* leniu,
      double* ru, int32_t* lenru);

  /** this callback is added by a patch, note that the final major is not called */
  int CallbackUsrmjr(int32_t* mMajor);

  /** with our snopt patch we can choose if we want cfs iterations as ncon (gradient evaluations called by counting Callback)
   * or my major (using the patched CallbackUsrmjr). The later does not call for the final solution.
   * When we select MAJOR, nObj is logged, this is more illustrative than nCon but cfs iterations for nObj make no sense */
  typedef enum { NCON = 0, MAJOR } Iteration;

  static Enum<Iteration> iteration;

  void SetMajor(int major);

private:
  /** Return the infinty value. According to snopt-manual this has to be at least 1e20 */
  double GetInfBound() const override
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
  void SolveProblem() override;
  
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
  inline void InitJacobians();
  
  /** if we have linear constraints we only have to set them up once;
   *  this is done in SolveProblem before the call to snopta_ */
  inline void SetupLinearConstraints();
  
  /** Helper function for setting an int32_t option
   *  currently known options are:
   *  minor_iterations_limit: maximum (total!) number of iterations of the simplex method or the QP algo
   *  superbasics_limit: default min{500, nonlin-vars + 1}, limit of storage allocated for suberbasics variables
   *  timing_level: default = 3, 0 suppresses output
   *  verify_level: in [-1,.., 3], finite difference checks for gradients
   *  major_print_level: in {0, 1, 11, 111, 1111, 11111}, controls amount of output
   *  minor_print_level: in {0, 1, 11}, controls amount of output for QP subproblem
   */
  void SetIntegerValue(const std::string& key, int32_t value);
  
  /** Helper function for setting a real valued option
   *  currently known options are:
   *  major_optimality_tolerance: final accuracy of dual variables
   *  major_feasibility_tolerance: how accurately should the nonlin constraints be fulfilled
   *  minor_feasibility_tolerance: tolerance with which all variables satisfy upper/lower bounds
   *  linesearch_tolerance: in [0, 1], default 0.9, higher is less accurate
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
  int f_evals = 0;

  /** number of function evaluations */
  int g_evals = 0;
  
  /** nMajor set by patched callback - final nMajor needs to be incremented manually */
  int major = 0;

  /** indicates the usage of names for F rows and cols, value of 1 indicates default names generated by snopt */
  int32_t nxname;
  int32_t nFname;
    
  /** snopt variables for output */
  int32_t iSumm; // choose 6 for output to stdout
  int32_t iPrint; // choose 6 for output to stdout
  int32_t iSpecs;
  
  /** some variables that snopt needs for info passing
   *  they must all have size at least 500 */
  int32_t minrw;
  int32_t miniw;
  int32_t mincw;
  int32_t lenrw;
  int32_t leniw;
  int32_t lencw;
  
  int32_t n;    // number of variables
  int32_t nF;   // is the number of problem functions (objective and constraints, linear and nonlinear)
  int32_t nG;   // number of nonzero nonlinear gradient values
  int32_t lenG; // dimension of the coordinate arrays iGfun, jGvar
  int32_t nA;   // number of nonzero linear gradient values
  int32_t lenA; // dimension of the coordinate arrays iAfun, jAvar
  
  int32_t nS; // number of superbasic variables
  int32_t nInf;
  double sInf;
  
  // 0 = cold start
  // 1 = same as 0 but more meaningful when a basis file is given
  // 2 = warm start: xstate and Fstate define a valid starting point
  int32_t Start;
  
  /** SnOpt is hard to kill, so if we get the request to stop the optimization, we set this trigger
   *  to really end the optimization! */
  bool stop;
 
  /** reports the result of the call to snOptA */
  int32_t INFO;
  
  /** return value of snOptA */
  int32_t EXIT;
  
  const int32_t ObjRow; // which row of F do we use as objective? -> by definition, it will be 1! NOTE: Me must add one to mesh with fortran
  double ObjAdd; // add this value to objective when outputting
  
  /** number of linear constraints */
  int lin_constraints;
  
  /** number of nonlinear constraints */
  int nonlin_constraints;
  
  StdVector<double> rw;
  StdVector<int32_t> iw;
  StdVector<char> cw; // must be 8 times larger than lencw!!!
  
  StdVector<double> x;     // start values
  StdVector<double> xlow;  // holds the lower bounds on x
  StdVector<double> xupp;   // holds the upper bounds on x
  StdVector<double> xmul;
  StdVector<int32_t> xstate; // is a set of initial states for each x  (0,1,2,3,4,5)
  
  StdVector<double> F;     // nonlinear objective and constraint terms
  StdVector<double> Flow;  // holds the lower bounds on F
  StdVector<double> Fupp;   // holds the upper bounds on F
  StdVector<double> Fmul;   // set of initial values for the dual variables
  StdVector<int32_t> Fstate;
  
  // G and A always have n colums!
  StdVector<double> G;  // nonlinear objective and constraint gradient terms
  StdVector<int32_t> iGfun; // (iGfun(k),jGvar(k)), k = 1,2,...,nG, define the
  StdVector<int32_t> jGvar; // coordinates of the nonzero nonlinear problem derivatives
  
  StdVector<double> A;     // linear objective and constraint gradient terms
  StdVector<int32_t> iAfun; // (iAfun(k),jAvar(k)), k = 1,2,...,nA, define the
  StdVector<int32_t> jAvar; // coordinates of the nonzero linear problem derivatives
 
  /** interface for eval_grad_f and eval_jac_g requires a StdVector */
  StdVector<double> gradhelper;
  
  /** this is the size of the objective gradients */
  unsigned int n_obj_grad;

  /** For iteration_ == NCOM, we do not know when the major iterations and define it as the last function evaluation
   * before a new gradient step is evaluated. Such we can do the CommitIteration() only "post mortem" :(
   * Not necessary when iteration_ == MAJOR. */
  bool perform_commit_iteration_;
  
  /** filename of snopt output file */
  std::string outfilename;

  Iteration iteration_  ;
};

} // end of namespace
#endif /*SNOPT_HH_*/
