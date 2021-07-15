#ifndef SCPIPBASE_HH_
#define SCPIPBASE_HH_

#include <iostream>
#include <string>

#include "General/Exception.hh"

// This Interface is written to be used in Finite Element Simulation software CFS++
// of chair for sensorics at the university of Erlangen in Germany.
//
// The double usage in the C++SCPIP interface is realized via the USE_4_CFS compile
// time flag.
#ifdef USE_4_CFS
  #include "General/Enum.hh"
  #include "Utils/StdVector.hh"

  using namespace CoupledField;
#else
  #include "BasicEnum.hh"
  #include "BasicException.hh"
  #include "BasicStdVector.hh"
#endif

  /** <p>This class is an abstract interface for the SCPIP optimizer by
   * Christian Zillober. The interface is an adaption of the IPOPT structure.
   * Hence it should be easy to switch between these optimizers.</p>
   * 
   * <p><b>References</b>:
   * <ul>
   *   <li>SCPIP is written by Christian Zillober 
   *      (http://www.mathematik.uni-wuerzburg.de/~zillober). This is the
   *       core implementation of the optimizer written in Fortran. You
   *       have to contact Christian Zillober to get this code as it its not
   *       public available. </li>
   *   <li>From IPOPT I use the interface as template and took the example implementation.
   *       It is open source and they ask to give the following reference:  
   *       A. Waechter and L. T. Biegler, 
   *       On the Implementation of a Primal-Dual Interior Point Filter Line Search Algorithm 
   *       for Large-Scale Nonlinear Programming, Mathematical Programming 106(1), pp. 25-57, 2006</li>
   *   <li>SCPIPBase.hh/cc, SCPIP_HS071.hh/cc and scpip30.hh is written by Fabian Wein,
   *       http://www.lse.e-technik.uni-erlangen.de/webneu/_rubric/index.php?rubric=Fabian+Wein
   *       But most naming and comments are copied and pasted from SCPIP and IPOPT without
   *       giving any further references!!</li>
   * </ul>
   * </p>
   * 
   * <p><b>Features:</b>
   * This version offers the features of SCPIP only partially with an object oriented interfaces,
   * as for most options a selection is hard coded (spstrat and linsys). But an extension
   * is not that difficult. Dynamic memory allocation is implemented. 
   * 
   * On the plus side the constraints are now handled like in IPOPT and in the form 
   * g(x) = c and h(x) <= d or h(x) >= d where c and d are constants that need not be 0.
   * 
   * The following IPOPT options are implemented:
   * obj_scaling_factor and 'lp_scaling_method with the value user-scaling.
   * </p>
   *  */
  class SCPIPBase
  {
  public:
    SCPIPBase();

    virtual ~SCPIPBase(); 

    /** does the real construction */
    void Initialize();

    /** Solves the optimization problem. Can take a while!
     * @param fromWarmstart starts scpip wie ierr=-4, see SCPIP docu
     * Returns the last ierr value which can be decoded via ToString()
     * To be overwritten in FeasSCP! */
    virtual int solve_problem(bool fromWarmstart = false);

    /** Determines information about the problem dimensions. Taken from IPOPT.
     * @param n return here the total number of design variables
     * @param m return here the sum of all equality and inequality constraints
     * @param nnz_jac_g return the size of the jacobian matrix. For dense gradients
     *        this is m * n. If and only if nnz_jac_g is returnes smaller than m * n
     *        then get_sparsity_pattern() is called next.
     * @return true if you could provide the information.
     * @see get_sparsity_pattern() */
    virtual bool get_nlp_info(int& n, int& m, int& nnz_jac_g) = 0;

    /** This method does not exist in IPOPT. As SCPIP divides between equality and
     * inequality constraints but this interface hides the difference additional information
     * about the sparsity pattern of the constraint gradients is needed.
     * If you have dense Jacobians there is no need to overload this method.
     * @param m for verification the total number of constraints
     * @param jac_g_dim a vector of size m where all entries are to be set with the number
     *        of nonzero entries in the respective constraint. The maximum value is n for each
     *        entry.
     * @return the sum over jac_g_dim = n * m in the dense case. -1 if you cannot provie this information
     * @see get_nlp_info() */
    virtual int get_sparsity_pattern_size(int m, int* jac_g_dim);

    /** overload this method to return the information about the bound
     *  on the variables and constraints. The value that indicates
     *  that a bound does not exist is specified in the parameters
     *  nlp_lower_bound_inf and nlp_upper_bound_inf.  By default,
     *  nlp_lower_bound_inf is -1e19 and nlp_upper_bound_inf is
     *  1e19. (see TNLPAdapter) */
    virtual bool get_bounds_info(int n, double* x_l, double* x_u,
        int m, double* g_l, double* g_u) = 0;

    /** overload this method to return the starting point. The bools
     *  init_x and init_lambda are both inputs and outputs. As inputs,
     *  they indicate whether or not the algorithm wants you to
     *  initialize x and lambda respepctively. If, for some reason, the
     *  algorithm wants you to initialize these and you cannot, set
     *  the respective bool to false.
     */
    virtual bool get_starting_point(int n, double* x) = 0;

    /** overload this method to return the value of the objective function */
    virtual bool eval_f(int n, const double* x, double& obj_value) = 0;

    /** overload this method to return the vector of the gradient of
     *  the objective w.r.t. x */
    virtual bool eval_grad_f(int n, const double* x, double* grad_f) = 0;

    /** overload this method to return the vector of constraint values */
    virtual bool eval_g(int n, const double* x, int m, double* g) = 0;

    /** overload this method to return the jacobian of the  constraints.
     * In contrast to IPOPT does SCPIP support active set and w.r.t SCPIP active sets
     * are implemented. The memory layout becomes quite complicated for active sets,
     * therefore we request the full data here.
     * It is possible to extend this method by a vector active indicating the active
     * constraints, just note that SCPIP internally makes a difference of equality and
     * inequality constraints and the existing vector active is only for inequality but
     * this interface mixes both up as with IPOPT.
     * @param n number of design variables. Sparse constraint jacobians have less
     * @param x the design vector of size n
     * @param nele_jac the nnz pattern of all gradients. if all dense this is m*n */
    virtual bool eval_jac_g(int n, const double* x, int m, int nele_jac, double* values) = 0;

    /** This method is called when the algorithm is complete, so one can
     * store/write the solution. Hence this method is directly called before
     * SolveProblem() returns! 
     * @param status not that this are the SCPIP ierr values! only for 0 good!
     * @param n the design space for x, x_l and z_U
     * @param x the final design variables   
     * @param x_L the lagrange multipliers of the lower bound
     * @param x_U the lagrange multipliers of the upper bound
     * @param m the number of equality and inequality constraints.
     * @param g the lagrange multupliers for the constraints
     * @param obj_value the objective value of the last iterate */ 
    virtual void finalize_solution(int status, int n, const double* x, const double* z_L, 
        const double* z_U, int m, const double* g, 
        const double* lambda, double obj_value) = 0;

    /** Intermediate Callback method for the user. 
     * This method is called every SCPIP return, hence multiple times per iteration.
     * This default implementation does PrintInfo(std::cout) on next_iter.
     * @param iter the current iteration this is info[20-1] and only new if next_iter
     * @param next_iter only true if first time or iter value has changed */  
    virtual bool intermediate_callback(int iter, bool next_iter);

    /** this method allows to do objective and constraint scaling. 
     *  It is only called if "nlp_scaling_method" is set to "user-scaling". This
     * is called after SolveProblem() prior to the eval_*() methods. */ 
    virtual bool get_scaling_parameters(double& obj_scaling, bool& use_g_scaling, int m, double* g_scaling)
    {
      throw Exception("get_scaling_parameters() is called because of 'nlp_scaling_method' but not implemented");
    }


    /** Set a string value -> option */
    void SetStringValue(const std::string& key, const std::string& value);

    void SetIntegerValue(const std::string& key, int value);

    void SetNumericValue(const std::string& key, double value);             

    /** This containts the most important exit states. Make it readable with ToString().
     * It resembles the notation of IPOPT but clearly not the numbering */
    typedef enum { Solve_Succeeded = 0, Maximum_Iterations_Exceeded = 1, LineSearch_Max_Iter = 22, 
      Infeasible = 23, Subproblem_Max_Iter = 31, User_Requested_Stop = 1024,
      Gradients_Return_False = 1025} StateCode; 

      /** This is convenient class serving status/info values like IPOPT does it */
      class Statistic
      {
      public:
        Statistic(SCPIPBase* base)
        {
          this->base = base;
        }

        int IterationCount() { return base->info[20-1]; }

        double FinalObjective() { return base->f_org; }

      private:
        SCPIPBase* base;
      };

      Statistic* Statistics() { return statistic_; }

  protected:

    /** This struct helps converting/normalizing the constraints */
    struct ConstraintShift
    {
      /** the number of the constrain is the index within m. From 0! */
      int number;

      /** is this an equality constraint or an inequality constraint */
      bool equal;

      /** the original value! */
      double lower;

      /** the original value */
      double upper;

      /** org value - shift = 0 */
      double shift;

      /** only for unequal. is 1/-1 AFTER shifting */
      double factor;

      /** shortcut to and redundant with grad_j_size */
      int nnz;

      /** for debugging purpose */
      std::string ToString() const;
    };

    /** Generates an message out of the last status (ierr). */
    std::string ToString(int ierr);

    /** Dumps the info block (info and rinfo). For debugging */
    void PrintInfo(std::ostream& os);

    /** Set the gradient structure according the current active set!
     * Is called on initialization and every time the active set changes.
     * IERN, IECN, EQRN and EQCN are filled.
     * Assumes IERN, IECN, EQRN and EQCN to be properly sized and
     * IELENG, IELPAR and EQLENG, EQLPAR to be set..
     * The default implementation calls SetDenseConstraintGradientPattern()
     * so nothing needs to be done for dense Jacobians.
     * @see PrepareConstraintPattern()
     * @see SetDenseConstraintGradientPattern() */
    virtual void SetConstraintSparsityPattern() {
      SetDenseConstraintGradientPattern();
    }

    /**-- From the IPOPT interface --- */
    /** The number of nonzeros of the constraints jacobian */ 
    int nnz_jac_g;

    /** the number of constraints = mie + meq*/
    int m;

    /** The number of function evaluations */
    int f_evals;

    /** The number of gradient evaluations */
    int grad_evals;

    /** Do we want to to obj scaling? Can come from 'obj_scaling_factor' or
     * from get_scaling_parameters() if 'nlp_scaling_method'='user-scaling'.
     * If provided via both ways the values must be equal!  */
    bool use_obj_scaling;

    /** The obj_scaling value foe objective and gradient. 
     * @see use_obj_scaling */
    double obj_scaling;

    /** This is the original objective value */
    double f_org_unscaled;

    /** do constraint scaling? must be set with 'nlp_scaling_method'='user-scaling' 
     * and is then read via get_scaling_parameters() */
    bool use_g_scaling;

    /** @see use_g_scaling */       
    StdVector<double> g_scaling;

    /** this is our unsorted gradient stuff and unscaled gradient stuff! */
    StdVector<double> g_unscaled;   

    /** this is our unsorted gradient stuff - to be converted to the 
     * equality and inequality stuff from SCPIP.
     * The m constraint results */
    StdVector<double> g;   

    /** The nnz_jac_g IPOPT like constraint gradients */
    StdVector<double> jac_g;

    /** Number of nnz per constraint. Each is n in the dense case */
    StdVector<int> jac_g_size;

    /** We use this only to store the combined gradiend lagrange multipliers
     * for finalize_solution. It is small, hence better than missusing g */
    StdVector<double> y_g;

    /** SCPIP constraints are normalized to 0, but the interface is more general.
     * Here we keep all information to do the normalization */
    StdVector<ConstraintShift> shift;

    /** -- Variables for the SCPIP call --- */

    /** the number of design variables */
    int n;

    /** the number of inequality constraints */
    int mie;

    /** the number of equality constraints */
    int meq;

    /** the number of entries (nnz) of all inequality contraint jacobians.
     * Depends on the active set. */
    int ieleng;

    int ielpar;

    /** the number of entries (nnz) of all equality contraint jacobians.
     * Does not depend on the active set */
    int eqleng;

    /** Dimension of inequality dependent arrays H ORG, Y IE, ACTIVE */
    int iemax;

    /** Dimension of equality dependent arrays G ORG, Y EQ */
    int eqmax;

    /** The design space */
    StdVector<double> x;

    /** lower bound of design space */
    StdVector<double> x_l;

    /** upper bound of design space */
    StdVector<double> x_u;

    /** Contains the objective function value of the last computed iterate.
     * In case of scaling this is scaled! */
    double f_org;

    /** contains the values of the inequality constraints at the last computed iterate */
    StdVector<double> h_org;

    /** contains the values of the equality constraints at the last computed iterate */
    StdVector<double> g_org;

    /** gradient of the objective function */
    StdVector<double> df;

    /** lagrange multipliers of inequality constraints */
    StdVector<double> y_ie;

    /** lagrange multipliers of equality constraints */
    StdVector<double> y_eq;

    /** lagrange multipliers of lower design bounds */
    StdVector<double> y_l;

    /** lagrange multipliers of upper design bounds */
    StdVector<double> y_u;

    /** integer control variables */
    StdVector<int> icntl;

    /** real control variables */
    StdVector<double> rcntl;

    /** integer info */
    StdVector<int> info;

    /** real info */
    StdVector<double> rinfo;

    /** output unit number */
    int nout;

    /** real working array */
    StdVector<double> r_scp;

    /** real working sub array */
    StdVector<double> r_sub;

    /** integer working array */
    StdVector<int> i_scp;

    /** integer working sub array */
    StdVector<int> i_sub; 

    /** indicates active inequality constraints */
    StdVector<int> active;

    /** is constant set to reverse communication */
    int mode;

    /** this return value controls the reverse communication */
    int ierr;

    /** inequality constraints row numbers.*/
    StdVector<int> iern;

    /** inequality constraints column numbers. */
    StdVector<int> iecn;

    /** inequality constraint derivative values */
    StdVector<double> iederv;

    /** equality constraints row numbers. */
    StdVector<int> eqrn;

    /** equality constraints column numbers. */
    StdVector<int> eqcn;

    /** The derivative values of the equality constraints.
     * Is once called eqderv in the SCPIP manual but this is wrong. */
    StdVector<double> eqcoef;

    /** number of active constraints in the subproblem */
    int mactiv;

    /** integer subproblem working array */
    StdVector<int> spiw;

    /** real subproblem working array */
    StdVector<double> spdw;

    /** subproblem solution strategy */
    int spstrat;

    /** linear solver to be used. here dense cholesky solver as we are dense! */
    int linsys;


    /** for some strange reason there is a problem with std::abs(). */
    double abs(double val)
    {
      return val > 0 ? val : -1.0 * val;
    }

    /** Initialize all variables, so we can check that they are really set! */
    void InitVariables();

    /** reserve the fixed space with constant size */
    void AllocateFixed();

    /** reserve the space with is dynamic and where we have to call abstract methods.
     * This is executed in PostInit().
     * To be overwritten in FeasSCP */
    virtual void AllocateProblem();

    /** allocates dynamically based on the info values. Does nothing when there is 
     * no change! */
    void AllocateDynamic(); 

    /** Sets the default parameters. This can be overwritten with Set*Value() */
    void SetDefaultParameters();

    /** Evaluates function and gradient values by calling the abstract
     * IPOPT like methods. Gradients need special care! */
    void EvaluateFunctionValues(); 

    /** Evaluates the objective and constraint gradients */
    bool EvaluateGradients();

    /** This prepares everything what is necessary for SetConstraintSparsityPattern().
     * It uses active and is to be called potentially multiple times when the active
     * set changes.
     * @param inital_call sets also the equality stuff and initializes active and other
     * fixed variables. */
    void PrepareConstraintPattern(bool initial_call);

    /** Sets the dense constraint gradients structure according to the current active set. */
    void SetDenseConstraintGradientPattern();

    /** Helper which copies the element from the ipopt unsorted gradient
     * to the the sorted ones and does eventually a negation.
     * @param ipopt the index within the unsorted jac_g
     * @param scpip the index within iederv or eqcoef
     * @param cs multiplies the value with factor with can be -1 in inequality case */ 
    void CopyConstraintGradient(const double* ipopt, double* scpip, ConstraintShift& cs);

    /** Sets the enums wich do the settings stuff */
    void SetEnums();

    /** Dumps the active array */
    void DumpActive();

    /** calls the virtual finalize_solution. Does the constraint sorting stuff! */
    void CallFinalizeSolution();

    /** The enums for controlling settings via text */

    /** Identifier of the icntl index. */
    Enum<int> icntl_;

    /** The real control values */
    Enum<int> rcntl_;

    /** icntl settings */
    Enum<int> opt_meth_;
    Enum<int> output_level_;
    Enum<int> converge_;

    /** The trivial statistics object */
    Statistic* statistic_;

    /** Do we have to call get_scaling_parameters() ? Determined by 'nlp_scaling_method' set to
     * 'user-scaling'. This is set via SetStringValue() and get_scaling_parameters() is called
     * in SolveProblem(). Note, that the scalings are always correct (default 1.0) */
    bool call_scale_parameters_;
  };

#endif /*SCPIPBASE_HH_*/
