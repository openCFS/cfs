  #ifndef BASEOPTIMIZER_HH_
#define BASEOPTIMIZER_HH_

#include <iosfwd>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Optimization/Optimization.hh"

namespace CoupledField {
template <class TYPE> class StdVector;
}  // namespace CoupledField


namespace CoupledField
{
  class Condition;
  class Timer;
  
  /** This is the base class of the optimizer tools.
   * Note that for SCPIP we have multiple inheritance */
  class BaseOptimizer
  {
    /** IPOPT would inherint from BaseOptimizer if there wouldn't be that 
     * SharedPointer stuff. Hence it is a friend */
    friend class IPOPT;
    
  public:

    /* snopt needs separated evalutations for the linear and nonlinear constraint gradients
     * all other optimizers always need everything */
    typedef enum { ALL, LINEAR, NONLINEAR } GradientType;


    /** cass PostInit() afterwards!
     * @param optimization this is the actual optimization problem
     * @param pn to hold the complete "optimizer" element!. Not NULL! */
    BaseOptimizer(Optimization* optimization, PtrParamNode pn, Optimization::Optimizer type);

    virtual ~BaseOptimizer();

    /** call this after the constructor */
    virtual void PostInit() {};

    /** This solves the complete Optimization problem by using
     * the CalcObjective, CalcConstraint(), ... from the Optimization problem
     * class. Note that there are implementations in the other base
     * classes of some optimizers */
    void SolveOptimizationProblem();

    /** Defines LogFileLine() */
    virtual void LogFileHeader(Optimization::Log& log);
    
    /** called by Optimization::CommitIteration(), to be overwritten to add optimizer
     * specific data. Shall match LogFileHeader().Don't add a new-line here!! */
    virtual void LogFileLine(std::ofstream* out, PtrParamNode iteration);
    
    /** optionally adds some general information after initialization. */
    virtual void ToInfo(PtrParamNode pn) { };

    /** Evaluates the objective. In the autoscale case checks for old value.
     * @param cfs_scale if true use cfs scaling values, not the optimizer values.
     *        Do it also for the gradient!!
     * @return the objective */
    Double EvalObjective(int n, const Double* x, bool cfs_scale);

    /* Evaluates the objective gradient. In the autoscale case checks for old evaluation first.
     * Also does an EvalObjective() implicit!
     * @param cfs_scale @see EvalObjective()
     * @return true if within autoscale tolerance - false if restart necessary */
    bool EvalGradObjective(int n, const Double* x, bool cfs_scale, StdVector<Double>& grad_f);

    /** Evaluates the constraints or rather passes them on to the optimization
     * @param cfs_scale @see EvalObjective()
     * @param normalize transform to g <= 0 constraint. */
    void EvalConstraints(int n, const Double* x, int m, bool cfs_scale, Double* g, bool normalize);
    
    /** helper for EvalConstraints() assuming the design to be set
     * @param normalize @see EvalConstraints()
     * @param direct_call if false we don't touch timers as we assume to be called by  EvalConstraints() */
    Double EvalConstraint(Condition* g, bool cfs_scale, bool normalize, bool direct_call = true);

    /** Evaluates the constraint gradients
     * @param cfs_scale @see EvalObjective()
     * @param nonlin_only snopt makes a difference between linear and nonlinear constraints and only
     *  need evaluation for the nonlinear part */
    void EvalGradConstraints(int n, const Double* x, int m, int nentries, bool cfs_scale, bool normalize,
        StdVector<Double>& values, GradientType grtype = ALL);

    /** Helper vor EvalGradConstraint()
     * Called directly by FeasPP
     * @return nnz of constraint
     * @param direct_call @see EvalConstraint() */
    int EvalGradConstraint(Condition* g, int start, bool cfs_scale, bool normalize, StdVector<Double>& values, bool direct_call = true);


    /** Return the infinity value (here for ipopt) */
    virtual Double GetInfBound() const { return 1e19; }

    Optimization* optimization;

    /** standard optimiers (snopt, scpip, ...) just call the BaseOptimizer::Eval*() functions where
     * the optimizer_timer_ is paused. Special optimizers like EvaluateOnly are more direct and need to pause themselves */
    boost::shared_ptr<Timer> GetOptimierTimer() { return optimizer_timer_; }

    /** returns the eval_[grad]_obj or eval_[grad]_const_timer_ or NULL if none is running */
    boost::shared_ptr<Timer> GetRunnungEvalTimer();

  protected:

    /** This is the specific SolveProblem() implementation. */
    virtual void SolveProblem() = 0;

    /** Call this in the optimizer constructor when you have manual_scaling. */
    void PostInitScale(Double manual_scaling, bool no_autoscale = false);


    /** Provide Upper and Lower bounds to the optimizer.
     * Note that snopt is able to do sparse linear abs functions like slope constraints by setting upper and lower bounds */
    void GetBounds(int n, Double* x_l, Double* x_u, int m, Double* g_l, Double* g_u);
    
    /** If the actual optimizer is able to handle active sets return here the total number of active
     * constraints (including equality constraints which usually always active).
     * @return < 0 does not implement active sets, >= 0 the current active set */
    virtual int GetCurrenActiveSetSize() const { return -1; }

    /** Combines a design_in with an objective */
    struct DesignMemory
    {
      explicit DesignMemory(int id, Double val) : design_id(id), value(val), gradient_design_id(id) {}
      int design_id;
      /** is either the objective or a scaling */
      Double value;
      int gradient_design_id;
    };
    
    
    /** The problem with autoscale is that we need the values before initializing the external
     * optimizer and on the other side we want to reuse this prior calcuated values for the optimizer!
     * This becomes even more an issue if we do restarted autoscale when a tolerance band is left.<br>
     * This class holds the autoscaling or manual scaling parameters for either objective or a 
     * constraint (therefore multiple instances of Scale). It does not alter the scalar or gradients
     * but it keeps the designs for reuse. Note that an external optimizer might do objective/ scalar
     * constraint evaluations and gradient evaluations for individual design sets! Therefore we store
     * such many designs. But they are small against the full PDE. */
    class Scale
    {
    public:
      /** This sets all value and prepares everything. Note that you have to do PostInit()! */
      Scale(BaseOptimizer* base, PtrParamNode autoscale, Double manual_scale, bool no_autoscale);
      
      void PostInit();
      
      /** This is an potentially expensive method and calculates the autoscale again
       * for the restarted case. Otherwise you don't need it. */
      void CalcAutoscale();
      
      /** Checks the current gradients in the design space against the tolerance.
       * Calculates opt_scaling if target is given!
       * @return true if not active or no tolerance or within tolerance */
      bool CheckScaling(int n, StdVector<Double>& grad);

      /** Did we do autoscale? Interesting for iteration-0 commit */
      bool DoAutoscale() 
      {
        return autoscale_;
      }
      
      std::string ToString();
      
      /** out target for the autoscaled gradient. Not 0.0 means we do autoscale */
      Double target;
      
      /** The tolerance as fraction before we rescale. 0.0 is only initial autoscale or no autoscale  */
      Double tol;
      
      /** do we do logarithmic scaling */
      bool logscale;

      /** The optimal scaling for the design */
      DesignMemory opt_scaling;
      
      /** The scaling and the design when it was set (which opt_scaling for autoscale) */
      DesignMemory scaling;
      
      /** max(abs(grad)) */
      DesignMemory current;
      
    private:
      
      bool autoscale_;
      
      BaseOptimizer* base_;
    };

    /** out type */
    Optimization::Optimizer type_;

    /** Info Node base */
    PtrParamNode info_;
    
    Scale* objective;
    
    /** This flag indicates a scaling error and request for restart */
    bool restart_requested;

    /** Determine the time spent in the external optimizer.
     * This is SolveProblem minus all evaluations */
    boost::shared_ptr<Timer> optimizer_timer_;
    boost::shared_ptr<Timer> eval_obj_timer_;
    boost::shared_ptr<Timer> eval_grad_obj_timer_;
    boost::shared_ptr<Timer> eval_const_timer_;
    boost::shared_ptr<Timer> eval_grad_const_timer_;

    /** this is the link to the general optimization where we can find autoscale. Is not NULL */
    PtrParamNode gen_opt_pn_;

    /** this is the link to the general optimization where we can find autoscale. Is not NULL */
    PtrParamNode this_opt_pn_;

  private:

    /** we need the adjoint for gradients only. In case of a line search the state problems are sufficient.
     * For the harmonic case we need to do the adjoint with the state as as for multiple frequences the system is reassembled.
     * Also note the multiple sequence issue! */
    bool SolveAdjointProblemsIfNeeded(int n, const Double* x, bool cfs_scale);
    
    /** Here we store the objective value for a design. */
    DesignMemory design_;
    
  };

} // end of namespace
#endif /* BASEOPTIMIZER_HH_*/
