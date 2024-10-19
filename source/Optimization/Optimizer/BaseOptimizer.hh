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
    /** same situation as with IPOPT*/
    friend class SGP;
    
  public:

    /* snopt needs separated evalutations for the linear and nonlinear constraint gradients
     * all other optimizers always need everything */
    typedef enum { ALL, LINEAR, NONLINEAR } GradientType;

    /** call PostInit() afterwards!
     * Child classes shall optimizer_timer_  either in the constructor or after PostInit().
     * @param optimization this is the actual optimization problem
     * @param pn to hold the complete "optimizer" element!. Not NULL! */
    BaseOptimizer(Optimization* optimization, PtrParamNode pn, Optimization::Optimizer type);

    virtual ~BaseOptimizer();

    /** call this after the constructor. Make sure optimizer_timer_ is stopped here and if you don't implement it, stop it in the constructor */
    virtual void PostInit();

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
    
    /** optionally adds some general information after initialization and after solution */
    virtual void ToInfo(PtrParamNode pn);

    /** Evaluates the objective. In the autoscale case checks for old value.
     * @param cfs_scale if true use cfs scaling values, not the optimizer values.
     *        Do it also for the gradient!!
     * @return the objective */
    double EvalObjective(int n, const double* x, bool cfs_scale);

    /* Evaluates the objective gradient. In the autoscale case checks for old evaluation first.
     * Also does an EvalObjective() implicit!
     * @param cfs_scale @see EvalObjective()
     * @return true if within autoscale tolerance - false if restart necessary */
    bool EvalGradObjective(int n, const double* x, bool cfs_scale, StdVector<double>& grad_f);

    /** Evaluates the constraints or rather passes them on to the optimization
     * @param cfs_scale @see EvalObjective()
     * @param normalize transform to g <= 0 constraint. */
    void EvalConstraints(int n, const double* x, int m, bool cfs_scale, double* g, bool normalize);
    
    /** helper for EvalConstraints() assuming the design to be set
     * @param normalize @see EvalConstraints()
     * @param direct_call if false we don't touch timers as we assume to be called by  EvalConstraints() */
    double EvalConstraint(Condition* g, bool cfs_scale, bool normalize, bool direct_call = true, Excitation* ev_only_excite = NULL);

    /** Evaluates the constraint gradients
     * @param cfs_scale @see EvalObjective()
     * @param nonlin_only snopt makes a difference between linear and nonlinear constraints and only
     *  need evaluation for the nonlinear part */
    void EvalGradConstraints(int n, const double* x, int m, int nentries, bool cfs_scale, bool normalize,
        StdVector<double>& values, GradientType grtype = ALL);

    /** Helper for EvalGradConstraint()
     * Called directly by FeasPP
     * @return nnz of constraint
     * @param direct_call @see EvalConstraint() */
    int EvalGradConstraint(Condition* g, int start, bool cfs_scale, bool normalize, StdVector<double>& values, bool direct_call = true);

    /** return the objective value. useful for multi objective */
    double GetObjectiveValue() const { return design_.value; }

    /** Return the infinity value (here for ipopt) */
    virtual double GetInfBound() const { return 1e19; }

    PtrParamNode GetInfoNode() { return info_; }

    Optimization* optimization;

    /** returns the eval_[grad]_obj or eval_[grad]_const_timer_ or NULL if none is running */
    boost::shared_ptr<Timer> GetRunningEvalTimer();

    /** validate that no main timer within optimization is running. Make virtual to add optimizer local timers.
     * Has internal asserts().
     * To be called by assert in Optimization, ErsatzMaterial, Design, ... */
    bool ValidateTimers();

    /** this is a helper to call Optimization::CommitIteration() which switches of the optimizer_timer and does timer validation.
     * This shall be protected but we need it for IPOPT */
    void CommitIteration();

    /** add some descriptive properties for python. Default is nothing added */
    virtual void DescribeProperties(StdVector<std::pair<std::string, std::string> >& map) const { }

    /** allow Python to set a parameter for the current optimizer.
     * Only very few Optimizers will implement the feature (e.g. OCM)
     * @param tuple of two strings: parameter, value
     * @see Optimization::PythoGetOptimizerProperties() */
    virtual void PythonSetProperty(pyObject* args)
    {
      throw "optimizer_set_property() not implement by " + Optimization::optimizer.ToString(type_);
    }

  protected:

    /** This is the specific SolveProblem() implementation. */
    virtual void SolveProblem() = 0;

    /** Call this in the optimizer constructor when you have manual_scaling. */
    void PostInitScale(double manual_scaling, bool no_autoscale = false);

    /** Provide Upper and Lower bounds to the optimizer.
     * Note that snopt is able to do sparse linear abs functions like slope constraints by setting upper and lower bounds */
    void GetBounds(int n, double* x_l, double* x_u, int m, double* g_l, double* g_u);
    
    /** Provide upper and lower bounds on constraints to the optimizer.
     * This function is also called by GetBounds(), but we also want to allow that it can be called independently from GetBounds() */
    void GetConstraintsBounds(int m, double* g_l, double* g_u);

    /** If the actual optimizer is able to handle active sets return here the total number of active
     * constraints (including equality constraints which usually always active).
     * @return < 0 does not implement active sets, >= 0 the current active set */
    virtual int GetCurrenActiveSetSize() const { return -1; }

    /** Combines a design_in with an objective */
    struct DesignMemory
    {
      explicit DesignMemory(int id, double val) : design_id(id), value(val), gradient_design_id(id) {}
      int design_id;
      /** is either the objective or a scaling */
      double value;
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
      Scale(BaseOptimizer* base, PtrParamNode autoscale, double manual_scale, bool no_autoscale);
      
      void PostInit();
      
      /** This is an potentially expensive method and calculates the autoscale again
       * for the restarted case. Otherwise you don't need it. */
      void CalcAutoscale();
      
      /** Checks the current gradients in the design space against the tolerance.
       * Calculates opt_scaling if target is given!
       * @return true if not active or no tolerance or within tolerance */
      bool CheckScaling(int n, StdVector<double>& grad);
      
      /** Did we do autoscale? Interesting for iteration-0 commit */
      bool DoAutoscale()
      {
        return autoscale_;
      }

      void ToInfo(ParamNode* in);

      std::string ToString();
      
      /** out target for the autoscaled gradient. Not 0.0 means we do autoscale */
      double target;
      
      /** The tolerance as fraction before we rescale. 0.0 is only initial autoscale or no autoscale  */
      double tol;
      
      /** do we do logarithmic scaling */
      bool logscale;

      /** do we do scaling by a manually given factor */
      double manual;

      /** The optimal scaling for the design */
      DesignMemory opt_scaling;
      
      /** The scaling and the design when it was set (which opt_scaling for autoscale) */
      DesignMemory scaling;
      
      /** max(abs(grad)) */
      DesignMemory current;
      
    private:

      bool autoscale_;

      BaseOptimizer* base_;
    }; // end Scale

    /** Tuned is meant for the projection filter (Heaviside, tanh) with tune beta.
     * It takes the value to update move_limit. See Dunning & Wein, Achieving near
     * binary topology optimization solutions via automatic projection parameter increase */
    struct Tuned
    {
    public:
      /** Will search for a Tune with beta. Throw error, if not found.
       * @param pn is 'tuned' element
       * @param max and divider are default values to be overwritten by pn */
      Tuned(PtrParamNode pn, double* value, double max, double divider, BaseOptimizer* base);

      void ToInfo(PtrParamNode info);

      static const std::string MOVE_LIMIT; // initialized with "move_limit"

      /** query beta and calculate move limit. Call right after Tune::Update()  */
      void Update();

      double max = -1;
      double min = 1e-4;
      double divider = -1;

    private:
      double *value = nullptr;
      Tune* tune = nullptr;
      BaseOptimizer* base = nullptr;
    };

    /** out type */
    Optimization::Optimizer type_;

    /** Info Node base  for Optimizer */
    PtrParamNode info_;

    /** when set, ocm and (dumas_)mma calculate beta driven move limits */
    Tuned* tuned = nullptr;
    
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

    /** this is the link the specific optimizers param node */
    PtrParamNode this_opt_pn_;

    /** Here we store the objective value for a design. */
    DesignMemory design_;

  private:

    /** we need the adjoint for gradients only. In case of a line search the state problems are sufficient.
     * If the design changed, the state problems are solved and the objective is evaluated.
     * For the harmonic case we need to do the adjoint with the state as as for multiple frequences the system is reassembled.
     * Also note the multiple sequence issue! */
    bool SolveAdjointProblemsIfNeeded(int n, const double* x, bool cfs_scale);
    
  };

} // end of namespace
#endif /* BASEOPTIMIZER_HH_*/
