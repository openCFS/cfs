  #ifndef BASEOPTIMIZER_HH_
#define BASEOPTIMIZER_HH_

#include "Utils/StdVector.hh"
#include <string>
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{
  class Optimization;
  class ParamNode;
  class ParamNode;
  class Timer;
  class Condition;
  
  /** This is the base class of the optimizer tools.
   * Note that for SCPIP we have multiple inheritance */
  class BaseOptimizer
  {
    /** IPOPT would inherint from BaseOptimizer if there wouldn't be that 
     * SharedPointer stuff. Hence it is a friend */
    friend class IPOPT;
    
  public:
    /** @param optimization this is the actual optimization problem
     * @param pn to hold the complete "optimizer" element!. Not NULL! */
    BaseOptimizer(Optimization* optimization, PtrParamNode pn); 

    virtual ~BaseOptimizer();

    /** This solves the complete Optimization problem by using
     * the CalcObjective, CalcConstraint(), ... from the Optimization problem
     * class. Note that there are implementations in the other base
     * classes of some optimizers */
    void SolveOptimizationProblem();

    /** Defines LogFileLine() */
    virtual std::string LogFileHeader();
    
    /** called by Optimization::CommitIteration(), to be overwritten to add optimizer
     * specific data. Shall match LogFileHeader().Don't add a new-line here!! */
    virtual void LogFileLine(std::ofstream* out, PtrParamNode iteration);

  protected:

    /** This is the specific SolveProblem() implementation. */
    virtual void SolveProblem() = 0;
    
    /** Call this in the optimizer constructor when you have manual_scaling. */
    void PostInit(double manual_scaling, bool no_autoscale = false);

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
     * @param cfs_scale @see EvalObjective() */
    void EvalConstraints(int n, const double* x, int m, bool cfs_scale, double* g);
    
    /** Evaluates the constraint gradients and does reordering if necessary
     * @param cfs_scale @see EvalObjective() */
    void EvalGradConstraints(int n, const double* x, int m, int nentries, bool cfs_scale, StdVector<double>& values);
    
    /** helper function for snopt-optimizer which separates linear and nonlinear constraint gradients */
    int EvalGradConstraints(Condition* g, int start, bool cfs_scale, StdVector<double>& values);
    
    /** Provide Upper and Lower bounds to the optimizer */
    void GetBounds(int n, double* x_l, double* x_u, int m, double* g_l, double* g_u);
    
    /** Return the infinty value (here for ipopt) */
    virtual double GetInfBound() 
    {
      return 1e19;
    }
    
    /** If the actual optimizer is able to handle active sets return here the total number of active
     * constraints (including equality constraints which usually always active).
     * @return < 0 does not implement active sets, >= 0 the current active set */
    virtual int GetCurrenActiveSetSize() const
    {
      return -1;
    }

    /** Combines a design_in with an objective */
    struct DesignMemory
    {
      explicit DesignMemory(int id, double val) : design_id(id), value(val) {}
      int design_id;
      /** is either the objective or a scaling */
      double value;
    };
    
    
    /** The problem with autoscale is that we need the values before initializing the external
     * optimizer and on the other side we want to reuse this prior calcuated values for the optimizer!
     * This becomes even more an issue if we do restarted autoscale when a tolerance band is left.<br>
     * This class holds the autoscaling or manual scaling parameters for either objective or a 
     * constraint (therefore multiple instances of Scale). It does not alter the scalar or gradients
     * but it keeps the designs for reuse. Note that an external optimizer might do objective/ scalar
     * constraint evaluations and gradient evaluations for individual design sets! Therefore we store
     * such many designs. Bu they are small against the full PDE. */
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
      
      std::string ToString();
      
      /** out target for the autoscaled gradient. Not 0.0 means we do autoscale */
      double target;
      
      /** The tolerance as fraction before we rescale. 0.0 is only initial autoscale or no autoscale  */
      double tol;
      
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
    
    Optimization* optimization;

    /** Info Node base */
    PtrParamNode info_;
    
    Scale* objective;
    
    /** This flag indicates a scaling error and request for restart */
    bool restart_requested;

    /** Determine the time spent in the external optimizer.
     * This is SolveProblem minus all evaluations */
    Timer* timer_;
    
  private:
    /** Here we store the objective value for a design. */
    DesignMemory design_;
    
    PtrParamNode optimizer_pn_;
    
  };

} // end of namespace
#endif /* BASEOPTIMIZER_HH_*/
