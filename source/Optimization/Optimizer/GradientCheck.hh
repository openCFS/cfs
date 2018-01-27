#ifndef GRADIENTCHECK_HH_
#define GRADIENTCHECK_HH_

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Optimization/Optimizer/BaseOptimizer.hh"

namespace CoupledField
{
  class DesignElement;
  class DesignSpace;
  class Optimization;

  /** Similar to EvaluateOnly only the first step is evaluated. This is done with finite differences
   * and compared to the gradient methods from the problem */
  class GradientCheck : public BaseOptimizer
  {
  public:
    GradientCheck(Optimization* optimization, PtrParamNode pn);
    
    ~GradientCheck() {}

  protected:

    /** Solves a lot of problems :) */
    void SolveProblem();

  private:
    /** Performs the finite difference evaluation. Does not store. Resets design */
    double PerformFiniteDifferenceEval(DesignElement* de, double current_objective, PtrParamNode cg);
    
    /** our h-spacing from XML or default */
    double h;
    
    /** our order == 1 or 2 */
    double order;
    
    /** special result index for finite differences. -1 = no index!
     * This is where the information is stored in the special result section of the design element */
    int finite_diff_result_index_;

    /** special result index for (fd_grad - alg_grad)/(fd_grad). -1 = no index! */
    int error_result_index_;
    
    /** Shortcut */
    DesignSpace* design;
  };


} // end of namespace

#endif /* GRADIENTCHECK_HH_ */
