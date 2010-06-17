#ifndef EVALUATEONLY_HH_
#define EVALUATEONLY_HH_

#include "Optimization/Optimizer/BaseOptimizer.hh"

namespace CoupledField
{
  class Optimization;
  class ParamNode;

  /** This class only evaluates the initial design. This is mostly for
   * debugging and export purpose. A further feature is, that we can handle multiple frequencies,
   * such we can do a sweep for the objective value! */ 
  class EvaluateOnly : public BaseOptimizer
  {
  public:
    EvaluateOnly(Optimization* optimization, PtrParamNode pn);
    
    ~EvaluateOnly() {}
    
  protected:

    /** Evaluates the objective function and it's gradient. Same for all constraints
     * and their gradients and the outputs functions once with the inital guess..
     * This is e.g. to test hand made designs and is nicer than max_iterations = 1.
     * A further benefit is, that we can handle multiple frequencies and sweep over the objectives! */
    void SolveProblem();
  };


} // end of namespace



#endif /*EVALUATEONLY_HH_*/
