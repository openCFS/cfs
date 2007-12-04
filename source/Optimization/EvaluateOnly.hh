#ifndef EVALUATEONLY_HH_
#define EVALUATEONLY_HH_

#include "Optimization/BaseOptimizer.hh"

namespace CoupledField
{
  class Optimization;
  class ParamNode;

  /** This class only evaluates the iniial design. This is mostly for
   * debugging and export purpose */ 
  class EvaluateOnly : public BaseOptimizer
  {
  public:
    EvaluateOnly(Optimization* optimization, ParamNode* pn);
    
    ~EvaluateOnly() {}
    
    /** Evaluates the objective function and it's gradient. Same for all constraints
     * and their gradients and the outputs functions once with the inital guess..
     * This is e.g. to test hand made designs and is nicer than max_iterations = 1 */
    void SolveProblem();
  };


} // end of namespace



#endif /*EVALUATEONLY_HH_*/
