#ifndef SGPHolder_HH_
#define SGPHolder_HH_

#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "Optimization/Optimizer/SGP.hh"
#include "sgp/sgp.hpp"

namespace CoupledField
{
  /** Copied from IPOPT! This class wraps our SGP implementation to separate the
   * SGP SmaprtPointer implementation from the Optimization stuff */
  class SGPHolder : public BaseOptimizer
  {
  public:
    SGPHolder(Optimization* optimization, PtrParamNode pn);
    virtual ~SGPHolder() {};
    
  protected:

    void SolveProblem();
    
    void PostInitScale();

  private:
    /**
     * Helper function to extract core stiffness tensor from cfs
     */
    Matrix<double> GetCoreStiffnessTensor();

    Vector<double> GetElemVols();

    /** design space */
    DesignSpace* space_;

    /** number of all finite elements */
    unsigned int n_elems_;

    SGP* sgp_;

  };
} // end of namespace

#endif /*SGPHolder_HH_*/
