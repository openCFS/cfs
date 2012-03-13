#ifndef IPOPTHOLDER_HH_
#define IPOPTHOLDER_HH_

#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "Optimization/Optimizer/IPOPT.hh"

namespace CoupledField
{
  /** This class wrappes our IPOPT implementation to seperate the
   * IPOPT SmaprtPointer implementation from the Optimization stuff */
  class IPOPTHolder : public BaseOptimizer
  {
  public:
    IPOPTHolder(Optimization* optimization, PtrParamNode pn) : BaseOptimizer(optimization, pn, Optimization::IPOPT_SOLVER)
    {
      ipopt_ = new IPOPT(optimization, this, pn);
    }
    
    virtual ~IPOPTHolder() 
    { 
    }
    
  protected:

    void SolveProblem()
    {
      do
      {
        if(restart_requested)
          ipopt_->Init();
        ipopt_->SolveProblem();
      }
      while(restart_requested);  
    }
    
  private:
    SmartPtr<IPOPT> ipopt_;
  };
} // end of namespace

#endif /*IPOPTHOLDER_HH_*/
