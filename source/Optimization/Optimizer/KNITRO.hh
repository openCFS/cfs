#ifndef KNITRO_HH_
#define KNITRO_HH_

#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "knitro.h"
namespace CoupledField
{

class Optimization;


/** Embeds the commercial KNITRO optimizer. See the original documentation. This example uses the names from reverseCommExample.c */
class KNITRO : public BaseOptimizer
{
public:
  /** @param optimization the problem we optimize
   * @param pn here we can have options - might be NULL! */
  KNITRO(Optimization* optimization, PtrParamNode pn);

  virtual ~KNITRO();

protected:

  /** @see BaseOptimizer::SolveProblem() */
  void SolveProblem();

  /** @see BaseOptimizer::GetInfBound() */
  double GetInfBound() const { return KTR_INFBOUND; }
  
  /** @see BaseOptimizer::LogFileHeader() */
  void LogFileHeader(Optimization::Log& log)
  {
    log.AddToHeader("abs_feasibility");
    log.AddToHeader("abs_optimality");
  }

  /** @see BaseOptimizer::LogFileLine() */
  void LogFileLine(std::ofstream* out, PtrParamNode iteration);

private:

  /** KNITRO initialization, called by constructor */
  void Init(PtrParamNode pn);

  void SetParameter(const std::string& name, int value);

  void SetParameter(const std::string& name, double value);

  /** kind of this pointer */
  KTR_context* kc;
  int n;
  int m;
  int nnzJ;
  
  // vectors given in KTR_solve()
  StdVector<double> x;
  StdVector<double> lambda;
  StdVector<double> c;
  StdVector<double> objGrad;
  StdVector<double> jac;
  
  // scalars given in KTR_solve()
  double obj; // objective value


};

} // end of namespace
#endif /*KNITRO_HH_*/
