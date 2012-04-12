#ifndef FEASPP_HH_
#define FEASPP_HH_


#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "Optimization/Optimizer/FeasSubProblem.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{

class Approximation;
class FeasPP;

/** this is a C++ implementation of MMA using IPOPT as subproblem solver. A feature is
 * to forward feasibility constraints to the subproblem, see thesis of Sonja Lehmann.
 * The only dependency is IPOPT and we use the dependencies of the IPOPT CFS solver for
 * the FeasPP subsolver IPOPT. Hence USE_IPOPT needs to be enabled to build FeasPP.
 */

class FeasPP : public BaseOptimizer
{
public:
  FeasPP(Optimization* optimization, PtrParamNode pn);

  virtual ~FeasPP();

  /** initialize problem
   * @see BaseOptimizer::PostInit() */
  void PostInit();

  /** number of total design variables */
  unsigned int n;

  /** number of constraint functions */
  unsigned int m;

  Approximation* obj; // created in PostInit()

  StdVector<Approximation*> constr;

  /** the cfs design variable */
  StdVector<double> x_outer;

  /** lower asymptote */
  StdVector<double> L;

  /** upper asymptote */
  StdVector<double> U;

  /** are there constraints which are not to be approximated?
      The objective can be queried by its own flag.  */
  bool non_approx_constraints;

  bool approx_linear;

  bool approx_feasibility;

protected:

  void SolveProblem();

private:

  /** setup and solve the subproblem by ipopt */
  void SolveSubProblem();

  SmartPtr<FeasSubProblem> ipopt;
};



/** this is either the approximation of a function or the function itself */
class Approximation
{
public:
  /** @param constraint_idx -1 for objective */
  Approximation(FeasPP* feas_pp, int constraint_idx);

  /** evaluate function value or first derivative according to the MMA approximation */
  double Eval(const double* x_inner, bool gradient, StdVector<double>* out = NULL);

  /** helper for logging */
  std::string ToString();

  /** the gradient of the real (outer) function) */
  StdVector<double> outer_grad;

  /** this stores the sparsity pattern, might be dense :) */
  StdVector<unsigned int> pattern;

  /** bounds */
  double lower;
  double upper;

  /** function value of the outer function */
  double outer_val;

  /** shall this function be approximated */
  bool approximate;
private:

  /** Remember do reset the condition container via Done()!
   * @see constraint_idx. */
  Function* GetFunction();

  Condition* GetCondition();

  /** We cannot store a function pointer because of the local constraints where we need an index.
   * -1 is for objective */
  int constraint_idx;

  FeasPP* common;
};



} // end of namespace


#endif /* FEASPP_HH_ */
