#ifndef SNOPT_HH_
#define SNOPT_HH_

#include "Optimization/BaseOptimizer.hh"
#include "Optimization/SnOptBase.hh"

namespace CoupledField
{
class Optimization;
class ParamNode;

/** Interface class to the SnOpt-Optimizer
 *  @see http://www.sbsi-sol-optimize.com/asp/sol_product_snopt.htm for details */ 
class SnOpt : public BaseOptimizer, public SnOptBase
{
public:
  /** @param optimization the problem we optimize
   * @param pn here we can have options - might be NULL! */
  SnOpt(Optimization* optimization, ParamNode* pn);
  
  ~SnOpt();

  /** mainly calls sninit_ (which must be called prior to anything else from snopt)
   *  and resizes the data vectors according to the problem dimensions
   *  Also set snopt options to values defined in xml */
  void Init();

  /** Does as the name suggests: calls snopta_ to solve the optimization problem */ 
  void SolveProblem();
  
  /** translates the snopt call to a call of cfs-functions */
  int Callback(int* Status, const int n,
      double* x, int* needF, int* nF, double* F,
      int* needG, int* nG, double* G,
      char* cu, int* lencu,
      int* iu, int* leniu,
      double* ru, int* lenru);

private:
  /** Return the infinty value. According to snopt-manual this has to be at least 1e20 */
  virtual double GetInfBound() 
  {
    return 1e20;
  }

  /** set problem parameters */
  inline bool get_nlp_info(int& n, int& nF, int& nG);
  
  /** set the box constraints and the bounds for the rest of the constraints */
  bool get_bounds_info(int n, double* x_l, double* x_u, int m, double* g_l, double* g_u);
  
  bool get_starting_point(int n, double* x);
  
  /** evaluates the cost function */
  bool eval_f(int n, const double* x, double& obj_value);
  
  /** evaluates the gradient of the cost function */
  bool eval_grad_f(int n, const double* x, double* grad_f);
  
  /** evaluate the constraints */
  bool eval_g(int n, const double* x, int m, double* g);
  
  /** evaluate the gradient of the constraints */
  bool eval_jac_g(int n, const double* x, int m, int nele_jac, double* values);
  
  /** Use SnOpt-function snmema_ to determine the amount of memory needed.
   *  snmema_ just gives the absolute minimal values, so we overestimate by a 
   *  certain factor. This is atm hardcoded.
   *  If it should not be sufficient for a certain problem, one could put it in xml */
  void AdjustWorkArrayMemory();
  
  /** convenience function to collect all option modifications; 
   *  must be called after get_nlp_info! */
  void SetSnOptOptions();
 
  /** The optimizer ParamNode - away from the constructor to support restart */
  ParamNode* optimizer_pn_;
};

} // end of namespace
#endif /*SNOPT_HH_*/
