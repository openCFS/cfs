#ifndef IPOPT_HH_
#define IPOPT_HH_

#include "Utils/StdVector.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

// references to the lib
#include "coin-or/IpTNLP.hpp"
#include "coin-or/IpIpoptApplication.hpp"

namespace CoupledField
{
  
    
class Optimization;
class BaseOptimizer;

using namespace Ipopt;

/** <p>Ipopt (inner point optimization) is (c) IBM and written by Carl Laird and
 *  Andreas Waechter (main author). It is free at least for academic use. 
 *  It requires a solver which is normally pardiso.</p>
 * <p>Note that the namespace of the api is Ipopt, not to be confused with this class
 * IPOPT.</p> 
 * <p>This class is a clone from the example MyNLP class</p>
 * <p>This is a complete optimizer. Calling SolveProblem() does all - the interface is
 * via the overloaded TMLP methods which finally leads to domain->SolveProblem() calls and
 * evaluation of SIMP methods or that like.</p> 
 * <p>This compiles with Ipopt version v4.3.3.</p>*/ 
class IPOPT : public TNLP
{
public:
  /** @param optimization the problem we optimize
   * @param pn here we can have options - might be NULL! */
  IPOPT(Optimization* optimization, BaseOptimizer* base, PtrParamNode pn);
  
  virtual ~IPOPT();

  /** This is the actual initialization. Might be called multiple times for restart */
  void Init();

  /** Solves the problem. All stuff, including evaluations of the state problem is done
   * within this method. ipopt calls via the overloaded TNLP method the corresponding methods
   * in optimization.
   * @throws exception when not ok! */ 
  void SolveProblem();

  /** Method to return some info about the nlp */
  bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                    Index& nnz_h_lag, IndexStyleEnum& index_style);

  /** Method to return the bounds for my problem */
  bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                       Index m, Number* g_l, Number* g_u);

  /** Method to return the starting point for the algorithm */
  bool get_starting_point(Index n, bool init_x, Number* x,
                          bool init_z, Number* z_L, Number* z_U,
                          Index m, bool init_lambda,
                          Number* lambda);

  /** Method to return the objective value */
  bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  /** Method to return the gradient of the objective */
  bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  /** Method to return the constraint residuals */
  bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  bool eval_jac_g(Index n, const Number* x, bool new_x,
                  Index m, Index nele_jac, Index* iRow, Index *jCol,
                  Number* values);

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  /*virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values);*/

   /** This is called for each major iteration and we use it to write the
    * current state to the result files */
   bool intermediate_callback(AlgorithmMode mode,
                              Index iter, Number obj_value,
                              Number inf_pr, Number inf_du,
                              Number mu, Number d_norm,
                              Number regularization_size,
                              Number alpha_du, Number alpha_pr,
                              Index ls_trials,
                              const IpoptData* ip_data,
                              IpoptCalculatedQuantities* ip_cq);

   /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
   void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq);

  /** this method allows to do constraint scaling. It is only called by IPOPT if
   *  "nlp_scaling_method" is set to "user-scaling". This is done by this class 
   * automatically if there is a constraint scaling provided. */  
  bool get_scaling_parameters(Number& obj_scaling, 
                              bool& use_x_scaling, Index n, Number* x_scaling,
                              bool& use_g_scaling, Index m, Number* g_scaling);

  /** count objective evaluations, different from major iterations */
  int nObj = 0;

private:
  /** Reference to the problem. We could get it globally but this way it is more explicit */
  Optimization* opt_;

  /** our BaseOptimizer from which we would also inherit it there wouldn't be that SmartPtr stuff */
  BaseOptimizer* base_;
  
  /** The optimizer ParamNode - away from the constructor to support restart */
  PtrParamNode optimizer_pn_;
  
  /** We handle the ipopt framework within this class */
  SmartPtr<IpoptApplication> app;


};


} // end of namespace
#endif /*IPOPT_HH_*/
