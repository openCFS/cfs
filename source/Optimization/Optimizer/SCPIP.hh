#ifndef SCPIP_HH_
#define SCPIP_HH_

#include <assert.h>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "Optimization/Optimizer/SCPIPBase.hh"

namespace CoupledField
{
class Optimization;


  /** This is an implementation of the C++ wrapper SCPIPBase of the
   * MAA implementation SCPIP by Ch. Zillober */
  class SCPIP : public BaseOptimizer, public SCPIPBase
  {
    public:
      /** @param optimization the problem we optimize
       * @param pn here we can have options - might be NULL!
       * @param type given here to have it transparent with FeasSCP */
      SCPIP(Optimization* optimization, PtrParamNode pn, Optimization::Optimizer type = Optimization::SCPIP_SOLVER);
      
      virtual ~SCPIP();
    
      void PostInit() override;

      void ToInfo(PtrParamNode pn)  override;

    protected:

      /** Solves the problem. All stuff, including evaluations of the state problem is done
       * within this method. 
       * @throws exception when not ok! */ 
      void SolveProblem() override;

      /** overload BaseOptimier method */
      int GetCurrenActiveSetSize() const override
      {
        assert(mactiv >= 0);
        return mactiv;
      }

      /** Overloads the default implementation to allow sparse jacobians */
      void SetConstraintSparsityPattern() override;

      /** This implementation uses get_sparsity_pattern_size to determine nnz_jac_g */
      bool get_nlp_info(int& n, int& m, int& nnz_jac_g) override;
    
      int get_sparsity_pattern_size(int m, int* jac_g_dim)  override;

      /** overload this method to return the information about the bound
       *  on the variables and constraints. The value that indicates
       *  that a bound does not exist is specified in the parameters
       *  nlp_lower_bound_inf and nlp_upper_bound_inf.  By default,
       *  nlp_lower_bound_inf is -1e19 and nlp_upper_bound_inf is
       *  1e19. (see TNLPAdapter) */
      bool get_bounds_info(int n, double* x_l, double* x_u,
                           int m, double* g_l, double* g_u)  override;
    
      /** overload this method to return the starting point. The bools
       *  init_x and init_lambda are both inputs and outputs. As inputs,
       *  they indicate whether or not the algorithm wants you to
       *  initialize x and lambda respectively. If, for some reason, the
       *  algorithm wants you to initialize these and you cannot, set
       *  the respective bool to false.
       */
      bool get_starting_point(int n, double* x) override;
    
      /** overload this method to return the value of the objective function */
      bool eval_f(int n, const double* x, double& obj_value)  override;
    
      /** overload this method to return the vector of the gradient of
       *  the objective w.r.t. x */
      bool eval_grad_f(int n, const double* x, double* grad_f)  override;
    
      /** overload this method to return the vector of constraint values */
      bool eval_g(int n, const double* x, int m, double* g)  override;
    
      /** overload this method to return the jacobian of the
       *  constraints. The vectors iRow and jCol only need to be set
       *  once. The first call is used to set the structure only (iRow
       *  and jCol will be non-NULL, and values will be NULL) For
       *  subsequent calls, iRow and jCol will be NULL. */
      bool eval_jac_g(int n, const double* x, int m, int nele_jac, double* values)  override;
    
      void finalize_solution(int status, int n, const double* x, const double* z_L, 
                             const double* z_U, int m, const double* g, 
                             const double* lambda, double obj_value)  override;

      /** This method is called every time SCPIP returns. It is more common
       * to IPOPT when next_iter is true! */
      bool intermediate_callback(int iter, bool next_iter)  override;

      /** This are the virtual constraint indices of inequality constraints */
      StdVector<unsigned int> ie_idx;

      /** This are the virtual constraint indices of equality constraints */
      StdVector<unsigned int> eq_idx;
  };

} // end of namespace

#endif // SCPIP_HH

