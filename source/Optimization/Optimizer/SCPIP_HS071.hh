#ifndef SCPIP_HS071_HH_
#define SCPIP_HS071_HH_

#ifdef USE_4_CFS
  #include "Optimization/SCPIPBase.hh"
  using namespace CoupledField;
#else
  #include "SCPIPBase.hh"
#endif

  /** C++ Example NLP for interfacing a problem with IPOPT.
   *  HS071_NLP implements a C++ example of problem 71 of the
   *  Hock-Schittkowsky test suite. This example is designed to go
   *  along with the tutorial document and show how to interface
   *  with IPOPT through the TNLP interface. 
   *
   * Problem hs071 looks like this
   *
   *     min   x1*x4*(x1 + x2 + x3)  +  x3
   *     s.t.  x1*x2*x3*x4                   >=  25
   *           x1**2 + x2**2 + x3**2 + x4**2  =  40
   *           1 <=  x1,x2,x3,x4  <= 5
   *
   *     Starting point:
   *        x = (1, 5, 5, 1)
   *
   *     Optimal solution:
   *        x = (1.00000000, 4.74299963, 3.82114998, 1.37940829)
   *
   *
   */
  class SCPIP_HS071 : public SCPIPBase
  {
  public:
    SCPIP_HS071();

    ~SCPIP_HS071();

    bool get_nlp_info(int& n, int& m, int& nnz_jac_g);

    /** overload this method to return the information about the bound
     *  on the variables and constraints. The value that indicates
     *  that a bound does not exist is specified in the parameters
     *  nlp_lower_bound_inf and nlp_upper_bound_inf.  By default,
     *  nlp_lower_bound_inf is -1e19 and nlp_upper_bound_inf is
     *  1e19. (see TNLPAdapter) */
    bool get_bounds_info(int n, double* x_l, double* x_u,
        int m, double* g_l, double* g_u);

    /** overload this method to return the starting point. The bools
     *  init_x and init_lambda are both inputs and outputs. As inputs,
     *  they indicate whether or not the algorithm wants you to
     *  initialize x and lambda respectively. If, for some reason, the
     *  algorithm wants you to initialize these and you cannot, set
     *  the respective bool to false.
     */
    bool get_starting_point(int n, double* x);

    /** overload this method to return the value of the objective function */
    bool eval_f(int n, const double* x, double& obj_value);

    /** overload this method to return the vector of the gradient of
     *  the objective w.r.t. x */
    bool eval_grad_f(int n, const double* x, double* grad_f);

    /** overload this method to return the vector of constraint values */
    bool eval_g(int n, const double* x, int m, double* g);

    /** overload this method to return the jacobian of the
     *  constraints. The vectors iRow and jCol only need to be set
     *  once. The first call is used to set the structure only (iRow
     *  and jCol will be non-NULL, and values will be NULL) For
     *  subsequent calls, iRow and jCol will be NULL. */
    bool eval_jac_g(int n, const double* x, int m, int nele_jac, double* values);

    void finalize_solution(int status, int n, const double* x, const double* z_L, 
        const double* z_U, int m, const double* g, 
        const double* lambda, Double obj_value);
  };

#endif /*SCPIP_HS071_HH_*/
