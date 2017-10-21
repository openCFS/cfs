#include <iostream>
#include <cmath>
#include <cassert>

#ifdef USE_4_CFS
  #include "Optimization/SCPIP_HS071.hh"
  using namespace CoupledField;
#else
  #include "SCPIP_HS071.hh"
#endif

SCPIP_HS071::SCPIP_HS071()
{
}

SCPIP_HS071::~SCPIP_HS071()
{
}

bool SCPIP_HS071::get_starting_point(int n, double* x)
{
  assert(n == 4);

  // initialize to the given starting point
  x[0] = 1.0;
  x[1] = 5.0;
  x[2] = 5.0;
  x[3] = 1.0;

  return true;
}

bool SCPIP_HS071::get_bounds_info(int n, double* x_l, double* x_u,
                            int m, double* g_l, double* g_u)
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == 4);
  assert(m == 2);

  // the variables have lower bounds of 1
  for (int i=0; i<4; i++) {
    x_l[i] = 1.0;
  }

  // the variables have upper bounds of 5
  for (int i=0; i<4; i++) {
    x_u[i] = 5.0;
  }

  // the first constraint g1 has a lower bound of 25
  g_l[0] = 25;
  // the first constraint g1 has NO upper bound, here we set it to 2e19.
  // Ipopt interprets any number greater than nlp_upper_bound_inf as
  // infinity. The default value of nlp_upper_bound_inf and nlp_lower_bound_inf
  // is 1e19 and can be changed through ipopt options.
  g_u[0] = 2e19;

  // the second constraint g2 is an equality constraint, so we set the
  // upper and lower bound to the same value
  g_l[1] = g_u[1] = 40.0;

  return true;
}


bool SCPIP_HS071::get_nlp_info(int& n, int& m, int& nnz_jac_g)
{
  // The problem described in HS071_NLP.hpp has 4 variables, x[0] through x[3]
  n = 4;

  // one equality constraint and one inequality constraint
  m = 2;

  // in this example the jacobian is dense and contains 8 nonzeros
  nnz_jac_g = 8;
  
  return true;
}

bool SCPIP_HS071::eval_f(int n, const double* x, double& obj_value)
{
  assert(n == 4);

  obj_value = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];

  std::cout << "eval_f(" << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << ") -> " << obj_value << std::endl;

  return true;
}

bool SCPIP_HS071::eval_grad_f(int n, const double* x, double* grad_f)
{
  assert(n == 4);

  grad_f[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
  grad_f[1] = x[0] * x[3];
  grad_f[2] = x[0] * x[3] + 1;
  grad_f[3] = x[0] * (x[0] + x[1] + x[2]);

  std::cout << "eval_grad_f(" << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << ") -> "
                              << grad_f[0] << ", " << grad_f[1] << ", " << grad_f[2] << ", " << grad_f[3] << std::endl;

  return true;
}

bool SCPIP_HS071::eval_g(int n, const double* x, int m, double* g)
{
  assert(n == 4);
  assert(m == 2);

  g[0] = x[0] * x[1] * x[2] * x[3];
  g[1] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];

  std::cout << "eval_g(" << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << ") -> g_[0] = " << g[0] << " g_[1] = " << g[1] << std::endl;


  return true;
}

bool SCPIP_HS071::eval_jac_g(int n, const double* x, 
                        int m, int nele_jac, double* values)
{
  assert(values != NULL);

  // return the values of the jacobian of the constraints

  values[0] = x[1]*x[2]*x[3]; // 0,0
  values[1] = x[0]*x[2]*x[3]; // 0,1
  values[2] = x[0]*x[1]*x[3]; // 0,2
  values[3] = x[0]*x[1]*x[2]; // 0,3

  values[4] = 2*x[0]; // 1,0
  values[5] = 2*x[1]; // 1,1
  values[6] = 2*x[2]; // 1,2
  values[7] = 2*x[3]; // 1,3

  std::cout << "eval_jac_g(" << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << ") -> " 
                             << values[0] << ", " << values[1] << ", " << values[2] << ", " << values[3] 
                             << " -> " << values[4] << ", " << values[5] << ", " << values[6] << ", " << values[7] << std::endl;
  return true;
}                        

void SCPIP_HS071::finalize_solution(int status, int n, const double* x, const double* z_L, 
                                    const double* z_U, int m, const double* g, 
                                    const double* lambda, double obj_value)
{
  std::cout << "SCPIP finished: f=" << obj_value;
  for(int i = 0; i < m; i++)
    std::cout << " + " << lambda[i] << "*" << g[i];
  std::cout << std::endl;  
}
 


int main(int argc, char* argv[])
{
  SCPIP_HS071 scpip;
  scpip.Initialize();
  scpip.SetStringValue("output_level", "more_detailed_and_intermediate_results");
  // scpip.SetStringValue("optimization_method", "sequential_convex_programming");
  scpip.SolveProblem();
}
