#ifndef FEASSUBPROBLEM_HH_
#define FEASSUBPROBLEM_HH_

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "MatVec/Vector.hh"

// references to the lib
#include "coin/IpTNLP.hpp"
#include "coin/IpIpoptApplication.hpp"

namespace CoupledField
{
class FeasPP;

using namespace Ipopt;

/** see IPOPT */
class FeasSubProblem : public TNLP
{
public:
  FeasSubProblem(FeasPP* feas_pp, PtrParamNode pn);

  virtual ~FeasSubProblem();

  /** This is the actual initialzation. Might be called multiple times for restart */
  void Init();

  /** @param the results and status is written added there
   * @return ony for "" we have solve succeeded!
   * @param refine factors the originally set tolerances. In case linesearch returned to the old point in the next iteration
   * the subproblem is solved again then with tighter tolerances. */
  std::string SolveProblem(PtrParamNode in, double refine);

  /** to be called when the first problem is solved */
  void EnableWarmstart();

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
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL) */
  bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values);

   /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
   void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq);

   /** design from finalize_solution */
   Vector<double> x_final;

   /** the Lagrange multipliers. Called y in the feasibility paper. This is the transformed lambda!
    * @see MMAApproximation::TransformMultiplyer() */
   Vector<double> lambda;

   /** This is the original, not-transfomed lambda!! We need it for the warmstart! */
   Vector<double> orig_lambda;

   /** this are the original tolerances, modified in the refine case */
   double org_tol;
   double org_constr_viol_tol;
   double org_acceptable_constr_viol_tol;

private:

  FeasPP* feas_pp;

  PtrParamNode ipopt;

  /** temporary array */
  StdVector<double> tmp_;

  /** data for warmstart. Set after solving the problem */
  StdVector<double> z_l_;
  StdVector<double> z_u_;

  /** We handle the ipopt framework within this class */
  SmartPtr<IpoptApplication> app;

};

} // end of namespace

#endif /* FEASSUBPROBLEM_HH_ */
