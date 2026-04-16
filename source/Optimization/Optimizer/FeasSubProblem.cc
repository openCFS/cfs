#include "Optimization/Optimizer/FeasSubProblem.hh"
#include "Optimization/Optimizer/FeasPP.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Condition.hh"
#include "MatVec/Matrix.hh"
#include "General/Exception.hh"
#include "Utils/StdVector.hh"
#include "Utils/ToolsFull.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "coin-or/IpSolveStatistics.hpp"

#include <cstring>


using namespace CoupledField;

EXTERN_LOG(feasPP)

FeasSubProblem::FeasSubProblem(FeasPP* feas_pp, PtrParamNode pn)
{
  this->feas_pp = feas_pp;
  this->ipopt = pn;
  this->org_tol = -1.0;
  this->org_constr_viol_tol = -1.0;
  this->org_acceptable_constr_viol_tol = -1.0;

  Init();
}

void FeasSubProblem::Init()
{
  // smart pointer!
  app = new IpoptApplication();
  // Initialize the IpoptApplication and process the options
  app->Initialize();

  // check for optional paramters
  if(ipopt != NULL)
  {
    ParamNodeList list;
    list = ipopt->GetListByVal("option", "type", "string");

    for(unsigned int i = 0; i < list.GetSize(); i++)
      app->Options()->SetStringValue(list[i]->Get("name")->As<std::string>(), list[i]->Get("value")->As<std::string>());

    list = ipopt->GetListByVal("option", "type", "integer");
    for(unsigned int i = 0; i < list.GetSize(); i++)
      app->Options()->SetIntegerValue(list[i]->Get("name")->As<std::string>(), list[i]->Get("value")->As<Integer>());

    list = ipopt->GetListByVal("option", "type", "real");
    for(unsigned int i = 0; i < list.GetSize(); i++)
      app->Options()->SetNumericValue(list[i]->Get("name")->As<std::string>(), list[i]->Get("value")->As<Double>());
  }

  app->Options()->GetNumericValue("tol", org_tol, "");
  app->Options()->GetNumericValue("constr_viol_tol", org_constr_viol_tol, "");
  app->Options()->GetNumericValue("acceptable_constr_viol_tol", org_acceptable_constr_viol_tol, "");
}

FeasSubProblem::~FeasSubProblem()
{}


std::string FeasSubProblem::SolveProblem(PtrParamNode in, double refine)
{
  // check for warmstart
  bool no_ws = false;
  if(ipopt != NULL && ipopt->HasByVal("option", "name", "warm_start_init_point"))
    no_ws = ipopt->GetByVal("option", "name", "warm_start_init_point")->As<std::string>() == "no";
  if(orig_lambda.GetSize() > 0 && !no_ws)
    app->Options()->SetStringValue("warm_start_init_point", "yes");
  else
    app->Options()->SetStringValue("warm_start_init_point", "no");

  app->Options()->SetNumericValue("tol", org_tol * refine);
  app->Options()->SetNumericValue("constr_viol_tol", org_constr_viol_tol * refine);
  app->Options()->SetNumericValue("acceptable_constr_viol_tol", org_acceptable_constr_viol_tol * refine);

  // this is where the work is done!
  ApplicationReturnStatus status = app->OptimizeTNLP(this);

  std::string msg = "not set";
  switch(status)
  {
  case Solve_Succeeded:
    break;

  case NonIpopt_Exception_Thrown:
    msg = "(CFS) exception thrown when performing IPOPT call. Try cfs with '-f'";
    break;

  case Restoration_Failed:
    msg = "Restoration failed";
    break;

  case Insufficient_Memory:
    msg = "Insufficient memory";
    break;

  case Infeasible_Problem_Detected:
    msg = "Infeasible problem detected";
    break;

  case Maximum_Iterations_Exceeded:
  {
    int value = 0;
    app->Options()->GetIntegerValue("max_iter", value, "");
    msg = "Maximum subrpoblem iterations exceeded (" + boost::lexical_cast<std::string>(value) + ")";
    break;
  }

  case Feasible_Point_Found:
    msg = "Feasible point found";
    break;

  default:
    msg = (status < 0 ? "error: " : "warning: ") + std::to_string(status);
    break;
  }

  in->Get("iters")->SetValue(app->Statistics()->IterationCount());

  if(status != Solve_Succeeded)
  {
    in->Get("ipopt_msg")->SetValue(msg);
    in->Get("converged")->SetValue("no"); // be very strict and fail on warnings!
  }

  return status == Solve_Succeeded ? "" : msg;
}


bool FeasSubProblem::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  n = feas_pp->n;
  m = feas_pp->m;
  assert((int) feas_pp->constr.GetSize() == m);

  nnz_jac_g = 0;
  for(int i = 0; i < m; i++)
    nnz_jac_g += feas_pp->constr[i]->jac_pattern.GetSize();

  nnz_h_lag = feas_pp->hessian->nnz();

  // We use the standard C index style for row/col entries
  index_style = C_STYLE;

  LOG_TRACE2(feasPP) << ":get_nlp_info n <- " << n << "; m <- " << m << " nnz_jac_g <- "
                    << nnz_jac_g << " nnz_h_lag <- " << nnz_h_lag << " index_style <- "
                    << (index_style == C_STYLE ? "C_STYLE" : "FORTRAN_STYLE");
  return true;
}

bool FeasSubProblem::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
{
  LOG_TRACE2(feasPP) << "get_bounds_info: n = " << n << "; m = " << m;
  assert(n == (int) feas_pp->optimization->GetDesign()->full_data.GetSize()); // TODO include aux design!

  // in Svanberg's original paper there are dynamic bounds given in (8)
  double asy = feas_pp->dynamic_design_bounds ? 0.9 : 1.0;
  double var = feas_pp->dynamic_design_bounds ? 0.1 : 0.0;

  for(int i = 0; i < n; i++)
  {
    BaseDesignElement* de = feas_pp->optimization->GetDesign()->full_data[i];

    // because we might use design bound constraints we must not use de.GetLowerBound(), ..

    x_l[i] = std::max(feas_pp->lower_bound[de->GetType()], asy * feas_pp->L[i] + var * feas_pp->x_outer[i]);
    x_u[i] = std::min(feas_pp->upper_bound[de->GetType()], asy * feas_pp->U[i] + var * feas_pp->x_outer[i]);

    LOG_DBG3(feasPP) << "FSP::gbi: i=" << i << " dt=" << de->GetType() << " asy=" << asy << " var=" << var << " lb=" << feas_pp->lower_bound[de->GetType()] << " ub=" << feas_pp->upper_bound[de->GetType()] << " l=" << feas_pp->L[i]
                      << " u=" << feas_pp->U[i] << " x=" << feas_pp->x_outer[i] << " -> x_l=" << x_l[i] << " x_u=" << x_u[i];
  }

  for(int i = 0; i < m; i++)
  {
    MMAApproximation* g = feas_pp->constr[i];

    g_l[i] = g->lower;
    g_u[i] = g->upper;
  }

  LOG_DBG2(feasPP) << "FSP::gbi: x_l=" << StdVector<double>::ToString(n, x_l);
  LOG_DBG2(feasPP) << "FSP::gbi: x_u=" << StdVector<double>::ToString(n, x_u);
  LOG_DBG2(feasPP) << "FSP::gbi: g_l=" << StdVector<double>::ToString(m, g_l);
  LOG_DBG2(feasPP) << "FSP::gbi: g_u=" << StdVector<double>::ToString(m, g_u);

  return true;
}

bool FeasSubProblem::get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda, Number* ptr_lambda)
{
  LOG_TRACE2(feasPP) << "FSP::gsp: n = " << n << "; m = " << m;

  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the others if
  // you wish.

  assert(n == (int) feas_pp->n);
  assert(n == (int) feas_pp->x_outer.GetSize());

  assert(init_x == true);
  memcpy(x, feas_pp->x_outer.GetPointer(), sizeof(double) * n);
  LOG_DBG3(feasPP) << "FSP::gsp: x=" << StdVector<double>::ToString(n, x);

  if(init_z)
  {
    assert(!z_l_.IsEmpty() && !z_u_.IsEmpty());
    memcpy(z_L, z_l_.GetPointer(), sizeof(double) * n);
    memcpy(z_U, z_u_.GetPointer(), sizeof(double) * n);
  }

  if(init_lambda)
  {
    assert(orig_lambda.GetSize() > 0);
    memcpy(ptr_lambda, orig_lambda.GetPointer(), sizeof(double) * m);
  }
  return true;
}

bool FeasSubProblem::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert((int) feas_pp->x_outer.GetSize() == n);
  MMAApproximation* f = feas_pp->obj;

  if(!f->approximate)
    feas_pp->optimization->GetDesign()->ReadDesignFromExtern(x);

  obj_value = f->Evaluate(x, MMAApproximation::FUNC);

  LOG_DBG3(feasPP) << "eval_f: new_x = " << new_x << " obj_value= << " << obj_value << " x=" << StdVector<double>::ToString(n, x);
  return true;
}

bool FeasSubProblem::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  MMAApproximation* f = feas_pp->obj;

  assert((int) feas_pp->x_outer.GetSize() == n);
  assert((int) f->jac_pattern.GetSize() == n); // assume dense!

  if(!f->approximate)
  {
    feas_pp->optimization->GetDesign()->ReadDesignFromExtern(x);
    feas_pp->optimization->GetDesign()->Reset(DesignElement::COST_GRADIENT);
  }

  tmp_.Resize(n);

  f->Evaluate(x, MMAApproximation::GRAD, &tmp_);
  memcpy(grad_f, tmp_.GetPointer(), sizeof(double) * n);

  return true;
}

bool FeasSubProblem::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  LOG_DBG3(feasPP) << "eval_g: n = " << n << "; new_x = " << new_x << "; m = " << m;

  if(feas_pp->non_approx_constraints)
      feas_pp->optimization->GetDesign()->ReadDesignFromExtern(x);

  assert(m == (int) feas_pp->m);
  for(int i = 0; i < m; i++)
    g[i] = feas_pp->constr[i]->Evaluate(x, MMAApproximation::FUNC);

  feas_pp->optimization->constraints.view->Done(); // reset local constraint to global mode

  return true;
}

bool FeasSubProblem::eval_jac_g(Index n, const Number* x, bool new_x,  Index m, Index nele_jac, Index* iRow, Index *jCol, Number* values)
{
  assert(m == (int) feas_pp->m);

  if(feas_pp->non_approx_constraints && values != NULL)
  {
     feas_pp->optimization->GetDesign()->ReadDesignFromExtern(x);
     feas_pp->optimization->GetDesign()->Reset(DesignElement::CONSTRAINT_GRADIENT);
  }

  unsigned int index = 0;
  for(int func = 0; func < m; func++)
  {
    MMAApproximation* f = feas_pp->constr[func];

    if(values != NULL)
    {
      tmp_.Resize(f->jac_pattern.GetSize());
      f->Evaluate(x, MMAApproximation::GRAD, &tmp_);
    }

    for(unsigned int j = 0; j < f->jac_pattern.GetSize(); j++)
    {
      // return the structure of the Jacobian of the constraints ?
      if(values == NULL)
      {
        iRow[index] = func;
        jCol[index] = f->jac_pattern[j];
        LOG_DBG3(feasPP) << "FSP:ejg f=" << f->ToString(true) << " iRow[" << index << "]=" << func << " jCol=" << jCol[index];
      }
      else
      {
        values[index] = tmp_[j];
      }
      index++;
    }
  }

  feas_pp->optimization->constraints.view->Done(); // reset local constraint to global mode

  return true;
}


bool FeasSubProblem::eval_h(Index n, const Number* x, bool new_x, Number obj_factor, Index m, const Number* lambda,
                    bool new_lambda, Index nele_hess, Index* iRow, Index* jCol, Number* values)
{
  assert(feas_pp->hessian != NULL);
  compressed_matrix<double>& H = *(feas_pp->hessian);
  assert((int) H.nnz() == nele_hess);

  if(values == NULL)
  {
    assert(iRow != NULL && jCol != NULL && x == NULL && lambda == NULL);

    // IPOPT does not use the classical CRS format but the same as we have in MMAApproximation::hess_pattern
    assert((int) H.size1() == n);

    int idx = 0;
    for(compressed_matrix<double>::const_iterator1 it = H.begin1(); it != H.end1(); ++it)
    {
      for(compressed_matrix<double>::const_iterator2 it2 = it.begin(); it2 != it.end(); ++it2)
      {
        iRow[idx] = it2.index1();
        jCol[idx] = it2.index2();
        idx++;
      }
    }
    assert(idx == nele_hess);
    LOG_DBG3(feasPP) << "FSP:eh structure: jCol=" << StdVector<int>::ToString(nele_hess, jCol);
    LOG_DBG3(feasPP) << "FSP:eh structure: iRow=" << StdVector<int>::ToString(nele_hess, iRow);
  }
  else
  {
    assert(iRow == NULL && jCol == NULL && x != NULL && lambda != NULL);

    if(feas_pp->non_approx_constraints && values != NULL)
    {
       feas_pp->optimization->GetDesign()->ReadDesignFromExtern(x);
       feas_pp->optimization->GetDesign()->Reset(DesignElement::CONSTRAINT_GRADIENT);
    }

    // evaluate the Hessian of the Lagrangian -> IPOPT docu (9)
    // obj_factor * H(obj) + sum_i lambda_i*H(g_i)

    // reset the values of H. We must not use H.clear() as this removes the pattern and for the second-order derivative tests, therefore overwrite
    for(compressed_matrix<double>::const_iterator1 it = H.begin1(); it != H.end1(); ++it)
      for(compressed_matrix<double>::const_iterator2 it2 = it.begin(); it2 != it.end(); ++it2)
        H(it2.index1(), it2.index2()) = 0.0;
    assert((int) H.nnz() == nele_hess); // still true?

    // objective
    MMAApproximation* f = feas_pp->obj;
    const Matrix<unsigned int>& S = f->hess_pattern;
    assert(S.GetNumRows() == 0 || S.GetNumCols() == 2);
    tmp_.Resize(S.GetNumRows());
    f->Evaluate(x, f->HESSIAN, &tmp_);
    for(unsigned int e=0, en = S.GetNumRows(); e < en; e++)
    {
      H(S(e,0), S(e,1)) = obj_factor * tmp_[e];
      LOG_DBG3(feasPP) << "FSP:eh H[" << S(e,0) << "][" << S(e,0) << "] = " << obj_factor << " * " << tmp_[e] << " -> " << H(S(e,0), S(e,1));
    }

    // constraints
    for(int g = 0; g < m; g++)
    {
      f = feas_pp->constr[g];
      Matrix<unsigned int>& S = f->hess_pattern; // it is necessary to shadow S, otherwise we change the pattern of the objective
      tmp_.Resize(S.GetNumRows());
      f->Evaluate(x, f->HESSIAN, &tmp_);
      for(unsigned int e=0, en = S.GetNumRows(); e < en; e++)
      {
        H(S(e,0), S(e,1)) += lambda[g] * tmp_[e];
        LOG_DBG3(feasPP) << "FSP:eh g=" << f->ToString(true) << " H[" << S(e,0) << "][" << S(e,0) << "] += " << lambda[g] << " * " << tmp_[e] << " -> " << H(S(e,0), S(e,1));
      }
    }
    feas_pp->optimization->constraints.view->Done(); // reset local constraint to global mode

    assert((int) H.nnz() == nele_hess);

    int i = 0;
    for(compressed_matrix<double>::const_iterator1 it = H.begin1(); it != H.end1(); ++it)
      for(compressed_matrix<double>::const_iterator2 it2 = it.begin(); it2 != it.end(); ++it2)
        values[i++] = *it2;

    assert(i == (int) H.nnz()); // H.value_data() can grow to large ?!

    LOG_DBG3(feasPP) << "FSP:eh of=" << obj_factor << " lambda=" << StdVector<double>::ToString(m, lambda) << " x=" << StdVector<double>::ToString(n, x) << " values=" << StdVector<double>::ToString(nele_hess, values);
  }

  return true;
}



void FeasSubProblem::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* ptr_lambda,
                              Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq)
{
  LOG_TRACE(feasPP) << "finalize_solution: status = " << status << "; n = " << n << "; m = " << m
                    << " obj_value = " << obj_value << " x_avg = " << Average(x, n) << " x_std_dev = "
                    << StandardDeviation(x, n) << " z_l_avg = " << Average(z_L, n) << " z_l_std_dev = "
                    << StandardDeviation(z_L, n) << " z_u_avg = " << Average(z_L, n) << " z_u_std_dev = "
                    << StandardDeviation(z_U, n);


  x_final.Resize(n);
  x_final.Fill(x, n);

  LOG_DBG3(feasPP) << "FSP:fs x=[" << x_final.ToString() << "]";

  z_l_.Resize(n);
  z_l_.Import(z_L, n);

  z_u_.Resize(n);
  z_u_.Import(z_U, n);

  // store original lambda for warmstart
  orig_lambda.Resize(m);
  orig_lambda.Fill(ptr_lambda, m);
  LOG_DBG3(feasPP) << "FSP:fs orig_lambda=" << orig_lambda.ToString();

  // calculate transformed lambda according to feasibility paper
  feas_pp->optimization->GetDesign()->ReadDesignFromExtern(x_final.GetPointer());
  lambda.Resize(m);
  for(int i = 0; i < m; i++)
    lambda[i] = feas_pp->constr[i]->TransformMultiplyer(orig_lambda[i]);


  LOG_DBG3(feasPP) << "FSP:fs lambda=" << lambda.ToString();

}

