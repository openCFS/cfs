#include "Optimization/Optimizer/IPOPT.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "General/Exception.hh"
#include "Utils/StdVector.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "coin-or/IpSolveStatistics.hpp"

#include "boost/lexical_cast.hpp"

using namespace CoupledField;

// declare class specific logging stream
DEFINE_LOG(ipopt, "ipopt")


IPOPT::IPOPT(Optimization* optimization, BaseOptimizer* base, PtrParamNode pn)
{
  LOG_TRACE(ipopt) << "Initialize IPOPT";
  this->opt_ = optimization;
  this->base_ = base;
  this->optimizer_pn_ = pn->Get(Optimization::optimizer.ToString(Optimization::IPOPT_SOLVER), ParamNode::PASS);
  
  double manual_scaling = optimizer_pn_ != NULL && optimizer_pn_->Has("obj_scaling_factor") ?
      optimizer_pn_->Get("obj_scaling_factor")->Get("value")->As<Double>() : 1.0;
  base->PostInitScale(manual_scaling);
  Init();
}

void IPOPT::Init()
{
  LOG_TRACE(ipopt) << "Init: restart=" << base_->restart_requested;
  
  // smart pointer!
  app = new IpoptApplication();
  // Initialize the IpoptApplication and process the options
  app->Initialize();

  app->Options()->SetIntegerValue("max_iter", opt_->GetMaxIterations() - opt_->GetCurrentIteration());
  LOG_TRACE2(ipopt) << "set max_iter to " << opt_->GetMaxIterations() << " - " 
                    << opt_->GetCurrentIteration();
  
  // up to now we don't have hessian  
  app->Options()->SetStringValue("hessian_approximation", "limited-memory");

  // handle restart case!
  if(base_->restart_requested)
  {
    base_->objective->CalcAutoscale();
  }
  std::cout << base_->objective->ToString() << std::endl;
  
  // do scaling via get_scaling_parameters() only if we do constraint scaling 
  // otherwise IPOPT is strange
  bool g_scale = false;
  for(int i = 0; i < opt_->constraints.view->GetNumberOfActiveConstraints(); i++)
    if(opt_->constraints.view->Get(i)->manual_scaling_value != 1.0) g_scale = true;

  if(g_scale)
  {
    // here we also set the obj_scaling_factor in get_scaling_parameters()
    app->Options()->SetStringValue("nlp_scaling_method", "user-scaling");
  }
  else
  {
    // no constraint scaling - hence better no get_scaling_parameters() to be called
    if(base_->objective->scaling.value != 1.0)
      app->Options()->SetNumericValue("obj_scaling_factor", base_->objective->scaling.value);
  }
    
  // check for optional paramters
  if(optimizer_pn_ != NULL)
  {
    ParamNodeList list;
    list = optimizer_pn_->GetListByVal("option", "type", "string");
    
    for(unsigned int i = 0; i < list.GetSize(); i++)
      app->Options()->SetStringValue(list[i]->Get("name")->As<std::string>(), list[i]->Get("value")->As<std::string>());

    list = optimizer_pn_->GetListByVal("option", "type", "integer");
    for(unsigned int i = 0; i < list.GetSize(); i++)
      app->Options()->SetIntegerValue(list[i]->Get("name")->As<std::string>(), list[i]->Get("value")->As<Integer>());

    list = optimizer_pn_->GetListByVal("option", "type", "real");
    for(unsigned int i = 0; i < list.GetSize(); i++)
    {
      // do not set obj_scaling_factor -> it is set via get_scaling_parameters() or before with maximation factor
      if(list[i]->Get("name")->As<std::string>() != "obj_scaling_factor")
        app->Options()->SetNumericValue(list[i]->Get("name")->As<std::string>(), list[i]->Get("value")->As<Double>());
    }
  }
}

IPOPT::~IPOPT()
{}


void IPOPT::SolveProblem()
{
  ApplicationReturnStatus status = app->OptimizeTNLP(this);

  if(status == Invalid_Number_Detected && base_->restart_requested)
  {
    string msg = "Invalid_Number_Detected in major " + std::to_string(base_->optimization->GetCurrentIteration()) + " do restart";
    base_->info_->Get(ParamNode::SUMMARY, ParamNode::APPEND)->SetValue(msg);
    return;
  }

  PtrParamNode in = opt_->DoStopOptimizationHelper(status >= 0 , "IPOPT status code " + std::to_string(status));

  if (status == Solve_Succeeded)
  {
    // Retrieve some statistics about the solve
    Index iter_count = app->Statistics()->IterationCount();
    Number final_obj = app->Statistics()->FinalObjective();
    
    std::cout << std::endl << "Problem solved in " << iter_count 
              << " iterations, final objective value is " << final_obj << std::endl; 
    return;
  }

  switch(status)
  {
    case NonIpopt_Exception_Thrown:
      opt_->DoStopOptimizationHelper(false,"non IPOPT exception occurred.");
      throw Exception("IPOPT stopped due to non-IPOPT exception. Try again with '-f'.");
         
    case Restoration_Failed:
      opt_->DoStopOptimizationHelper(false,"IPOPT: 'Restoration failed'");
      throw Exception("IPOPT stopped with 'Restoration failed'");
         
    case Insufficient_Memory:
      opt_->DoStopOptimizationHelper(false,"IPOPT: insufficient memory");
      throw Exception("IPOPT reports insufficient memory.");
         
    case Maximum_Iterations_Exceeded:
      opt_->DoStopOptimizationHelper(false,"Maximum iterations exceeded");
      break;
      
    default:
      if(status < 0)
        std::cout << std::endl << "Optimization likely failed: " << in->Get("reason/msg")->As<string>() << std::endl;
      // converged and status already set
      break;
  }
}


bool IPOPT::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  n = opt_->GetDesign()->GetNumberOfVariables();

  // arbitrary constraints ,
  m = opt_->constraints.view->GetNumberOfActiveConstraints();

  // up to now we have only dense constraint gradients. 
  // In practice one could make the non-matching part spare or
  // combine the sparse parts of constraints ... but this is future!
  nnz_jac_g = n * m; 

  // we have no hessian
  nnz_h_lag = 0;

  // We use the standard C index style for row/col entries
  index_style = C_STYLE;

  LOG_TRACE2(ipopt) << ":get_nlp_info n <- " << n << "; m <- " << m << " nnz_jac_g <- " 
                    << nnz_jac_g << " nnz_h_lag <- " << nnz_h_lag << " index_style <- "
                    << (index_style == C_STYLE ? "C_STYLE" : "FORTRAN_STYLE");
  return true;
}

bool IPOPT::get_bounds_info(Index n, Number* x_l, Number* x_u,
                            Index m, Number* g_l, Number* g_u)
{
  LOG_TRACE2(ipopt) << "get_bounds_info: n = " << n << "; m = " << m;
  
  base_->GetBounds(n, x_l, x_u, m, g_l, g_u);

  return true;
}

bool IPOPT::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{
  LOG_TRACE2(ipopt) << "get_starting_point: n = " << n << "; m = " << m;

  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the others if
  // you wish.
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // we initialize x in bounds, in the upper right quadrant
  opt_->GetDesign()->WriteDesignToExtern(x);

  return true;
}

bool IPOPT::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  int old_design = base_->objective->scaling.design_id;
  nObj++;

  // return the value of the objective function.
  try
  {
    obj_value = base_->EvalObjective(n, x, false); // IPOPT can scale by itself!
  }
  catch(Exception& e)
  {
    std::cerr << "CFS exception occured within a call from IPOPT:" << std::endl << e.what() << std::endl;
    throw IpoptException(e.what(), __FILE__, __LINE__); 
  }
    
  LOG_TRACE2(ipopt) << "eval_f: new_x = " << new_x << " design=" << base_->objective->scaling.design_id
                    << " needed_eval=" << (base_->objective->scaling.design_id != old_design)
                    << " -> obj_value = " << obj_value;  
  return true;
}

bool IPOPT::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  int old_design = base_->objective->scaling.design_id;
  // return the gradient of the objective function grad_{x} f(x)

  // TODO: do better :)
  StdVector<double> tmp(n);

  bool ok = base_->EvalGradObjective(n, x, false, tmp);

  for(unsigned int i = 0; i < tmp.GetSize(); i++)
    grad_f[i] = tmp[i];

  LOG_TRACE2(ipopt) << "eval_grad_f: new_x = " << new_x << " design=" << base_->objective->scaling.design_id
                    << " needed_eval=" << (base_->objective->scaling.design_id != old_design)
                    << " good=" << ok << " -> avg = " << Average(grad_f, n) << " std_dev = " 
                    << StandardDeviation(grad_f, n);
  return ok;
}

bool IPOPT::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  LOG_TRACE2(ipopt) << "eval_g: n = " << n << "; new_x = " << new_x << "; m = " << m;

  // we overwrite the design space, but we do this all the time - especially before eval_f
  opt_->GetDesign()->ReadDesignFromExtern(x);
  
  base_->EvalConstraints(n, x, m, false, g, false);

  return true;
}

bool IPOPT::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{
  if (values == NULL) 
  {
    // return the structure of the jacobian of the constraints
    
    // note that we are dense and have for every constraint n elements.
    // These are 0.0 if the design does not match!

    for(int row = 0; row < m; row++)
    {
      for(int col = 0; col < n; col++)
      {
        // the jacobian array/matrix is a long vector with n*m elements
        // see the IPOPT example in the docu to understand!
        int index = row * n + col;
        iRow[index] = row;
        jCol[index] = col;
      }
    }
  }
  // do evaluation
  else
  {
    StdVector<double> tmp(nele_jac);

    base_->EvalGradConstraints(n, x, m, nele_jac, false, false, tmp);

    for(unsigned int i = 0; i < tmp.GetSize(); i++)
      values[i] = tmp[i];
  }
  return true;
}

void IPOPT::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number obj_value, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq)
{
  LOG_TRACE(ipopt) << "finalize_solution: status = " << status << "; n = " << n << "; m = " << m 
                    << " obj_value = " << obj_value << " x_avg = " << Average(x, n) << " x_std_dev = " 
                    << StandardDeviation(x, n) << " z_l_avg = " << Average(z_L, n) << " z_l_std_dev = "
                    << StandardDeviation(z_L, n) << " z_u_avg = " << Average(z_L, n) << " z_u_std_dev = "
                    << StandardDeviation(z_U, n);

  // write a lagrange multiplier formulation
  std::cout << "IPOPT finished: f=" << obj_value;
  for(int i = 0; i < m; i++)
    std::cout << " + " << lambda[i] << "*" << g[i];
  std::cout << std::endl;  

  // write the last gradients
 // std::cout << "Objective gradient -> " << current_[0].ToString() << std::endl;
  //for(int i = 0; i < m; i++)
  //  std::cout << "jac_g[" << i << "] -> " << current_[i+1].ToString() << std::endl;
  
  // we do not sture final x as we stored the history already ?! 
}

bool IPOPT::intermediate_callback(AlgorithmMode mode,Index iter, Number obj_value,
     Number inf_pr, Number inf_du, Number mu, Number d_norm,Number regularization_size,
     Number alpha_du, Number alpha_pr, Index ls_trials, const IpoptData* ip_data,
     IpoptCalculatedQuantities* ip_cq)
{
  LOG_TRACE2(ipopt) << "intermediate_callback: mode = " << mode << "; iter = " << iter 
                    << " obj_value = " << obj_value << "; ls_trials = " << ls_trials;

  base_->CommitIteration();
  // break the ipopt calculations - e.g. if our relative change is smaller than given in xml
  return opt_->DoStopOptimization() ? false : true;
}     

bool IPOPT::get_scaling_parameters(Number& obj_scaling, bool& use_x_scaling, 
              Index n, Number* x_scaling, bool& use_g_scaling, Index m, Number* g_scaling)
{
  // this method is only called if nlp_scaling_method is set to user-scaling!

  // scaling from the xml file
  obj_scaling = base_->objective->scaling.value;

  // we don't do x_scaling
  use_x_scaling = false;

  use_g_scaling = false; // preliminary

  for(int i = 0; i < m; i++)
  {
    Condition* g = opt_->constraints.view->Get(i);
    g_scaling[i] = g->DoObjectiveScaling() ? obj_scaling : g->manual_scaling_value;
    if(g_scaling[i] != 1.0) use_g_scaling = true;

    LOG_TRACE(ipopt) << "get_scaling_parameters: g=" << g->type.ToString(g->GetType()) << " scaling=" << g_scaling[i];
  }
  opt_->constraints.view->Done();

  // note, that when use_g_scaling = false we have a different IPOPT behavior than with
  // not calling get_scaling_parameters().
                   
  return true;                 
}              

