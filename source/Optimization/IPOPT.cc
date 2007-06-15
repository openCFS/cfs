#include "Optimization/IPOPT.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Condition.hh"
#include "General/exception.hh"
#include "Utils/StdVector.hh"
#include "Utils/VecStat.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "IpSolveStatistics.hpp"

#include "boost/lexical_cast.hpp"

using namespace CoupledField;

// declare class specific logging stream
DECLARE_LOG(ipopt)
DEFINE_LOG(ipopt, "ipopt")


IPOPT::IPOPT(Optimization* optimization, ParamNode* pn)
{
  LOG_TRACE(ipopt) << "Initialize IPOPT";
  this->optimization_ = optimization;
  
  // smart pointer!
  app = new IpoptApplication();
  // Initialize the IpoptApplication and process the options
  app->Initialize();

  app->Options()->SetIntegerValue("max_iter", optimization_->GetMaxIterations());
  
  // up to now we don't have hessian  
  app->Options()->SetStringValue("hessian_approximation", "limited-memory");

  // save the obj scaling for get_scaling_parameters() -> is 1.0 if not set
  obj_scaling_factor_ = pn != NULL && pn->Has("option", "name", "obj_scaling_factor") ?
    pn->Get("option", "name", "obj_scaling_factor")->Get("value")->AsDouble() : 1.0;
  // handle maximation  
  obj_scaling_factor_ *= optimization_->GetObjective()->task == Optimization::MAXIMIZE ? -1.0 : 1.0;  
  
  // do scaling via get_scaling_parameters() only if we do constraint scaling 
  // otherwise IPOPT is strange
  bool g_scale = false;
  for(unsigned int i = 0; i < optimization->constraints.GetSize(); i++)
    if(optimization->constraints[i].scaling != 1.0) g_scale = true;    

  if(g_scale)
  {
    // here we also set the obj_scaling_factor in get_scaling_parameters()
    app->Options()->SetStringValue("nlp_scaling_method", "user-scaling");
  }
  else
  {
    // no constraint scaling - hence better no get_scaling_parameters() to be called
    if(obj_scaling_factor_ != 1.0)
      app->Options()->SetNumericValue("obj_scaling_factor", obj_scaling_factor_);
  }
  
  // check for optional paramters
  if(pn != NULL)
  {
    StdVector<ParamNode*> list = pn->GetList("option", "type", "string");
    for(unsigned int i = 0; i < list.GetSize(); i++)
      app->Options()->SetStringValue(list[i]->Get("name")->AsString(), list[i]->Get("value")->AsString());

    list = pn->GetList("option", "type", "integer");
    for(unsigned int i = 0; i < list.GetSize(); i++)
      app->Options()->SetIntegerValue(list[i]->Get("name")->AsString(), list[i]->Get("value")->AsInt());

    list = pn->GetList("option", "type", "real");
    for(unsigned int i = 0; i < list.GetSize(); i++)
    {
      // do not set obj_scaling_factor -> it is set via get_scaling_parameters() or before with maximation factor
      if(list[i]->Get("name")->AsString() != "obj_scaling_factor")
        app->Options()->SetNumericValue(list[i]->Get("name")->AsString(), list[i]->Get("value")->AsDouble());
    }
  }
  
  // initialize the VecStat where we store the initial and current gradients statistics
  // the 0-element is for the function, then for the gradients
  initial_.Resize(optimization->constraints.GetSize()+1); 
  current_.Resize(optimization->constraints.GetSize()+1);
  
}

IPOPT::~IPOPT()
{}


void IPOPT::SolveProblem()
{
  ApplicationReturnStatus status = app->OptimizeTNLP(this);

  if (status == Solve_Succeeded) {
    // Retrieve some statistics about the solve
    Index iter_count = app->Statistics()->IterationCount();
    printf("\n\n*** The problem solved in %d iterations!\n", iter_count);

    Number final_obj = app->Statistics()->FinalObjective();
    printf("\n\n*** The final value of the objective function is %e.\n", final_obj);
    return;
  }

  switch(status)
  {
    case NonIpopt_Exception_Thrown:
         throw Exception("IPOPT stopped due to non-IPOPT exception. Try again with '-f'.");
         
    case Restoration_Failed:
         throw Exception("IPOPT stopped with 'Restautation failed'");  
         
    case Insufficient_Memory:
         throw Exception("IPOPT reports insufficient memory.");
         
     default:
        // positive is warning
        // Maximum_Iterations_Exceeded == -1
       if(status < Maximum_Iterations_Exceeded) 
         throw Exception("IPOPT reported an error during optimization "
                          + boost::lexical_cast<std::string>(status));
       // else is no bad error and exits this void method :)  
  }
}


bool IPOPT::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  n = optimization_->GetDesign()->data.GetSize();

  // arbitrary constraints ,
  m = optimization_->constraints.GetSize();

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

  const StdVector<DesignElement>* data = &optimization_->GetDesign()->data;
  // we know only simp up to now!
  for(int i = 0; i < n; i++)
  {
    x_l[i] = (*data)[i].GetLowerBound();
    x_u[i] = (*data)[i].GetUpperBound();
  }

  if(m != (int) optimization_->constraints.GetSize()) 
    throw Exception("mixed up number of constraints");

  if(n != (int) optimization_->GetDesign()->data.GetSize())
    throw Exception("mixed up design size");  
    
  for(int i = 0; i < m; i++)
  {
    Condition* g = &optimization_->constraints[i];
    g_l[i] = g_u[i] = g->value;
    // see ipopt documentation 
    if(g->GetType() == Condition::LOWER_BOUND) g_u[i] = 1e19;
    if(g->GetType() == Condition::UPPER_BOUND) g_l[i] = -1e19;
    LOG_TRACE2(ipopt) << "set bounds: g[" << i << "]_u=" << g_u[i] << " g[" << i << "]_l="
                      << g_l[i] << ": g: " << g->ToString() << " type: " 
                      << Condition::type.ToString(g->GetType()) << " value: " << g->value
                      << " parameter: " << g->parameter; 
  }
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
  optimization_->GetDesign()->WriteDesignToExtern(x);

  return true;
}

bool IPOPT::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  LOG_DBG(ipopt) << "eval_grad_f: nex_x = " << new_x << " NeedObjectiveEval = " 
                 << optimization_->NeedObjectiveEval(x) << " x_avg = " << Average(x, n) 
                 << " std_dev = " << StandardDeviation(x, n);;
  
  // return the value of the objective function.
  try
  {
    optimization_->GetDesign()->ReadDesignFromExtern(x);
    optimization_->SolveStateProblem();
    obj_value = optimization_->CalcObjective();
  }
  catch(Exception& e)
  {
    std::cerr << "CFS exception occured within a call from IPOPT:" << std::endl << e.what() << std::endl;
    throw IpoptException(e.what(), __FILE__, __LINE__); 
  }
  
  LOG_TRACE2(ipopt) << "eval_f: new_x = " << new_x << " x_avg = " << Average(x, n)
                    << " x_std_dev = " << StandardDeviation(x, n) << " -> obj_value = " << obj_value;  
  return true;
}

bool IPOPT::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  // return the gradient of the objective function grad_{x} f(x)
  LOG_TRACE2(ipopt) << "eval_grad_f-1: n = " << n << "; new_x = " << new_x
                    << " x_avg = " << Average(x, n) << " std_dev = " 
                    << StandardDeviation(x, n);;

  LOG_DBG(ipopt) << "eval_grad_f: nex_x = " << new_x << " NeedObjectiveEval = " 
                 << optimization_->NeedObjectiveEval(x) << " x_avg = " << Average(x, n) 
                 << " std_dev = " << StandardDeviation(x, n);;

  // check if we have to eval f first!
  if(optimization_->NeedObjectiveEval(x))
  { 
      LOG_TRACE2(ipopt) << "eval_grad_f-1.1:IPOPT design space is not equal! -> eval_f";
      Number obj_value;
      eval_f(n, x, new_x, obj_value);
  }    

  optimization_->CalcObjectiveGradient(grad_f);
  LOG_TRACE2(ipopt) << "eval_grad_f-2: -> avg = " << Average(grad_f, n) << " std_dev = " 
                    << StandardDeviation(grad_f, n);
  
  // store the current gradient statistics
  if(!initial_[0].IsInitialized())
  { 
    initial_[0].Calc(grad_f, n);
    std::cout << "Initial grad_f is " << initial_[0].ToString() << std::endl;
  }
  current_[0].Calc(grad_f, n);
                       
  return true;
}

bool IPOPT::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  LOG_TRACE2(ipopt) << "eval_g: n = " << n << "; new_x = " << new_x << "; m = " << m;

  // we overwrite the design space, but we do this all the time - especially befor eval_f
  optimization_->GetDesign()->ReadDesignFromExtern(x); 
     
  // iterate over all constraints
  for(int i = 0; i < m; i++) 
  {   
     double value =  optimization_->CalcConstraint(&optimization_->constraints[i]);
     g[i] = value;
     LOG_TRACE2(ipopt) << "eval_g-" << i << ": new_x = " << new_x << " x_avg = "  
                       << Average(x, n) << " std_dev = " << StandardDeviation(x, n)
                       << " -> " << optimization_->constraints[i].ToString() << " = " 
                       << value;
  }

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
    // note, that we have dense gradient pointers!
    // iterate over the gradients
    for(int c = 0; c < m; c++)
    {
      double* ptr = values + (c*n);
      optimization_->CalcConstraintGradient(&optimization_->constraints[c], ptr);
      LOG_TRACE2(ipopt) << "eval_jac_g-" << c << ": new_x=" << new_x << " x_avg= " 
                        << Average(x, n) << " std_dev = " << StandardDeviation(x, n)
                        << " -> avg = " << Average(ptr, n) << " std_dev = " << StandardDeviation(ptr, n);
 
      // pos 0 is for objective, then the constraints
      if(!initial_[c+1].IsInitialized()) 
      {
        initial_[c+1].Calc(ptr, n);
        std::cout << "Initial jac_g[" << c << "] -> " << initial_[c+1].ToString() << std::endl;
      }
      current_[c+1].Calc(ptr, n);
    }
  }
  return true;
}

void IPOPT::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number obj_value)
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
  std::cout << "Objective gradient -> " << current_[0].ToString() << std::endl;
  for(int i = 0; i < m; i++)
    std::cout << "jac_g[" << i << "] -> " << current_[i+1].ToString() << std::endl;
  
  // we do not sture final x as we stored the history already ?! 
}

bool IPOPT::intermediate_callback(AlgorithmMode mode,Index iter, Number obj_value,
     Number inf_pr, Number inf_du, Number mu, Number d_norm,Number regularization_size,
     Number alpha_du, Number alpha_pr, Index ls_trials, const IpoptData* ip_data,
     IpoptCalculatedQuantities* ip_cq)
{
  LOG_TRACE2(ipopt) << "intermediate_callback: mode = " << mode << "; iter = " << iter 
                    << " obj_value = " << obj_value << "; ls_trials = " << ls_trials;

  optimization_->CommitIteration();
  // break the ipopt calculations - e.g. if our relative change is smaller than given in xml
  return optimization_->IsMinimumReached() ? false : true;
}     

bool IPOPT::get_scaling_parameters(Number& obj_scaling, bool& use_x_scaling, 
              Index n, Number* x_scaling, bool& use_g_scaling, Index m, Number* g_scaling)
{
  // this method is only called if nlp_scaling_method is set to user-scaling!

  // scaling from the xml file
  obj_scaling = obj_scaling_factor_;

  // we don't do x_scaling
  use_x_scaling = false;

  use_g_scaling = false;
  for(int i = 0; i < m; i++)
  {
    Condition* g = &optimization_->constraints[i];
    g_scaling[i] = g->scaling;
    if(g->scaling != 1.0) use_g_scaling = true;    
  }

  // note, that when use_g_scaling = false we have a different IPOPT behaviour than with
  // not calling get_scaling_parameters().
  LOG_TRACE(ipopt) << "get_scaling_parameters: obj_scaling -> " << obj_scaling << " use_g_scaling -> "
                   << use_g_scaling;
                   
  return true;                 
}              

