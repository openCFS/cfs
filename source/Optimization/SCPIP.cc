#include "Optimization/SCPIP.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Condition.hh"
#include "General/exception.hh"
#include "Utils/StdVector.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

using namespace CoupledField;

// declare class specific logging stream
DECLARE_LOG(scpip)
DEFINE_LOG(scpip, "scpip")

SCPIP::SCPIP(Optimization* optimization, ParamNode* optimizer_pn)
 : BaseOptimizer(optimization, optimizer_pn)
{
  LOG_TRACE(scpip) << "Initialize SCPIP";

  // reduce to our actual ParamNode
  ParamNode* scpip_pn =  optimizer_pn->Get(Optimization::optimizer.ToString(Optimization::SCPIP_SOLVER), false);
  
  // init scaling - to be set in SolveProblem() for restarted autoscale
  double manual_scaling = 1.0;
  if(scpip_pn != NULL) scpip_pn->Get<double>("option", "name", "obj_scaling_factor", manual_scaling, "value", false);
  
  //double manual_scaling = scpip_pn != NULL && optimizer_pn->Has("option", "name", "obj_scaling_factor") ?
  //    scpip_pn->Get("option", "name", "obj_scaling_factor")->Get("value")->AsDouble() : 1.0;
  PostInit(manual_scaling); // does autoscale    
  std::cout << objective->ToString() << std::endl;
  
  Initialize();

  SetIntegerValue("max_iter", optimization->GetMaxIterations());
  
  // SCPIP has not automated scaling we might disable with 1.0
  SetNumericValue("obj_scaling_factor", objective->scaling.value);
  
  // check for optional parameters
  if(scpip_pn != NULL)
  {
    StdVector<ParamNode*> list;
    
    list = scpip_pn->GetList("option", "type", "string");
    for(unsigned int i = 0; i < list.GetSize(); i++)
      SetStringValue(list[i]->Get("name")->AsString(), list[i]->Get("value")->AsString());

    list = scpip_pn->GetList("option", "type", "integer");
    for(unsigned int i = 0; i < list.GetSize(); i++)
      SetIntegerValue(list[i]->Get("name")->AsString(), list[i]->Get("value")->AsInt());

    list = scpip_pn->GetList("option", "type", "real");
    for(unsigned int i = 0; i < list.GetSize(); i++)
    {
      // do not set obj_scaling_factor -> it was set before
      if(list[i]->Get("name")->AsString() != "obj_scaling_factor")
        SetNumericValue(list[i]->Get("name")->AsString(), list[i]->Get("value")->AsDouble());
    }
  }
}

SCPIP::~SCPIP()
{}


void SCPIP::SolveProblem()
{
  do
  {
    if(restart_requested)
    {  
      // set the scaling here to handle restart
      objective->CalcAutoscale();
      LOG_TRACE(scpip) << "Restart SCPIP with scaling " << objective->scaling.value;
      std::cout << "Restart SCPIP with scaling " << objective->scaling.value << std::endl;
      
      // adjust the number of iterations
      int max_iter = std::max(optimization->GetMaxIterations() - optimization->GetCurrentIteration(), 0);
      SetIntegerValue("max_iter", max_iter);
      LOG_TRACE2(scpip) << "set max_iter to " << max_iter << " (" << optimization->GetMaxIterations()
                        << " - " << optimization->GetCurrentIteration() << ")";
      if(max_iter == 0)
      {
        std::cout << std::endl << "maximum number of iterations is 0 " << std::endl;
        return;
      }
    }
    // SCPIP has not automated scaling we might disable with 1.0
    SetNumericValue("obj_scaling_factor", objective->scaling.value);
    
    // call the base solver !!!!!!!!!!!
    int status = SCPIPBase::SolveProblem();

    switch(status)
    {
    case Solve_Succeeded:
      std::cout << std::endl << "The problem solved in " 
      << Statistics()->IterationCount() << " iterations!" << std::endl;
      std::cout << "The final value of the objective function is "
      << Statistics()->FinalObjective() << std::endl;
      break;          
      
    case Subproblem_Max_Iter:
    case Maximum_Iterations_Exceeded:
    case User_Requested_Stop:
      std::cout << std::endl << ToString(status) << std::endl;
      break; 
    
    case Gradients_Return_False:
      if(restart_requested) 
        break; // else throw exception
      
    default:
      throw Exception(ToString(status));
    }
  }
  while(restart_requested);
}

bool  SCPIP::get_nlp_info(int& n, int& m, int& nnz_jac_g)
{
  n = optimization->GetDesign()->data.GetSize();

  // arbitrary constraints ,
  m = optimization->constraints.GetSize();

  // up to now we have only dense constraint gradients. 
  // In practice one could make the non-matching part spare or
  // combine the sparse parts of constraints ... but this is future!
  nnz_jac_g = n * m; 

  return true;
}

bool SCPIP::get_bounds_info(int n, double* x_l, double* x_u,
                            int m, double* g_l, double* g_u)
{
  const StdVector<DesignElement>* data = &optimization->GetDesign()->data;

  for(int i = 0; i < n; i++)
  {
    x_l[i] = (*data)[i].GetLowerBound();
    x_u[i] = (*data)[i].GetUpperBound();
  }

  assert(m == (int) optimization->constraints.GetSize());
  assert(n == (int) optimization->GetDesign()->data.GetSize());
    
  // normalization to =0 and <=0 constraints is done SCPIPBase   
    
  for(int i = 0; i < m; i++)
  {
    Condition* g = &optimization->constraints[i];
    // handle as in IPOPT 
    g_l[i] = g_u[i] = g->value;    
    if(g->GetType() == Condition::LOWER_BOUND) g_u[i] = 1e19;
    if(g->GetType() == Condition::UPPER_BOUND) g_l[i] = -1e19;
  }
  return true;
}

bool SCPIP::get_starting_point(int n, double* x)
{
  // we initialize x in bounds, in the upper right quadrant
  optimization->GetDesign()->WriteDesignToExtern(x);

  return true;
}

bool SCPIP::eval_f(int n, const double* x, double& obj_value)
{
  obj_value = EvalObjective(n, x);

  return true;
}

bool SCPIP::eval_grad_f(int n, const double* x, double* grad_f)
{
  // restart_requested handled in intermediate_callback
  return EvalGradObjective(n, x, grad_f); 
}

bool SCPIP::eval_g(int n, const double* x, int m, double* g)
{
  // we overwrite the design space, but we do this all the time - especially before eval_f
  optimization->GetDesign()->ReadDesignFromExtern(x); 
     
  // iterate over all constraints
  for(int i = 0; i < m; i++) 
     g[i] = optimization->CalcConstraint(&optimization->constraints[i]);

  return true;
}

bool SCPIP::eval_jac_g(int n, const double* x, int m, int nele_jac, double* values)
{
  // the gradients are dense in SCPIPBase
  assert(values != NULL);

  // note, that we have dense gradient pointers!
  // iterate over the gradients
  for(int c = 0; c < m; c++)
  {
    double* ptr = values + (c*n);
    optimization->CalcConstraintGradient(&optimization->constraints[c], ptr);
    
    LOG_TRACE2(scpip) << "eval_jac_g-" << c << ": x_avg= " << Average(x, n) 
                      << " std_dev = " << StandardDeviation(x, n)
                      << " -> avg = " << Average(ptr, n) << " std_dev = " 
                      << StandardDeviation(ptr, n);
  }
  return true;
}


void SCPIP::finalize_solution(int status, int n, const double* x, const double* z_L, 
                             const double* z_U, int m, const double* g, 
                             const double* lambda, double obj_value)
{
  LOG_TRACE2(scpip) << "finalize_solution: status = " << status << "; n = " << n << "; m = " << m 
                    << " obj_value = " << obj_value << " x_avg = " << Average(x, n) << " x_std_dev = " 
                    << StandardDeviation(x, n) << " z_l_avg = " << Average(z_L, n) << " z_l_std_dev = "
                    << StandardDeviation(z_L, n) << " z_u_avg = " << Average(z_L, n) << " z_u_std_dev = "
                    << StandardDeviation(z_U, n) << "restart_requested = " << restart_requested;

  // save this iteration, otherwise it might be lost
  optimization->CommitIteration();  
  
  if(restart_requested) return;
  
  std::cout << "SCPIP finished: f=" << obj_value;
  for(int i = 0; i < m; i++)
    std::cout << " + " << lambda[i] << "*" << g[i];
  std::cout << std::endl;  
}
 
bool SCPIP::intermediate_callback(int iter, bool next_iter)
{
  if(iter == 0 || !next_iter) return true;

  optimization->CommitIteration();

  // break the optimization - e.g. if our relative change is smaller than given in xml
  // or we have to do a restart
  LOG_TRACE2(scpip) << "intermediate_callback iter=" << iter << " next_iter=" << next_iter
                    << " restart_requested=" << restart_requested << " min reached="
                    << optimization->IsMinimumReached();
  
  return (restart_requested || optimization->IsMinimumReached()) ? false : true; 
}     

bool SCPIP::get_scaling_parameters(double& obj_scaling, bool& use_g_scaling, int m, double* g_scaling)
{
  // this method is only called if nlp_scaling_method is set to user-scaling!
  obj_scaling = objective->scaling.value;

  use_g_scaling = false;
  for(int i = 0; i < m; i++)
  {
    Condition* g = &optimization->constraints[i];
    g_scaling[i] = g->scaling;
    if(g->scaling != 1.0) use_g_scaling = true;    
  }
  
  LOG_TRACE(scpip) << "get_scaling_parameters: obj_scaling -> " << obj_scaling << " use_g_scaling -> "
                   << use_g_scaling;
                   
  return true;                 
}              
