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

SCPIP::SCPIP(Optimization* optimization, ParamNode* pn)
{
  LOG_TRACE(scpip) << "Initialize SCPIP";
  this->optimization_ = optimization;
  
  Initialize();

  SetIntegerValue("max_iter", optimization_->GetMaxIterations());
  
  // the maximation is done via scaling. We scale with -1 if "obj_scaling_factor" 
  // is not given. If it is given we switch sign
  // save the obj scaling for get_scaling_parameters();
  obj_scaling_factor_ = pn != NULL && pn->Has("option", "name", "obj_scaling_factor") ?
    pn->Get("option", "name", "obj_scaling_factor")->Get("value")->AsDouble() : 1.0;
  // handle maximation  
  obj_scaling_factor_ *= optimization_->GetObjective()->task == Optimization::MAXIMIZE ? -1.0 : 1.0;  

  // set the scaling if it is not 1.0 -> get_scaling_parameters() will be called
  if(obj_scaling_factor_ != 1.0) 
    SetNumericValue("obj_scaling_factor", obj_scaling_factor_);

  // check for optional paramters
  if(pn != NULL)
  {
    StdVector<ParamNode*> list = pn->GetList("option", "type", "string");
    for(unsigned int i = 0; i < list.GetSize(); i++)
      SetStringValue(list[i]->Get("name")->AsString(), list[i]->Get("value")->AsString());

    list = pn->GetList("option", "type", "integer");
    for(unsigned int i = 0; i < list.GetSize(); i++)
      SetIntegerValue(list[i]->Get("name")->AsString(), list[i]->Get("value")->AsInt());

    list = pn->GetList("option", "type", "real");
    for(unsigned int i = 0; i < list.GetSize(); i++)
    {
      // do not set obj_scaling_factor -> it was before with maximation.
      if(list[i]->Get("name")->AsString() != "obj_scaling_factor")
        SetNumericValue(list[i]->Get("name")->AsString(), list[i]->Get("value")->AsDouble());
    }
  }

  // initialize the VecStat where we store the initial and current gradients statistics
  // the 0-element is for the function, then for the gradients
  initial_.Resize(optimization->constraints.GetSize()+1); 
  current_.Resize(optimization->constraints.GetSize()+1);
}

SCPIP::~SCPIP()
{}


void SCPIP::SolveProblem()
{
  // call the base solver
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
    case User_Requested_Stop:
    case Maximum_Iterations_Exceeded:
         std::cout << std::endl << ToString(status) << std::endl;
         break; 
    
    default:
         throw Exception(ToString(status));
  }
}

bool  SCPIP::get_nlp_info(int& n, int& m, int& nnz_jac_g)
{
  n = optimization_->GetDesign()->data.GetSize();

  // arbitrary constraints ,
  m = optimization_->constraints.GetSize();

  // up to now we have only dense constraint gradients. 
  // In practice one could make the non-matching part spare or
  // combine the sparse parts of constraints ... but this is future!
  nnz_jac_g = n * m; 

  return true;
}

bool SCPIP::get_bounds_info(int n, double* x_l, double* x_u,
                            int m, double* g_l, double* g_u)
{
  const StdVector<DesignElement>* data = &optimization_->GetDesign()->data;

  for(int i = 0; i < n; i++)
  {
    x_l[i] = (*data)[i].GetLowerBound();
    x_u[i] = (*data)[i].GetUpperBound();
  }

  assert(m == (int) optimization_->constraints.GetSize());
  assert(n == (int) optimization_->GetDesign()->data.GetSize());
    
  // normalization to =0 and <=0 constraints is done SCPIPBase   
    
  for(int i = 0; i < m; i++)
  {
    Condition* g = &optimization_->constraints[i];
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
  optimization_->GetDesign()->WriteDesignToExtern(x);

  return true;
}

bool SCPIP::eval_f(int n, const double* x, double& obj_value)
{
  LOG_DBG(scpip) << "eval_grad_f: NeedObjectiveEval = " 
                 << optimization_->NeedObjectiveEval(x) << " x_avg = " << Average(x, n) 
                 << " std_dev = " << StandardDeviation(x, n);;
  
  // return the value of the objective function.
  optimization_->GetDesign()->ReadDesignFromExtern(x);
  optimization_->SolveStateProblem();
  obj_value = optimization_->CalcObjective();
  LOG_TRACE2(scpip) << "eval_f: x_avg = " << Average(x, n) << " x_std_dev = " 
                    << StandardDeviation(x, n) << " -> obj_value = " << obj_value;
  return true;
}

bool SCPIP::eval_grad_f(int n, const double* x, double* grad_f)
{
  // check if we have to eval f first! Such that the gradients have the
  // proper solution of the forward problem
  if(optimization_->NeedObjectiveEval(x))
  { 
      LOG_TRACE2(scpip) << "eval_grad_f: design space is not equal! -> eval_f";
      double obj_value;
      eval_f(n, x, obj_value);
  }    

  optimization_->CalcObjectiveGradient(grad_f);
  LOG_TRACE2(scpip) << "eval_grad_f: -> avg = " << Average(grad_f, n) << " std_dev = " 
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

bool SCPIP::eval_g(int n, const double* x, int m, double* g)
{
  // we overwrite the design space, but we do this all the time - especially before eval_f
  optimization_->GetDesign()->ReadDesignFromExtern(x); 
     
  // iterate over all constraints
  for(int i = 0; i < m; i++) 
     g[i] = optimization_->CalcConstraint(&optimization_->constraints[i]);

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
    optimization_->CalcConstraintGradient(&optimization_->constraints[c], ptr);
    
    LOG_TRACE2(scpip) << "eval_g-" << c << ": x_avg= " << Average(x, n) 
                      << " std_dev = " << StandardDeviation(x, n)
                      << " -> avg = " << Average(ptr, n) << " std_dev = " 
                      << StandardDeviation(ptr, n);
                      
    // pos 0 is for objective, then the constraints
    if(!initial_[c+1].IsInitialized()) 
    {
      initial_[c+1].Calc(ptr, n);
      std::cout << "Initial jac_g[" << c << "] -> " << initial_[c+1].ToString() << std::endl;
    }
    current_[c+1].Calc(ptr, n);
  }
  return true;
}


void SCPIP::finalize_solution(int status, int n, const double* x, const double* z_L, 
                             const double* z_U, int m, const double* g, 
                             const double* lambda, double obj_value)
{
  LOG_TRACE(scpip) << "finalize_solution: status = " << status << "; n = " << n << "; m = " << m 
                    << " obj_value = " << obj_value << " x_avg = " << Average(x, n) << " x_std_dev = " 
                    << StandardDeviation(x, n) << " z_l_avg = " << Average(z_L, n) << " z_l_std_dev = "
                    << StandardDeviation(z_L, n) << " z_u_avg = " << Average(z_L, n) << " z_u_std_dev = "
                    << StandardDeviation(z_U, n);

  std::cout << "SCPIP finished: f=" << obj_value;
  for(int i = 0; i < m; i++)
    std::cout << " + " << lambda[i] << "*" << g[i];
  std::cout << std::endl;  
  
  // write the last gradients
  std::cout << "Objective gradient -> " << current_[0].ToString() << std::endl;
  for(int i = 0; i < m; i++)
    std::cout << "jac_g[" << i << "] -> " << current_[i+1].ToString() << std::endl;
  
}
 
bool SCPIP::intermediate_callback(int iter, bool next_iter)
{
  //PrintInfo(std::cout);
  if(iter == 0 || !next_iter) return true;
  
  optimization_->CommitIteration();
  // break the optimization - e.g. if our relative change is smaller than given in xml
  return optimization_->IsMinimumReached() ? false : true;
}     

bool SCPIP::get_scaling_parameters(double& obj_scaling, bool& use_g_scaling, int m, double* g_scaling)
{
  // this method is only called if nlp_scaling_method is set to user-scaling!

  // Note, that this is a ugly copy&paste from IPOPT.cc!

  // scaling from the xml file
  obj_scaling = obj_scaling_factor_;

  // the maximation is done via scaling. We scale with -1 if "obj_scaling_factor" 
  // is not given. If it is given we switch sign
  obj_scaling *= optimization_->GetObjective()->task == Optimization::MAXIMIZE ? -1.0 : 1.0;

  use_g_scaling = false;
  for(int i = 0; i < m; i++)
  {
    Condition* g = &optimization_->constraints[i];
    g_scaling[i] = g->scaling;
    if(g->scaling != 1.0) use_g_scaling = true;    
  }
  
  LOG_TRACE(scpip) << "get_scaling_parameters: obj_scaling -> " << obj_scaling << " use_g_scaling -> "
                   << use_g_scaling;
                   
  return true;                 
}              
