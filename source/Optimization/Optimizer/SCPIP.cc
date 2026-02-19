#include <stddef.h>
#include <algorithm>
#include <iostream>
#include <string>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/Exception.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Optimizer/SCPIP.hh"
#include "Utils/StdVector.hh"
#include "Utils/tools.hh"

using namespace CoupledField;

// declare class specific logging stream
DEFINE_LOG(scpip, "scpip")

SCPIP::SCPIP(Optimization* optimization, PtrParamNode optimizer_pn, Optimization::Optimizer type) :
 BaseOptimizer(optimization, optimizer_pn, type)
{
  LOG_TRACE(scpip) << "Initialize SCPIP";

  SetIntegerValue("max_iter", optimization->GetMaxIterations());
  
  // check for optional parameters
  if(this_opt_pn_ != NULL)
  {
    ParamNodeList list;

    list = this_opt_pn_->GetListByVal("option", "type", "string");
    for(unsigned int i = 0; i < list.GetSize(); i++)
      SetStringValue(list[i]->Get("name")->As<std::string>(), list[i]->Get("value")->As<std::string>());

    list = this_opt_pn_->GetListByVal("option", "type", "integer");
    for(unsigned int i = 0; i < list.GetSize(); i++)
      SetIntegerValue(list[i]->Get("name")->As<std::string>(), list[i]->Get("value")->As<int>());

    list = this_opt_pn_->GetListByVal("option", "type", "real");
    for(unsigned int i = 0; i < list.GetSize(); i++)
    {
      SetNumericValue(list[i]->Get("name")->As<std::string>(), list[i]->Get("value")->As<double>());
    }
  }
}

SCPIP::~SCPIP()
{}

void SCPIP::PostInit()
{

  BaseOptimizer::PostInit();
  Initialize();
  opt_timer->Stop();
  PostInitScale(1.0); // does autoscale

}

void SCPIP::ToInfo(PtrParamNode pn)
{
  BaseOptimizer::ToInfo(pn);

  PtrParamNode pn_ = pn->Get("icntl");
  for(int i = 0; i < 13; i++)
    pn_->Get("i" + lexical_cast<std::string>(i + 1))->SetValue(icntl[i]);

  pn_ = pn->Get("rcntl");
  for(int i = 0; i < 7; i++)
    pn_->Get("r" + lexical_cast<std::string>(i + 1))->SetValue(rcntl[i]);
}

void SCPIP::SolveProblem()
{
  // if we did scale, we can easily commit a calculated initial (iter-0) configuration
  // otherwise we also try to make this
  assert(opt_timer->IsRunning());

  if(objective->DoAutoscale())
     CommitIteration();
  
  do
  {
    if(restart_requested)
    {  
      // set the scaling here to handle restart
      objective->CalcAutoscale();
      LOG_TRACE(scpip) << "Restart SCPIP with scaling " << objective->scaling.value;
      std::cout << "Restart SCPIP with scaling " << objective->scaling.value << std::endl;
      PtrParamNode in = optimization->optInfoNode->Get(ParamNode::PROCESS)->Get("iteration");
      in->Get("rescale")->SetValue(objective->scaling.value);
      
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
    
    // call the base solver, which is SCPIPBase::solve_problem() or FeasSCP::solve_problem()
    int status = solve_problem();
    
    PtrParamNode in = optimization->optInfoNode->Get(ParamNode::SUMMARY)->Get("break");
    
    switch(status)
    {
    case Solve_Succeeded:
      std::cout << std::endl << "Problem solved in " << Statistics()->IterationCount() 
                << " iterations, final objective value is " << Statistics()->FinalObjective() << std::endl; 
      std::cout << "The final value of the objective function is "
                << Statistics()->FinalObjective() << std::endl;
      in->Get("converged")->SetValue("yes");
      break;          

    case User_Requested_Stop:
      // "break" already set!
      assert(in->GetChildren().GetSize() > 0);
      break;

    case Subproblem_Max_Iter:
    case LineSearch_Max_Iter:
    case Maximum_Iterations_Exceeded:
    case Infeasible:
    case Gradients_Return_False:
      optimization->DoStopOptimizationHelper(false, "SCPIP: " + ToString(status));
      break;
      
    default:
      optimization->DoStopOptimizationHelper(false, "SCPIP: " + ToString(status));
      throw Exception(ToString(status));
    }
    
    if(!restart_requested)
      optimization->optInfoNode->Get(ParamNode::SUMMARY)->Get("problem")->SetValue(ToString(status));
  }
  while(restart_requested);
}


void SCPIP::SetConstraintSparsityPattern()
{
  assert(iern.GetSize() != 0 && iecn.GetSize() != 0 && iederv.GetSize() != 0);
  assert(eqrn.GetSize() != 0 && eqcn.GetSize() != 0 && eqcoef.GetSize() != 0);
  assert(m >= 0 && mie >= 0 && meq >= 0);

  // SCPIPBase has the constraints slit into equality and inequality constraints
  // beside this the ordering is not changed.
  int ie = 0; // summed up total inequality constraint
  int ie_active = 0; // summed up active inequality constraint
  (void) ie_active; // "use" the variable to avoid unused warning as we only use it in an assert
  int eq = 0;

  int ie_nnz = 0; // counter
  int eq_nnz = 0;

  assert(m == optimization->constraints.view->GetNumberOfActiveConstraints());

  ie_idx.Reserve(m);
  eq_idx.Reserve(m);

  for(int c = 0; c < m; c++)
  {
    Condition* g = optimization->constraints.view->Get(c);
    if(g->GetBound() != Condition::EQUAL)
    {
      assert(active[ie] == 1 || active[ie] == 0);
      if(active[ie] == 1)
      {
        StdVector<unsigned int>& pattern = g->GetSparsityPattern();
        for(unsigned int e = 0; e < pattern.GetSize(); e++)
        {
          iern[ie_nnz] = pattern[e] + 1; // fortran!
          iecn[ie_nnz] = ie +1; // fortran
          ie_nnz++;
        }
        ie_active++;
      }
      ie++;
      ie_idx.Push_back(c);
    }
    else // equality constraint
    {
      StdVector<unsigned int>& pattern = g->GetSparsityPattern();
      for(unsigned int e = 0; e < pattern.GetSize(); e++)
      {
        eqrn[eq_nnz] = pattern[e] + 1; // fortran!
        eqcn[eq_nnz] = eq +1; // fortran
        eq_nnz++;
      }
      eq++;
      eq_idx.Push_back(c);
    }
    LOG_DBG(scpip) << "SCSP g=" << g->ToString() << " ie=" << ie_nnz << " eq=" << eq_nnz;
  }

  optimization->constraints.view->Done(); // mandatory after traversing the view

  LOG_DBG(scpip) << "SCSP ie=" << ie_nnz << " eq=" << eq_nnz << " ieleng=" << ieleng << " eqleng=" << eqleng;

  assert(ie_nnz + eq_nnz == ieleng + eqleng);
  assert(mactiv == ie_active + eq);
  assert(ie + eq <= m);
}

bool  SCPIP::get_nlp_info(int& n, int& m, int& nnz_jac_g)
{
  n = optimization->GetDesign()->GetNumberOfVariables();

  // arbitrary constraints ,
  m = optimization->constraints.view->GetNumberOfActiveConstraints();

  // SCPIPBase will do this again, we use jac_g_size as "buffer"
  jac_g_size.Resize(m);
  nnz_jac_g = get_sparsity_pattern_size(m, m > 0 ? jac_g_size.GetPointer() : NULL);

  return true;
}

int SCPIP::get_sparsity_pattern_size(int m, int* jac_g_dim)
{
  assert(m == optimization->constraints.view->GetNumberOfActiveConstraints());

  int nnz = 0;
  for(int i = 0; i < m; i++)
  {
    Condition* g = optimization->constraints.view->Get(i);
    int local = g->GetSparsityPatternSize();
    jac_g_dim[i] = local;
    nnz += local;
    LOG_DBG(scpip) << "gsps: i=" << i << " g=" << g->ToString() << " sps=" << local << " nnz=" << nnz;
  }
  optimization->constraints.view->Done(); // mandatory after traversing the view

  return nnz;
}


bool SCPIP::get_bounds_info(int n, double* x_l, double* x_u,
                            int m, double* g_l, double* g_u)
{
  GetBounds(n, x_l, x_u, m, g_l, g_u);

  return true;
}

bool SCPIP::get_starting_point(int n, double* x)
{
  // we initialize x in bounds, in the upper right quadrant
  optimization->GetDesign()->WriteDesignToExtern(x);

  return true;
}

bool SCPIP::eval_f(int n, const double* x_org, double& obj_value)
{
  StdVector<double> x_srt;
  x_srt.Import(x_org, n);

  obj_value = EvalObjective(n, x_srt.GetPointer(), true); // always CFS scale!
  return true;
}

bool SCPIP::eval_grad_f(int n, const double* x_org, double* grad_f)
{
  assert(opt_timer->IsRunning());
  StdVector<double> x_srt;
  x_srt.Import(x_org, n);
  
  // restart_requested handled in intermediate_callback
  assert(grad_f == df.GetPointer());
  bool result = EvalGradObjective(n, x_srt.GetPointer(), true, df);

  assert(opt_timer->IsRunning());

  // do we have to write the initial iteration in the non-scale case?
  // SCPIP first does eval_f and then eval_grad_f
  if(optimization->GetCurrentIteration() == 0)
    CommitIteration();
  return result;
}

bool SCPIP::eval_g(int n, const double* x_org, int m, double* g)
{
  StdVector<double> x_srt;
  x_srt.Import(x_org, n);

  EvalConstraints(n, x_srt.GetPointer(), m, true, g, false); // we normalize in SCPIPBase

  return true;
}

bool SCPIP::eval_jac_g(int n, const double* x_org, int m, int nele_jac, double* values)
{
  StdVector<double> x_srt;
  x_srt.Import(x_org, n);

  // the gradients are dense in SCPIPBase
  assert(values != NULL);
  assert(jac_g.GetPointer() == values);
  
  EvalGradConstraints(n, x_srt.GetPointer(), m, nele_jac, true, false, jac_g);  // we normalize in SCPIPBase

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

  assert(opt_timer->IsRunning());
  // save this iteration, otherwise it might be lost
  CommitIteration();
  
  if(restart_requested) return;
  
  std::cout << "SCPIP finished: f=" << obj_value;
  for(int i = 0; i < m; i++)
    if(!optimization->constraints.view->Get(i)->IsVirtual()) // don't display blown up slopes
      std::cout << " + " << lambda[i] << "*" << g[i];
  std::cout << std::endl;  

  optimization->constraints.view->Done(); // swith slope constraints back to global
}
 
bool SCPIP::intermediate_callback(int iter, bool next_iter)
{
  if(iter == 0 || !next_iter) return true;

  CommitIteration();

  // break the optimization - e.g. if our relative change is smaller than given in xml
  // or we have to do a restart
  LOG_TRACE2(scpip) << "intermediate_callback iter=" << iter << " next_iter=" << next_iter
                    << " restart_requested=" << restart_requested;
  LOG_DBG(scpip) << "ic: mactiv=" << mactiv;
  LOG_DBG2(scpip) << "ic: active=" << active.ToString();
  
  return (restart_requested || optimization->DoStopOptimization()) ? false : true;
}     

