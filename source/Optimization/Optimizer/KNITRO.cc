#include <fstream>
#include "Optimization/Optimizer/KNITRO.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Condition.hh"
#include "General/exception.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/tools.hh"
#include "Utils/Timer.hh"


// declare class specific logging stream
DECLARE_LOG(knitro)
DEFINE_LOG(knitro, "knitro")


using std::string;

namespace CoupledField
{

KNITRO::KNITRO(Optimization* opt, PtrParamNode pn)
  : BaseOptimizer(opt, pn, Optimization::KNITRO_SOLVER)
{
  kc = KTR_new();
  if(kc == NULL)
    throw Exception("Failed to find a Ziena license for the KNITRO optimizer");
  
  // number of design variables
  n = optimization->GetDesign()->GetNumberOfVariables();
  // number of constraints
  m = optimization->constraints.view->GetNumberOfActiveConstraints();
  // number of non-sparse constraint gradients
  nnzJ = 0;
  for(int i = 0; i < m; i++)
    nnzJ += optimization->constraints.view->Get(i)->GetSparsityPattern().GetSize();
  optimization->constraints.view->Done(); // mandatory after traversing the view

  // allocate stuff
  x.Resize(n);
  lambda.Resize(m+n);
  c.Resize(std::max(m,1));
  objGrad.Resize(n);
  jac.Resize(std::max(nnzJ,1));


  BaseOptimizer::PostInitScale(1.0);
  // call the initialize function of KNITRO
  Init(pn);
}
  
KNITRO::~KNITRO()
{
  KTR_free(&kc); // frees the license
}

void KNITRO::Init(PtrParamNode pn)
{
  // initialized by the constructor.

  // set default parameters
  SetParameter("algorithm", KTR_ALG_BAR_CG);
  SetParameter("hessopt", KTR_HESSOPT_LBFGS); // actually we want a Hesse free algorithm but just in case
  SetParameter("outlev", 3); //  print basic information at each iteration
  SetParameter("maxit", optimization->GetMaxIterations());

  // get the parameters from the xml file and eventually overwrite the defaults
  if(pn->Has("knitro"))
  {
    ParamNodeList options = pn->Get("knitro")->GetList("option");

    for(unsigned int i = 0; i < options.GetSize(); i++)
    {
      if(options[i]->Get("type")->As<string>() == "integer")
        SetParameter(options[i]->Get("name")->As<string>(), options[i]->Get("value")->As<int>());
      else
        SetParameter(options[i]->Get("name")->As<string>(), options[i]->Get("value")->As<double>());
    }
  }

  // KNITRO copied the data within setup
  StdVector<double> xLoBnds(n); // set by GetBounds()
  StdVector<double> xUpBnds(n); // set by GetBounds()
  StdVector<double> xInitial(n);
  StdVector<int>    cType(std::max(m,1));
  StdVector<double> cLoBnds(std::max(m,1)); // set by GetBounds()
  StdVector<double> cUpBnds(std::max(m,1)); // set by GetBounds()
  StdVector<int>    jacIndexVars(std::max(nnzJ, 1));
  StdVector<int>    jacIndexCons(std::max(nnzJ, 1));

  GetBounds(n, xLoBnds.GetPointer(), xUpBnds.GetPointer(), m, cLoBnds.GetPointer(), cUpBnds.GetPointer());
  
  optimization->GetDesign()->WriteDesignToExtern(xInitial.GetPointer());

  // set constraint function type. We know no quadratic functions in CFS
  for(int i = 0; i < m; i++)
    cType[i] = optimization->constraints.view->Get(i)->IsLinear() ? KTR_CONTYPE_LINEAR : KTR_CONTYPE_GENERAL;
  optimization->constraints.view->Done(); // mandatory after traversing the view

  // Set the sparse pattern structure of the constraints
  int idx = 0;
  for(int c = 0; c < m; c++)
  {
    StdVector<unsigned int>& pattern = optimization->constraints.view->Get(c)->GetSparsityPattern();
    for(unsigned int p = 0; p < pattern.GetSize(); p++)
    {
      jacIndexVars[idx] = pattern[p];
      jacIndexCons[idx] = c;
      idx++;
    }
  }
  optimization->constraints.view->Done();
  assert(idx == nnzJ);

  // now set the parameters
  int objGoal = KTR_OBJGOAL_MINIMIZE;  // -1 multiplied in maximization case in BaseOptimizer
  int objType = KTR_OBJTYPE_GENERAL;   // we make in CFS no difference to linear and quadratic

  LOG_DBG3(knitro) << "Init: xLoBnds=" << xLoBnds.ToString();
  LOG_DBG3(knitro) << "Init: xUpBnds=" << xUpBnds.ToString();
  LOG_DBG3(knitro) << "Init: cType=" << cType.ToString();
  LOG_DBG3(knitro) << "Init: cLoBnds=" << cLoBnds.ToString();
  LOG_DBG3(knitro) << "Init: cUpBnds=" << cUpBnds.ToString();
  LOG_DBG3(knitro) << "Init: jacIndexVars=" << jacIndexVars.ToString();
  LOG_DBG3(knitro) << "Init: jacIndexCons=" << jacIndexCons.ToString();
  LOG_DBG3(knitro) << "Init: xInitial=" << xInitial.ToString();

  // initializing copies the values. No Hessian and no initial lambda
  int status = KTR_init_problem (kc, n, objGoal, objType, xLoBnds.GetPointer(), xUpBnds.GetPointer(),
                                 m, cType.GetPointer(), cLoBnds.GetPointer(), cUpBnds.GetPointer(),
                                 nnzJ, jacIndexVars.GetPointer(), jacIndexCons.GetPointer(),
                                 0, NULL, NULL, xInitial.GetPointer(), NULL);

  if(status != 0) EXCEPTION("initializing KNITRO failed with status " << status);

}

void KNITRO::SetParameter(const string& name, int value)
{
  int status = KTR_set_int_param_by_name (kc, name.c_str(), value);
  if(status != 0)
    EXCEPTION("Setting KNITRO integer option " << name << " to " << value << " failed");
}

void KNITRO::SetParameter(const string& name, double value)
{
  int status = KTR_set_double_param_by_name(kc, name.c_str(), value);
  if(status != 0)
    EXCEPTION("Setting KNITRO float option " << name << " to " << value << " failed");
}


void KNITRO::SolveProblem()
{
  // use reverse communication, is simpler to debug than callbacks

  int status = -1;

  int max_iter = optimization->GetMaxIterations();
  bool leave = false;

  int iter = 0; // to find the next iteration to do a commit

  while(!optimization->DoStopOptimization() && optimization->GetCurrentIteration() <= max_iter && !leave)
  {
    int evalStatus = 0; // one if we cannot evaluate functions/derivatives

    status = KTR_solve(kc, x.GetPointer(), lambda.GetPointer(), evalStatus,
                       &obj, c.GetPointer(), objGrad.GetPointer(), jac.GetPointer(),
                       NULL, NULL, NULL); // no Hessian

    LOG_DBG2(knitro) << "SP: s=" << status << " it=" << KTR_get_number_iters(kc) << " mji=" << KTR_get_number_major_iters(kc) << " mni=" << KTR_get_number_minor_iters(kc);

    switch(status)
    {
    case KTR_RC_EVALFC: // KNITRO WANTS obj AND c EVALUATED AT THE POINT x
    {
      obj = EvalObjective(n, x.GetPointer(), true);
      EvalConstraints(n, x.GetPointer(), m, true, c.GetPointer());
      int curr_iter = KTR_get_number_iters(kc);
      if(curr_iter != iter)
      {
        optimization->CommitIteration();
        iter = curr_iter;
      }
      break;
    }

    case KTR_RC_EVALGA: // KNITRO WANTS objGrad AND jac EVALUATED AT x
    {
      bool ego = EvalGradObjective(n, x.GetPointer(), true, objGrad);
      if(!ego) EXCEPTION("don't do autoscale with KNITRO ");
      EvalGradConstraints(n, x.GetPointer(), m, nnzJ, true, jac);
      break;
    }

    case KTR_RC_EVALHV:
    case KTR_RC_EVALH:
      EXCEPTION("Hessian for KNITRO not implemented");

    case KTR_RC_NEWPOINT:
      EXCEPTION("Don't know how to handle New Point in KNITRO");

    default:
      leave = true; // finished or error
    }
  }
  PtrParamNode in = optimization->optInfoNode->Get(ParamNode::SUMMARY)->Get("break");
  
  in->Get("feasibility_err_abs")->SetValue(KTR_get_abs_feas_error(kc));
  in->Get("feasibility_err_rel")->SetValue(KTR_get_rel_feas_error(kc));
  in->Get("optimality_err_abs")->SetValue(KTR_get_abs_opt_error(kc));
  in->Get("optimality_err_rel")->SetValue(KTR_get_rel_opt_error(kc));
  in->Get("majors")->SetValue(KTR_get_number_major_iters(kc));
  in->Get("minors")->SetValue(KTR_get_number_minor_iters(kc));
  in->Get("status")->SetValue(status);


  optimization->CommitIteration();

  if(status == 0)
  {
    in->Get("converged")->SetValue("yes");
    in->Get("reason/msg")->SetValue("KNITRO converged");
    return;
  }
  if(!leave)
  {
    assert(optimization->DoStopOptimization() || optimization->GetCurrentIteration() >= max_iter);
    return; // message written by Optimization
  }

  // error
  in->Get("converged")->SetValue("no");
  in->Get("reason/msg")->SetValue("KNITRO returned " + boost::lexical_cast<string>(status));

  EXCEPTION("KNITRO returned " << status);
}

void KNITRO::LogFileLine(std::ofstream* out, PtrParamNode iteration)
{
  *out << "\t" << KTR_get_abs_feas_error(kc) << "\t" << KTR_get_abs_opt_error(kc);
  
  iteration->Get("feasibility")->SetValue(KTR_get_abs_feas_error(kc));
  iteration->Get("optimality")->SetValue(KTR_get_abs_opt_error(kc));
}



} // namespace
