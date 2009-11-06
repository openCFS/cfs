#include "Optimization/SnOptInterface.hh"
#include "Optimization/SnOpt.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/Condition.hh"
#include "General/exception.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{

// declare class specific logging stream
DECLARE_LOG(snopt)
DEFINE_LOG(snopt, "snopt")

using std::endl;
using std::cout;

// use the pointer defined in SnOptBase.cc!
extern SnOptBase *static_snopt_base;

/** this is a global function that can be called in the snopt-interface
 *  it simply forwards the call to the SnOpt->Callback function which has access to the class members */
int SnOpt_C_Callback(int* Status, int *n,
    double* x, int* needF, int* nF, double* F,
    int* needG, int* nG, double* G,
    char* cu, int* lencu,
    int* iu, int* leniu,
    double* ru, int* lenru)
{
  return static_snopt_base->Callback(Status, *n, x, needF, nF, F, needG, nG, G, cu, lencu, iu, leniu, ru, lenru);
}

SnOpt::SnOpt(Optimization* opt, ParamNode* pn) :
  BaseOptimizer(opt, pn),
  SnOptBase(),
  optimizer_pn_(pn->Get(Optimization::optimizer.ToString(Optimization::SNOPT_SOLVER), false))
{
  LOG_TRACE(snopt) << "Initialize SnOpt";
  
  double manual_scaling(1.0);
  PostInit(manual_scaling);
  
  Init();
}

SnOpt::~SnOpt()
{
  LOG_TRACE(snopt) << "number of function evaluations = " << f_evals;
  LOG_TRACE(snopt) << "number of gradient evaluations = " << g_evals;
  
  snoptclose_(&iPrint);
}

void SnOpt::Init()
{
  LOG_TRACE(snopt) << "Init";
  
  // ================================================================== 
  // First,  sninit_ MUST be called to initialize optional parameters   
  // to their default values.                                           
  //  ==================================================================
  cw.resize(8*mincw);
  rw.resize(minrw);
  iw.resize(miniw);
  sninit_(&iPrint, &iSumm, &cw[0], &mincw, &iw[0], &miniw, &rw[0], &minrw, 8*mincw);
  
  // initialize problem
  get_nlp_info(n, nF, nG);
  
  // set the optimization options
  // call after get_nlp_info!! 
  SetSnOptOptions();
  
  x.resize(n, 1.0); // start value
  xlow.resize(n, -GetInfBound());
  xupp.resize(n,  GetInfBound());
  xmul.resize(n,  0.0);
  xstate.resize(n, 0);

  F.resize(nF, 0.0);
  Flow.resize(nF, -GetInfBound());
  Fupp.resize(nF,  GetInfBound());
  Fmul.resize(nF, 0.0);
  Fstate.resize(nF, 0);
  
  G.resize(nG, 0.0);
  iGfun.resize(lenG, 0);
  jGvar.resize(lenG, 0);

  for(int r = 0; r < nF; ++r)
  {
    for(int c = 0; c < n; ++c)
    {
      // the jacobian array/matrix is a long vector with n*m elements
      // see the IPOPT example in the docu to understand!
      int index = r * n + c;
      iGfun[index] = r+1;
      jGvar[index] = c+1;
    }
  }
  
  lenA = lenG;
  // FIXME this is not correct!
  iAfun.resize(lenA, 1);
  jAvar.resize(lenA, 1);
  A.resize(lenA, 0.0);
  /////////////////////////////////
  
  get_bounds_info(n, &xlow[0], &xupp[0], nF - 1, &Flow[1], &Fupp[1]);
  
  optimization->GetDesign()->WriteBoundsToExtern(&xlow[0], &xupp[0]);
}

void SnOpt::SolveProblem()
{
  LOG_TRACE(snopt) << "SolveProblem";
  
  // set start parameters here to allow restarted SolveProblem()
  get_starting_point(n, &x[0]);
  
  // this must always be called AFTER sninit!!!
  AdjustWorkArrayMemory();

  // before doing anything, we save the current state
  if(optimization->GetCurrentIteration() == 0)
  {
    optimization->SolveStateProblem();
    optimization->CommitIteration();
  }

  // this is needed for the call to snopt, but has no effect in our case
  int npname(5);
  char xnames[8] = "unused ";
  char Fnames[8] = "unused ";
  
  // here we could put a more suitable problem name...
  char Prob[9] = "FIXME   ";
  
  // ------------------------------------------------------------------
  // Go for it, using a Cold start.                                    
  // ------------------------------------------------------------------
  snopta_(
      &Start, &nF, &n, &nxname, &nFname,
      &ObjAdd, &ObjRow, Prob, SnOpt_C_Callback, snlog_,
      &iAfun[0], &jAvar[0], &lenA, &nA, &A[0],
      &iGfun[0], &jGvar[0], &lenG, &nG,
      &xlow[0], &xupp[0], xnames,  &Flow[0], &Fupp[0], Fnames,
      &x[0], &xstate[0], &xmul[0], &F[0], &Fstate[0], &Fmul[0],
      &INFO, &mincw, &miniw, &minrw, &nS, &nInf, &sInf,
      &cw[0], &lencw, &iw[0], &leniw, &rw[0], &lenrw,
      &cw[0], &lencw, &iw[0], &leniw, &rw[0], &lenrw,
      npname, 8*nxname, 8*nFname, 8*lencw, 8*lencw
      );
  
  LOG_TRACE(snopt) << "INFO = " << INFO << ", nS = " << nS;
  
  optimization->CommitIteration();
}

int SnOpt::Callback(int* Status, const int n,
    double* x, int* needF, int* nF, double* F,
    int* needG, int* nG, double* G,
    char* cu, int* lencu, int* iu, int* leniu, double* ru, int* lenru)
{
  LOG_TRACE(snopt) << "Callback";
  
  // if we want to stop, we have to set Status to a value < -1!
  if(optimization->DoStopOptimization() || stop)
  {
    stop = true;
    *Status = -10;
    cout << "user break!" << endl;
    return -10;
  }
  
  
  // it would be really nice at this point to know in which iteration
  // snopt currently is...
  /*
  int it(0);
  sngeti_(const_cast<char*>("Iter"), &it, &INFO,
        &cw[0], &lencw, &iw[0], &leniw, &rw[0], &lenrw, 4, 8*lencw );
  cout << "itn = " << it << endl;
  */
  
  // do we need a function evaluation?
  if(*needF > 0)
  {
    ++f_evals;
    
    // The objective row
    // we assume that ObjRow is equal to 1!!
    eval_f(n, x, F[ObjRow - 1]);
    
    // evaluate the constraints
    // we need to calculate (nF-1) constraints
    // we put them in F+1, since the first entry is the objective!
    // we assume that ObjRow is equal to 1!!
    eval_g(n, x, *nF-1, F + 1);
  }
  
  if(*needG > 0)
  {
    ++g_evals;
    
    // evaluate the gradient of the cost function
    // this again assumes that ObjRow is equal to 1!!
    eval_grad_f(n, x, G);
    
    // evaluate the constraint gradient
    // here we do the same pointer logic as for the function evaluations
    eval_jac_g(n, x, *nF - 1, *nG, G + n);
    
    // for every call to the gradient we also commit the iteration, because we have
    // yet to find a way to get the iteration counter from snopt...
    optimization->CommitIteration();
  }
  
  // everything went okay
  return 0;
}

bool SnOpt::get_nlp_info(int& n, int& nF, int &nG)
{
  LOG_TRACE(snopt) << "get_nlp_info";
  
  // number of design variables
  n = optimization->GetDesign()->GetNumberOfVariables();

  // number of problem functions (objective and constraints, linear and nonlinear)
  nF = optimization->objectives.data.GetSize() + optimization->constraints.GetSize();
  
  // number of nonzeros of jacobian
  // we use a dense jacobian, so every entry is != 0
  nG = n * nF;
  lenG = nG;

  LOG_TRACE(snopt) << "Problem dimensions: n = " << n << ", nF = " << nF << ", nG = " << nG;
  
  return true;
}

bool SnOpt::get_bounds_info(int n, double* xlow, double* xupp, int m, double* g_l, double* g_u)
{
  LOG_TRACE(snopt) << "get_bounds_info";
  
  GetBounds(n, xlow, xupp, m, g_l, g_u);
  return true;
}

bool SnOpt::get_starting_point(int n, double* x)
{
  LOG_TRACE(snopt) << "get_starting_point";
  
  // we initialize x in bounds, in the upper right quadrant
  optimization->GetDesign()->WriteDesignToExtern(x);

  return true;
}

bool SnOpt::eval_f(int n, const double* x, double &obj_value)
{
  LOG_TRACE(snopt) << "eval_f";
  
  obj_value = EvalObjective(n, x);
  
  return true;
}

bool SnOpt::eval_grad_f(int n, const double* x, double* grad_f)
{
  LOG_TRACE(snopt) << "eval_grad_f";
  
  // restart_requested handled in intermediate_callback
  bool result = EvalGradObjective(n, x, grad_f);

  return result;
}

bool SnOpt::eval_g(int n, const double* x, int m, double* g)
{
  LOG_TRACE(snopt) << "eval_g";
  
  EvalConstraints(n, x, m, g);
  
  return true;
}

bool SnOpt::eval_jac_g(int n, const double* x, int m, int nele_jac, double* values)
{
  LOG_TRACE(snopt) << "eval_jac_g";
  
  // the gradients are dense in snopt
  assert(values != NULL);
  
  EvalGradConstraints(n, x, m, nele_jac, values);

  return true;
}

void SnOpt::AdjustWorkArrayMemory()
{
  LOG_TRACE(snopt) << "AdjustWorkArrayMemory";
  LOG_TRACE(snopt) << "old values: lencw = " << lencw << ", leniw = " << leniw << ", lenrw = " << lenrw;
  
  // remember old values
  int tmpcw(lencw);
  int tmpiw(leniw);
  int tmprw(lenrw);
  
  // try to determine minimal amount of memory needed for this problem
  snmema_(&nF, &n, &nxname, &nFname, &nA, &nG, &mincw, &miniw, &minrw,
      &cw[0], &lencw, &iw[0], &leniw, &rw[0], &lenrw, 8*lencw);
  
  LOG_TRACE(snopt) << "After snmema_:";
  LOG_TRACE(snopt) << "min values: mincw = " << mincw << ", miniw = " << miniw << ", minrw = " << minrw;

  // enlarge the minimal values determined by snmema_ by this factor
  // this is because according to the snopt manual, the minimal values
  // might not be enough
  const double factor(1.2);
  
  // update lengths according to values obtained from snmema_
  if(lencw < mincw)
  {
    lencw = static_cast<int>(factor * static_cast<double>(mincw));
    cw.resize(8*lencw);
    snseti_("Total character workspace", &lencw, &iPrint, &iSumm, &INFO,
            &cw[0], &tmpcw, &iw[0], &tmpiw, &rw[0], &tmprw, 25, 8*lencw);
    LOG_TRACE(snopt) << "new value: lencw = " << lencw;
  }
  if(leniw < miniw)
  {
    leniw = static_cast<int>(factor * static_cast<double>(miniw));
    iw.resize(leniw);
    snseti_("Total integer workspace", &leniw, &iPrint, &iSumm, &INFO,
            &cw[0], &tmpcw, &iw[0], &tmpiw, &rw[0], &tmprw, 23, 8*lencw);
    LOG_TRACE(snopt) << "new value: leniw = " << leniw;
  }
  if(lenrw < minrw)
  {
    lenrw = static_cast<int>(factor * static_cast<double>(minrw));
    rw.resize(lenrw);
    snseti_("Total real workspace", &lenrw, &iPrint, &iSumm, &INFO,
            &cw[0], &tmpcw, &iw[0], &tmpiw, &rw[0], &tmprw, 20, 8*lencw);
    LOG_TRACE(snopt) << "new value: lenrw = " << lenrw;
  }
}

void SnOpt::SetSnOptOptions()
{
  LOG_TRACE(snopt) << "SetSnOptOptions";
  
  // check for maximization instead of minimization
  if(optimization->objectives.DoMaximize())
    SetIntegerValue("maximize", 1);
  
  // additionally, set the number of major iterations
  SetIntegerValue("major_iterations_limit", optimization->GetMaxIterations());

  // check for optional parameters from xml
  if(optimizer_pn_ != NULL)
  {
    StdVector<ParamNode*> list;
    
    list = optimizer_pn_->GetList("option", "type", "integer");
    for(unsigned int i = 0; i < list.GetSize(); i++)
      SetIntegerValue(list[i]->Get("name")->AsString(), list[i]->Get("value")->AsInt());
    
    list = optimizer_pn_->GetList("option", "type", "real");
    for(unsigned int i = 0; i < list.GetSize(); i++)
      SetNumericValue(list[i]->Get("name")->AsString(), list[i]->Get("value")->AsDouble());
  }
}

} // namespace
