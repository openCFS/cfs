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

/** global pointer to the snopt class to be used by the callback function */
SnOpt* static_snopt = NULL;

/** this is a global function that can be called in the snopt-interface
 *  it simply forwards the call to the SnOpt->Callback function which has access to the class members */
int SnOpt_C_Callback(int* Status, int *n,
    double* x, int* needF, int* nF, double* F,
    int* needG, int* nG, double* G,
    char* cu, int* lencu, int* iu, int* leniu, double* ru, int* lenru)
{
  return static_snopt->Callback(Status, *n, x, needF, nF, F, needG, nG, G, cu, lencu, iu, leniu, ru, lenru);
}

SnOpt::SnOpt(Optimization* opt, PtrParamNode pn) :
  BaseOptimizer(opt, pn),
  f_evals(0),       // number of function evaluations
  g_evals(0),       // number of gradient evaluations
  nxname(1),
  nFname(1),
  iSumm(6),         // variable for snoptbase output
  iPrint(1),       // variable for file output
  iSpecs(4),        // variable for file output
  DerOpt(1),
  minrw(500),
  miniw(500),
  mincw(500),
  lenrw(500),
  leniw(500),
  lencw(500),
  n(0),              // number of variables
  nF(0),            // number of cost functions (=1) and constraints
  nG(0),            // number of nonzeros in nonlinear gradient
  lenG(0),
  nA(0),            // number of nonzeros in linear gradient
  lenA(0),
  nS(0),          // number of superbasic variables
  nInf(0),
  sInf(0.0),
  Start(0),          // start with a cold start
  stop(false),     // snopt is a little bit hard to kill, so we need a stopper
  INFO(0),
  ObjRow(1),      // objective row, must be 1!!
  ObjAdd(0.0),
  optimizer_pn_(pn->Get(Optimization::optimizer.ToString(Optimization::SNOPT_SOLVER), ParamNode::PASS)),
  lin_constraints(0),
  nonlin_constraints(0)
{
  LOG_TRACE(snopt) << "Initialize SnOpt";
  
  static_snopt = this;
  
  PostInit(1.0); // BaseOptimizer
  
  Init();
}

SnOpt::~SnOpt()
{
  LOG_TRACE(snopt) << "number of function evaluations = " << f_evals;
  LOG_TRACE(snopt) << "number of gradient evaluations = " << g_evals;
  
  snclose_(&iPrint);
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
  get_nlp_info();
  
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
  
  // init Jacobian(s)
  initJacobians();

  LOG_TRACE(snopt) << "get_bounds_info";
  GetBounds(n, &xlow[0], &xupp[0], nF - 1, &Flow[1], &Fupp[1]);
  
  //for(unsigned int i = 0; i < Flow.size(); ++i)
    //std::cout << "Flow[" << i << "] = " << Flow[i] << ", Fupp[" << i << "] = " << Fupp[i] << std::endl;
  
  optimization->GetDesign()->WriteBoundsToExtern(&xlow[0], &xupp[0]);
  
  gradhelper.Resize(std::max(n, nG - n), 0.0);
}

void SnOpt::SolveProblem()
{
  LOG_TRACE(snopt) << "SolveProblem";
  
  // set start parameters here to allow restarted SolveProblem()
  // we initialize x in bounds, in the upper right quadrant
  LOG_TRACE(snopt) << "get_starting_point";
  optimization->GetDesign()->WriteDesignToExtern(&x[0]);
  
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

  // if we have linear constraints we have to setup the gradient
  if(lin_constraints > 0)
  {
    setupLinearConstraints();
  }  
  else
  {
    lenA = 0;
    nA = 0;
    assert(A.size() == 0);
    assert(iAfun.size() == 0);
    assert(jAvar.size() == 0);
  }
  
  snopta_(
      &Start, &nF, &n, &nxname, &nFname,
      &ObjAdd, &ObjRow, Prob, SnOpt_C_Callback,
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
  
  // do we need a function evaluation?
  if(*needF > 0)
  {
    ++f_evals;
    
    // The objective row
    // we assume that ObjRow is equal to 1!!
    assert(ObjRow == 1);
    eval_f(n, x, F[0]);
    
    // evaluate the constraints
    // we need to calculate (nF-1) constraints
    // we put them in F+1, since the first entry is the objective!
    // we assume that ObjRow is equal to 1!!
    eval_g(n, x, *nF - 1, &F[1]);
  }
  
  if(*needG > 0)
  {
    ++g_evals;
    
    // evaluate the gradient of the cost function
    // this again assumes that ObjRow is equal to 1!!
    eval_grad_f(n, x, G);
    
    // evaluate the constraint gradient
    // here we do the same pointer logic as for the function evaluations
    // input:
    // 1. n = number of design variables
    // 2. x = design variable values
    // 3. number of nonlinear constraints (for dense *nF-1)
    // 4. number of nonlinear constraint gradient values (for dense *nG-n)
    // 5. values of nonlinear gradient (always G+n)
    
    //assert(*nG-n == nonlin_constraints*n); // this is only valid if nonlin constr are DENSE!!!!
    
    eval_jac_g(n, x, nonlin_constraints, *nG - n, G + n);
    
    // for every call to the gradient we also commit the iteration, because we have
    // yet to find a way to get the iteration counter from snopt...
    optimization->CommitIteration();
  }
  
  // everything went okay
  return 0;
}

bool SnOpt::get_nlp_info()
{
  LOG_TRACE(snopt) << "get_nlp_info";
  
  // number of design variables
  n = optimization->GetDesign()->GetNumberOfVariables();

  // number of problem functions (objective and constraints, linear and nonlinear)
  nF = 1 + optimization->constraints.view->GetNumberOfActiveConstraints();
  
  // walk over all constraints and check for linearity
  nonlin_constraints = 0;
  lin_constraints = 0;

  for(int c = 0, nc = optimization->constraints.view->GetNumberOfActiveConstraints(); c < nc; ++c)
  {
    Condition* g = optimization->constraints.view->Get(c);
    assert(g->GetSparsityPattern().GetSize() > 0);

    if(g->IsLinear())
    {
      ++lin_constraints;
      nA += g->GetSparsityPattern().GetSize(); // zero slopes
    }
    else
    {
      ++nonlin_constraints;
      nG += g->GetSparsityPattern().GetSize();
    }
  }
  
  // number for nonlin gradient must be higher for objective function 
  nG += n;
  
  optimization->constraints.view->Done();

  LOG_TRACE(snopt) << "Problem dimensions: n = " << n << ", nF = " << nF << ", nG = " << nG << ", nA = " << nA;
  LOG_TRACE(snopt) << "Problem dimensions: lin_constraints = " << lin_constraints << ", nonlin_constraints = " << nonlin_constraints;
  
  assert(nA >= lin_constraints);
  
  return true;
}

bool SnOpt::eval_f(int n, const double* x, double &obj_value)
{
  LOG_TRACE(snopt) << "eval_f";
  obj_value = EvalObjective(n, x, true); // as with SCPIP we do always autoscale!
  return true;
}

bool SnOpt::eval_grad_f(int n, const double* x, double* grad_f)
{
  LOG_TRACE(snopt) << "eval_grad_f";
  assert(static_cast<int>(gradhelper.GetSize()) >= n);

  // restart_requested handled in intermediate_callback
  bool result = EvalGradObjective(n, x, true, gradhelper);

  for(int i = 0; i < n; i++)
  {
    //LOG_DBG3(snopt) <<   "old grad_f[" << i << "] = " << grad_f[i] 
    //                << ", new grad_f[" << i << "] = " << gradhelper[i];
    grad_f[i] = gradhelper[i];
  }

  return result;
}

bool SnOpt::eval_g(int n, const double* x, int m, double* g)
{
  LOG_TRACE(snopt) << "eval_g";
  
  //for(int i = 0; i < m; i++) LOG_DBG3(snopt) <<   "old G[" << i << "] = " << g[i]; 

  EvalConstraints(n, x, m, true, g);
  
  //for(int i = 0; i < m; i++) LOG_DBG3(snopt) <<   "new G[" << i << "] = " << g[i]; 

  return true;
}

bool SnOpt::eval_jac_g(int n, const double* x, int m, int nele_jac, double* pG)
{
  LOG_TRACE(snopt) << "eval_jac_g: n=" << n << ", m=" << m << ", nele_jac=" << nele_jac;
  assert(static_cast<int>(gradhelper.GetSize()) >= nele_jac);

  // the gradients are dense in snopt
  assert(pG != NULL);
  assert(m == nonlin_constraints);

  int nc(optimization->constraints.view->GetNumberOfActiveConstraints());
 
  EvalGradConstraints(n, x, nc, nele_jac, true, gradhelper, BaseOptimizer::NONLINEAR);
  
  for(int i = 0; i < nele_jac; i++)
  {
    // LOG_DBG3(snopt) <<   "old grad_G[" << i << "] = " << pG[i] << ", new grad_G[" << i << "] = " << gradhelper[i];
    pG[i] = gradhelper[i];
  }

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
  snmema_(&INFO, &nF, &n, &nxname, &nFname, &nA, &nG, &mincw, &miniw, &minrw,
      &cw[0], &lencw, &iw[0], &leniw, &rw[0], &lenrw, 8*lencw);
  
  LOG_TRACE(snopt) << "After snmema_:";
  LOG_TRACE(snopt) << "min values: mincw = " << mincw << ", miniw = " << miniw 
                   << ", minrw = " << minrw << ", INFO = " << INFO;

  // enlarge the minimal values determined by snmema_ by this factor
  // this is because according to the snopt manual, the minimal values
  // might not be enough
  const double factor(1.0);
  
  // update lengths according to values obtained from snmema_
  if(lencw < mincw)
  {
    lencw = static_cast<int>(factor * static_cast<double>(mincw));
    cw.resize(8*lencw);
    snseti_("Total character workspace", &lencw, &iPrint, &iSumm, &INFO,
            &cw[0], &tmpcw, &iw[0], &tmpiw, &rw[0], &tmprw, 25, 8*lencw);
    LOG_TRACE(snopt) << "new value: lencw = " << lencw << ", INFO = " << INFO;
  }
  
  if(leniw < miniw)
  {
    leniw = static_cast<int>(factor * static_cast<double>(miniw));
    iw.resize(leniw);
    snseti_("Total integer workspace", &leniw, &iPrint, &iSumm, &INFO,
            &cw[0], &tmpcw, &iw[0], &tmpiw, &rw[0], &tmprw, 23, 8*lencw);
    LOG_TRACE(snopt) << "new value: leniw = " << leniw << ", INFO = " << INFO;
  }
  
  if(lenrw < minrw)
  {
    lenrw = static_cast<int>(factor * static_cast<double>(minrw));
    rw.resize(lenrw);
    snseti_("Total real workspace", &lenrw, &iPrint, &iSumm, &INFO,
            &cw[0], &tmpcw, &iw[0], &tmpiw, &rw[0], &tmprw, 20, 8*lencw);
    LOG_TRACE(snopt) << "new value: lenrw = " << lenrw << ", INFO = " << INFO;
  }
}

void SnOpt::SetSnOptOptions()
{
  LOG_TRACE(snopt) << "SetSnOptOptions";
  
  const int minItLimit(900000);
  
  // must make sure these options are set
  bool setMinorItLimit(false);
  bool setItLimit(false);
  
  // additionally, set the number of major iterations
  SetIntegerValue("major_iterations_limit", optimization->GetMaxIterations());
  
  // snopt is a little bit of a resource hog, but we do not care
  // set superbasics_limit to maximum (which might be a gross overestimation,
  // but being told that the limit is too small after waiting 38 hours is just
  // plain stupid)
  SetIntegerValue("superbasics_limit", n);

  // check for optional parameters from xml
  if(optimizer_pn_ != NULL)
  {
    ParamNodeList list;
    list = optimizer_pn_->GetListByVal("option", "type", "integer");
    
    for(unsigned int i = 0; i < list.GetSize(); i++)
    {
      string name(list[i]->Get("name")->As<std::string>());
      int value(list[i]->Get("value")->As<int>());
      
      if(name == "superbasics_limit")
        EXCEPTION("superbasics_limit is set automatically, do not specify in xml-file!");
      
      if(name == "minor_iterations_limit")
      {
        value = value < minItLimit ? minItLimit : value;
        setMinorItLimit = true;
      }
      
      if(name == "iterations_limit")
      {
        value = value < minItLimit ? minItLimit : value;
        setItLimit = true;
      }
      
      SetIntegerValue(name, value);
    }
    
    list = optimizer_pn_->GetListByVal("option", "type", "real");
    
    for(unsigned int i = 0; i < list.GetSize(); i++)
      SetNumericValue(list[i]->Get("name")->As<std::string>(), list[i]->Get("value")->As<double>());
  }
  
  if(!setMinorItLimit) SetIntegerValue("minor_iterations_limit", minItLimit);
  if(!setItLimit) SetIntegerValue("iterations_limit", minItLimit);
}

void SnOpt::initJacobians()
{
  LOG_TRACE(snopt) << "initJacobians()";
  
  // create Jacobian
  G.resize(1, 0.0); // FIXME: G is not even used by snopt...
  iGfun.reserve(nG);
  jGvar.reserve(nG);

  A.resize(nA, 0.0);
  iAfun.reserve(nA);
  jAvar.reserve(nA);

  // it is very important that we use the same row counter for all the
  // constraints, linear and nonlinear!!!
  int cstr(1);    // fortran!

  int index = 0; // counter
  int indexlin = 0; // counter

  // for objective function(s) which in cfs is only one
  for(int e = 1; e <= n; ++e) // fortran
  {
    iGfun.push_back(cstr);
    jGvar.push_back(e);
    
    LOG_DBG3(snopt) << "cost function"
                    << "; iGfun[" << index << "] = " << iGfun[index]
                    << ", jGvar[" << index << "] = " << jGvar[index];
    
    ++index;
  }
  ++cstr;

  LOG_DBG3(snopt) << "Flow.size() = " << Flow.size() << ", Fupp.size() = " << Fupp.size();
  
  for(int c = 0, nc = optimization->constraints.view->GetNumberOfActiveConstraints(); c < nc; c++)
  {
    Condition* g = optimization->constraints.view->Get(c);
    const StdVector<unsigned int>& pattern = g->GetSparsityPattern();
    assert(pattern.GetSize() > 0);

    if(g->IsLinear()) // linear
    {
      for(unsigned int e = 0; e < pattern.GetSize(); ++e)
      {
        iAfun.push_back(cstr);
        jAvar.push_back(pattern[e] + 1);

        LOG_DBG3(snopt) << "g = " << g->ToString() 
                          << "; iAfun[" << indexlin << "] = " << iAfun[indexlin]
                          << ", jAvar[" << indexlin << "] = " << jAvar[indexlin];

        ++indexlin;
      }
    }
    else // nonlinear 
    {
      // up to now only dense so we do not need anything
      for(unsigned int e = 0; e < pattern.GetSize(); ++e)
      {
        iGfun.push_back(cstr);
        jGvar.push_back(pattern[e] + 1);
        
        LOG_DBG3(snopt) << "g = " << g->ToString()
                        << "; iGfun[" << index << "] = " << iGfun[index]
                        << ", jGvar[" << index << "] = " << jGvar[index];
        ++index;
      }
    }
    // increment constraint counter
    ++cstr;
  }

  assert(index <= nG);
  assert(indexlin <= nA);

  lenG = index;
  lenA = indexlin;
  
  LOG_DBG(snopt) << "lenG = " << lenG << ", nG = " << nG << "; lenA = " << lenA << ", nA = " << nA;
  LOG_TRACE(snopt) << "end of initJacobians()";
}

void SnOpt::setupLinearConstraints()
{
  // setup linear constraints
  StdVector<double> lincon(nA);

  EvalGradConstraints(n, &x[0], optimization->constraints.view->GetNumberOfActiveConstraints(), nA,
      true, lincon, BaseOptimizer::LINEAR);
  
  for(int i = 0; i < nA; i++)
  {
    A[i] = lincon[i];
    LOG_DBG3(snopt) << "A[" << i << "] = " << A[i];
  }
  
  //assert(count == lin_constraints);
}

void SnOpt::SetIntegerValue(const std::string& key, int value)
{
  string option;
  
  if(key == "major_iterations_limit")
  {
    option = "Major iterations limit";
  }
  else if(key == "minor_iterations_limit")
  {
    option = "Minor iterations limit";
  }
  else if(key == "iterations_limit")
  {
    option = "Iterations limit";
  }
  else if(key == "superbasics_limit")
  {
    option = "Superbasics limit";
  }
  else if(key == "timing_level")
  {
    option = "Timing level";
  }
  else if(key == "verify_level")
  {
    option = "Verify level";
  }
  else if(key == "major_print_level")
  {
    option = "Major print level";
  }
  else if(key == "minor_print_level")
  {
    option = "Minor print level";
  }
  else if(key == "print_frequency")
  {
    option = "Print frequency";
  }
  else if(key == "factorization_frequency")
  {
    option = "Factorization frequency";
  }
  
  if(!option.empty())
  {
    LOG_TRACE(snopt) << "adjusted " << key;
    // set the option
    snseti_(option.c_str(), &value, &iPrint, &iSumm, &INFO,
        &cw[0], &lencw, &iw[0], &leniw, &rw[0], &lenrw, option.size(), 8*lencw);
    
    // check INFO, if after setting an option INFO > 0 we have had an error!
    assert(INFO == 0);
  }
  else // nothing found
    EXCEPTION("trying to set an snopt-option that is unknown!");
}

void SnOpt::SetNumericValue(const std::string& key, double value)
{  
  string option;
  
  if(key == "major_optimality_tolerance")
  {
    option = "Major optimality tolerance";
  }
  else if(key == "major_feasibility_tolerance")
  {
    option = "Major feasibility tolerance";
  }
  else if(key == "minor_feasibility_tolerance")
  {
    option = "Minor feasibility tolerance";
  }
  else if(key == "linesearch_tolerance")
  {
    option = "Linesearch tolerance";
  }
     
  if(!option.empty())
  {
    LOG_TRACE(snopt) << "adjusted " << key;
    // set the option
    snsetr_(option.c_str(), &value, &iPrint, &iSumm, &INFO,
        &cw[0], &lencw, &iw[0], &leniw, &rw[0], &lenrw, option.size(), 8*lencw );
    
    // check INFO, if after setting an option INFO > 0 we have had an error!
    assert(INFO == 0);
  }
  else // nothing found
    EXCEPTION("trying to set an snopt-option that is unknown!");
}

} // namespace
