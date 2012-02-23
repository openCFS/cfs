#include <assert.h>
#include <stddef.h>
#include <algorithm>
#include <iostream>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/programOptions.hh"
#include "General/Enum.hh"
#include "General/exception.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Optimizer/SnOpt.hh"
#include "Optimization/Optimizer/SnOptInterface.hh"
#include "Utils/Timer.hh"
#include "Utils/tools.hh"


// declare class specific logging stream
DECLARE_LOG(snopt)
DEFINE_LOG(snopt, "snopt")


namespace CoupledField
{

using std::endl;
using std::cout;

/** global pointer to the snopt class to be used by the callback function */
SnOpt* static_snopt = NULL;

/** this is a global function that can be called in the snopt-interface
 *  it simply forwards the call to the SnOpt->Callback function which has access to the class members */
int SnOpt_C_Callback(integer* Status, integer* n,
    doublereal* x, integer* needF, integer* nF, doublereal* F,
    integer* needG, integer* nG, doublereal* G,
    char* cu, integer* lencu, integer* iu, integer* leniu, doublereal* ru, integer* lenru)
{
  return static_snopt->Callback(Status, *n, x, needF, nF, F, needG, nG, G, cu, lencu, iu, leniu, ru, lenru);
}

SnOpt::SnOpt(Optimization* opt, PtrParamNode pn) :
  BaseOptimizer(opt, pn, Optimization::SNOPT_SOLVER),
  f_evals(0),       // number of function evaluations
  g_evals(0),       // number of gradient evaluations
  nxname(1),
  nFname(1),
  iSumm(6),         // variable for file output, 6 == stdout
  iPrint(0),        // variable for file output, 6 == stdout
  iSpecs(0),        // variable for file input
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
  EXIT(0),
  ObjRow(1),      // objective row, must be 1!!
  ObjAdd(0.0),
  optimizer_pn_(pn->Get(Optimization::optimizer.ToString(Optimization::SNOPT_SOLVER), ParamNode::PASS)),
  lin_constraints(0),
  nonlin_constraints(0),
  perform_commit_iteration_(false), // we want to perform commit iteration only after the first function eval.
  outfilename(""),
  timer_(new Timer()),
  timer_callback_(new Timer())
{
  LOG_TRACE(snopt) << "Initialize SnOpt";
  
  info_->Get(ParamNode::SUMMARY)->Get("snopt_timer")->SetValue(timer_ );
  info_->Get(ParamNode::SUMMARY)->Get("snopt_callback_timer")->SetValue(timer_callback_ );
  
  static_snopt = this;
  
  BaseOptimizer::PostInit(1.0);
  
  Init();
}

SnOpt::~SnOpt()
{
  if (iPrint != 0 && iPrint != 6)
    snclose_(&iPrint);

  if (iSpecs != 0)
    snclose_(&iSpecs);
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

  // set the file for snopt, formerly fort.1
  setSnoptOutputFiles();
  
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

void SnOpt::setSnoptOutputFiles()
{
  outfilename = progOpts->GetSimName() + ".snopt";

  if(optimizer_pn_ != NULL && optimizer_pn_->Has("output"))
    outfilename = optimizer_pn_->Get("output")->As<std::string>();
  
  // we enable command line output
  if(outfilename == "commandline")
  {
    iPrint = 6;
  }
  else
  {
    if(outfilename == "silent")
    {
      iPrint = 0;
    }
    else
    { 
      iPrint = 15;
      snopenappend_(&iPrint, (char*) outfilename.c_str(), &INFO, outfilename.size());
    }
  }

  INFO = 0;
  SetIntegerValue("print_file", iPrint);
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

  // this is needed for the call to snopt, but has no effect in our case
  integer npname(5);
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
  
  timer_->Start();
  
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
  
  timer_->Stop();
  
  InfoXMLOutput();
  
  optimization->CommitIteration();
}

void SnOpt::InfoXMLOutput()
{
  // this is the snopt specific outout
  info_->Get(ParamNode::SUMMARY)->Get("snopt_exit/exit")->SetValue(EXIT);
  info_->Get(ParamNode::SUMMARY)->Get("snopt_exit/info")->SetValue(INFO);
  info_->Get(ParamNode::SUMMARY)->Get("evaluations/f_evals")->SetValue(f_evals);
  info_->Get(ParamNode::SUMMARY)->Get("evaluations/g_evals")->SetValue(g_evals);
  info_->Get(ParamNode::SUMMARY)->Get("snopt_outputfile")->SetValue(outfilename);

  
  std::string exitstring;
  switch(INFO)
  {
  case 1:
    exitstring = "finished successfully - optimality conditions satisfied";
    break;
  case 3:
    exitstring = "finished successfully - requested accuracy could not be achieved";
    break;
  case 11:
    exitstring = "the problem appears to be infeasible - infeasible linear constraints";
    break;
  case 13:
    exitstring = "the problem appears to be infeasible - nonlinear infeasibilities minimized";
    break;
  case 21:
    exitstring = "unbounded objective";
    break;
  case 31:
    exitstring = "resource limit error - iteration limit reached";
    break;
  case 32:
    exitstring = "resource limit error - major iteration limit reached";
    break;
  case 33:
    exitstring = "resource limit error - the superbasics limit is too small";
    break;
  case 41:
    exitstring = "terminated after numerical difficulties - current point cannot be improved";
    break;
  case 42:
    exitstring = "terminated after numerical difficulties - singular basis";
    break;
  case 43:
    exitstring = "terminated after numerical difficulties - cannot satisfy the general constraints";
    break;
  case 44:
    exitstring = "terminated after numerical difficulties - ill-conditioned null-space basis";
    break;
  case 51:
    exitstring = "error in the user-supplied functions - incorrect objective derivatives";
    break;
  case 52:
    exitstring = "error in the user-supplied functions - incorrect constraint derivatives";
    break;
  case 71:
    exitstring = "user requested termination - terminated during function evaluation";
    break;
  default:
    exitstring = "not yet documented";
    break;
  }

  // this is the break node used by all optimizers
  PtrParamNode summary = optimization->optInfoNode->Get(ParamNode::SUMMARY);
  summary->Get("break/converged")->SetValue(INFO == 1 || INFO == 3 ? "yes" : "no");
  summary->Get("problem")->SetValue("SNOPT: " +exitstring);

  // the old and not compatible stuff. Now common for iTop
  //info_->Get(ParamNode::SUMMARY)->Get("snopt_exit/string")->SetValue(exitstring);
}

int SnOpt::Callback(integer* Status, const integer n,
    doublereal* x, integer* needF, integer* nF, doublereal* F,
    integer* needG, integer* nG, doublereal* G,
    char* cu, integer* lencu, integer* iu, integer* leniu, doublereal* ru, integer* lenru)
{
  timer_->Stop();
  timer_callback_->Start();
  
  LOG_DBG(snopt) << "Callback: needF=" << *needF << " needG=" << *needG << " x_avg = " << Average(x, n) << " x_std_dev = " << StandardDeviation(x, n);
  LOG_DBG3(snopt) << " x=" << StdVector<double>::ToString(n, x);

  // when the last call was a gradient eval we interpret this as major and do a commit
  // with the OLD design and it's function evaluations. But only if there was really a change between the last commit.
  // the special cases are the first iteration if setupLinearConstraints() had a feasible design or on the last commit.
  if(perform_commit_iteration_ && !optimization->GetDesign()->CompareDesign(x))
  {
    optimization->CommitIteration();
    perform_commit_iteration_ = false; // to be reset when enough functions evaluations have been done
  }

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
    
    // when the linear constraints are not fulfilled there will immediately be gradient evaluation as the
    // first snopt action and we do not want to loose the initial design in the output.
    perform_commit_iteration_ = true; // done on the next Callback() call
  }
  
  timer_callback_->Stop();
  timer_->Start();

  // flush the output
  cout << std::flush;
  
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
    LOG_DBG3(snopt) <<   "old grad_f[" << i << "] = " << grad_f[i]
                    << ", new grad_f[" << i << "] = " << gradhelper[i];
    grad_f[i] = gradhelper[i];
  }

  return result;
}

bool SnOpt::eval_g(int n, const double* x, int m, double* g)
{
  LOG_TRACE(snopt) << "eval_g";
  
  // for(int i = 0; i < m; i++) LOG_DBG3(snopt) <<   "old G[" << i << "] = " << g[i];

  EvalConstraints(n, x, m, true, g);
  
  for(int i = 0; i < m; i++) LOG_DBG3(snopt) <<   "g[" << i << "] = " << g[i];

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
    LOG_DBG3(snopt) <<   "grad_g[" << i << "] = " << pG[i];
  }

  return true;
}

void SnOpt::AdjustWorkArrayMemory()
{
  LOG_TRACE(snopt) << "AdjustWorkArrayMemory";
  LOG_TRACE(snopt) << "old values: lencw = " << lencw << ", leniw = " << leniw << ", lenrw = " << lenrw;
  
  // remember old values
  integer tmpcw(lencw);
  integer tmpiw(leniw);
  integer tmprw(lenrw);
  
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
  
  integer iPrt = 0;
  integer iSum = 0;
  // update lengths according to values obtained from snmema_
  if(lencw < mincw)
  {
    lencw = static_cast<int>(factor * static_cast<double>(mincw));
    cw.resize(8*lencw);
    snseti_("Total character workspace", &lencw, &iPrt, &iSum, &INFO,
            &cw[0], &tmpcw, &iw[0], &tmpiw, &rw[0], &tmprw, 25, 8*lencw);
    LOG_TRACE(snopt) << "new value: lencw = " << lencw << ", INFO = " << INFO;
  }
  
  if(leniw < miniw)
  {
    leniw = static_cast<int>(factor * static_cast<double>(miniw));
    iw.resize(leniw);
    snseti_("Total integer workspace", &leniw, &iPrt, &iSum, &INFO,
            &cw[0], &tmpcw, &iw[0], &tmpiw, &rw[0], &tmprw, 23, 8*lencw);
    LOG_TRACE(snopt) << "new value: leniw = " << leniw << ", INFO = " << INFO;
  }
  
  if(lenrw < minrw)
  {
    lenrw = static_cast<int>(factor * static_cast<double>(minrw));
    rw.resize(lenrw);
    snseti_("Total real workspace", &lenrw, &iPrt, &iSum, &INFO,
            &cw[0], &tmpcw, &iw[0], &tmpiw, &rw[0], &tmprw, 20, 8*lencw);
    LOG_TRACE(snopt) << "new value: lenrw = " << lenrw << ", INFO = " << INFO;
  }
}

void SnOpt::SetSnOptOptions()
{
  LOG_TRACE(snopt) << "SetSnOptOptions";
  
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
        setMinorItLimit = true;
  
      
      if(name == "iterations_limit")
        setItLimit = true;
      
      SetIntegerValue(name, value);
    }
    
    list = optimizer_pn_->GetListByVal("option", "type", "real");
    
    for(unsigned int i = 0; i < list.GetSize(); i++)
      SetNumericValue(list[i]->Get("name")->As<std::string>(), list[i]->Get("value")->As<double>());
    
    list = optimizer_pn_->GetListByVal("option", "type", "string");
    
    for(unsigned int i = 0; i < list.GetSize(); i++)
      SetStringValue(list[i]->Get("name")->As<std::string>(), list[i]->Get("value")->As<std::string>());
  }
  
  if(!setMinorItLimit) SetIntegerValue("minor_iterations_limit", 5000);
  if(!setItLimit) SetIntegerValue("iterations_limit", 1000000);
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
    const unsigned int patternsize(pattern.GetSize());
    assert(patternsize > 0);

    if(g->IsLinear()) // linear
    {
      for(unsigned int e = 0; e < patternsize; ++e)
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
      for(unsigned int e = 0; e < patternsize; ++e)
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
  
  optimization->constraints.view->Done();

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
  
  // we have just evaluated the state problem, commit the (initial) iteration as snopt might
  // evaluate in its's callback already an updated  design
  perform_commit_iteration_ = true;
  // optimization->CommitIteration();

  for(int i = 0; i < nA; i++)
  {
    A[i] = lincon[i];
    LOG_DBG3(snopt) << "A[" << i << "] = " << A[i];
  }
  
  //assert(count == lin_constraints);
}

void SnOpt::SetIntegerValue(const std::string& key, integer value)
{
  string option;
  
  if(key == "major_iterations_limit")
  {
    option = "Major iterations limit";
  }
  else if(key == "print_file")
  {
    option = "Print file";
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
    EXCEPTION("trying to set an snopt-option (int) that is unknown!");
}

void SnOpt::SetNumericValue(const std::string& key, double value)
{  
  string option;
  
  if(key == "major_optimality_tolerance")
  {
    option = "Major optimality tolerance";
  }
  else if(key == "minor_optimality_tolerance")
  {
    option = "Minor optimality tolerance";
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
  else if(key == "difference_interval")
  {
    option = "Difference interval";
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
    EXCEPTION("trying to set an snopt-option (real) that is unknown!");
}

void SnOpt::SetStringValue(const std::string& key, const std::string& value)
{
  string option;

  if(key == "qpsolver")
  {
    option = "QPSolver ";
    option.append(value);
  }
  
  if(!option.empty())
  {
    LOG_TRACE(snopt) << "adjusted " << key;
    // set the option
    snset_(option.c_str(), &iPrint, &iSumm, &INFO,
        &cw[0], &lencw, &iw[0], &leniw, &rw[0], &lenrw, option.size(), 8*lencw );

    // check INFO, if after setting an option INFO > 0 we have had an error!
    assert(INFO == 0);
  }
  else // nothing found
    EXCEPTION("trying to set an snopt-option (string) that is unknown!");
}

} // namespace
