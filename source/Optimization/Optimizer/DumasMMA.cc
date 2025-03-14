#include "Optimization/Optimizer/DumasMMA.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "General/Exception.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include <MMASolver.h>
#include <GCMMASolver.h>

using namespace CoupledField;

// declare class specific logging stream
DEFINE_LOG(dumas, "dumas")

DumasMMA::DumasMMA(Optimization* optimization, PtrParamNode pn, Optimization::Optimizer dumas) : BaseOptimizer(optimization, pn, dumas)
{
  assert(dumas == Optimization::DUMAS_MMA || dumas == Optimization::DUMAS_GCMMA);

  if(this_opt_pn_)
  {
    asyminit = this_opt_pn_->Get("asyminit")->As<double>();
    asymdec  = this_opt_pn_->Get("asymdec")->As<double>();
    asyminc  = this_opt_pn_->Get("asyminc")->As<double>();

    ai_lin_aux_z = this_opt_pn_->Get("ai_lin_aux_z")->As<double>();
    c_lin_aux_y  = this_opt_pn_->Get("c_lin_aux_y")->As<double>();
    d_quad_aux_y = this_opt_pn_->Get("d_quad_aux_y")->As<double>();

    if(this_opt_pn_->Get("move_limit")->As<string>() == "tuned")
    {
      tuned = new Tuned(this_opt_pn_->Get("tuned", ParamNode::PASS), &move_limit, 1, 5, this);
      move_limit = tuned->max; // set by Tuned
    }
    else
      move_limit = this_opt_pn_->Get("move_limit")->As<double>();
  }
  Optimization* opt = this->optimization;

  int m = (int) opt->constraints.view->GetNumberOfActiveConstraints();
  int n = (int) opt->design->GetNumberOfVariables();
  assert(n >= (int) opt->design->data.GetSize()); // there might be auxiliary variables

  xval.Resize(n);
  dfdx.Resize(n);
  g.Resize(m);
  dgdx.Resize(m*n);
  xmin.Resize(n);
  xmax.Resize(n);

  if(m == 0)
    throw Exception("dumas requires at least one constraint");

  if(dumas == Optimization::DUMAS_MMA)
  {
    mma = new MMASolver(n,m, ai_lin_aux_z, c_lin_aux_y, d_quad_aux_y);
    mma->SetAsymptotes(asyminit, asymdec, asyminc);
  }
  else
  {
    gcmma = new GCMMASolver(n,m, ai_lin_aux_z, c_lin_aux_y, d_quad_aux_y);
    gcmma->SetAsymptotes(asyminit, asymdec, asyminc);
  }

  optimizer_timer_->Stop(); // we have no own PostInit()
}

DumasMMA::~DumasMMA()
{
  delete mma;
  delete gcmma;
}

void DumasMMA::PostInit()
{
  PostInitScale(1.0); // Overwritten if scale is set in the xml. The target for Niels is 10
  BaseOptimizer::PostInit();
}

void DumasMMA::ToInfo(PtrParamNode pn)
{
  BaseOptimizer::ToInfo(pn);

  pn->Get("move_limit")->SetValue(tuned ? "tuned" : std::to_string(move_limit));
  if(tuned)
    tuned->ToInfo(pn->Get("tuned_move_limit"));
  pn->Get("asyminit")->SetValue(asyminit);
  pn->Get("asymdec")->SetValue(asymdec);
  pn->Get("asyminc")->SetValue(asyminc);
  pn->Get("a0_lin_aux_z")->SetValue(1.0);
  pn->Get("ai_lin_aux_z")->SetValue(ai_lin_aux_z);
  pn->Get("c_lin_aux_y")->SetValue(c_lin_aux_y);
  pn->Get("d_quad_aux_y")->SetValue(d_quad_aux_y);
}

void DumasMMA::LogFileHeader(Optimization::Log& log)
{
  BaseOptimizer::LogFileHeader(log);

  if(gcmma)
    log.AddToHeader("gc_search");
}

void DumasMMA::LogFileLine(std::ofstream* out, PtrParamNode iteration)
{
  BaseOptimizer::LogFileLine(out, iteration);

  if(gcmma)
  {
    if(out) *out << " \t" << gc_search;
    iteration->Get("gc_search")->SetValue(gc_search);
  }
}

void DumasMMA::DescribeProperties(StdVector<std::pair<std::string, std::string> >& map) const
{
  map.Push_back(std::make_pair("move_limit", std::to_string(move_limit)));
  map.Push_back(std::make_pair("asymdec", std::to_string(asymdec)));
  map.Push_back(std::make_pair("asyminc", std::to_string(asyminc)));
}


void DumasMMA::SolveProblem()
{
  Optimization* opt = this->optimization;
  unsigned int m = (int) opt->constraints.view->GetNumberOfActiveConstraints();
  unsigned int n = (int) opt->design->GetNumberOfVariables();
  assert(xval.GetSize() == n);
  assert(n >= opt->design->data.GetSize());

  // design and f_i for gcmma
  Vector<double> xnew(gcmma ? n : 0);
  Vector<double> gnew(gcmma ? m : 0);

  // read initial design
  opt->design->WriteDesignToExtern(xval.GetPointer());

  // start with iteration 0 which is the initial design
  assert(optimization->GetCurrentIteration() == 0);

  // run at least once such that we have observes and constraints for the stopping criteria evaluated
  do
  {
    // calc gradients to store the results in data[element]...
    // only the gradients are needed for the calculation of the next iteration
    double f = EvalObjective(n, xval.GetPointer(), true); // only for gcmma and logging/stopping
    EvalGradObjective(n,  xval.GetPointer(), true, dfdx); // we scale
    EvalConstraints(n, xval.GetPointer(), m, true, g.GetPointer(), true); // normalize f_i <= 0, we do not give the constraint bound to the solver
    EvalGradConstraints(n, xval.GetPointer(), m, n*m, true, true, dgdx); // normalize

    LOG_DBG(dumas) << "SP: it=" << optimization->GetCurrentIteration() << " f=" << f << " gcmma=" << (gcmma != nullptr);
    // store iteration 0
    if(optimization->GetCurrentIteration() == 0)
    {
      optimizer_timer_->Stop();
      optimization->SolveStateProblem();
      optimizer_timer_->Start();
    }

    // Set outer move limits
    for(unsigned int i=0; i<n; i++)
    {
      BaseDesignElement* de = opt->design->GetDesignElement(i);
      xmax[i] = std::min(de->GetUpperBound(), xval[i] + move_limit);
      xmin[i] = std::max(de->GetLowerBound(), xval[i] - move_limit);
    }

    assert(dfdx.GetSize() > 0);
    assert(g.GetSize() > 0 && dgdx.GetSize() > 0);

    if(mma)
      mma->Update(xval.GetPointer(),dfdx.GetPointer(),g.GetPointer(),dgdx.GetPointer(),xmin.GetPointer(),xmax.GetPointer());
    else
    {
      // see https://github.com/jdumas/mma/blob/master/tests/toy.cpp
      assert(xnew.GetSize() == n); // for gcmma only
      gcmma->OuterUpdate(xnew.GetPointer(), xval.GetPointer(), f, dfdx.GetPointer(), g.GetPointer(), dgdx.GetPointer(), xmin.GetPointer(), xmax.GetPointer());

      // test xnew
      opt->design->ReadDesignFromExtern(xnew, false); // don't write this test design to density.xml
      double fnew = EvalObjective(n, xnew.GetPointer(), true);
      EvalConstraints(n, xnew.GetPointer(), m, true, gnew.GetPointer(), true);

      bool conserv = gcmma->ConCheck(fnew, gnew.GetPointer());
      for(gc_search = 1; !conserv && gc_search < 15; gc_search++) // we searched one line above :)
      {
        // Inner iteration update
        gcmma->InnerUpdate(xnew.GetPointer(), fnew, gnew.GetPointer(), xval.GetPointer(), f,
          dfdx.GetPointer(), g.GetPointer(), dgdx.GetPointer(), xmin.GetPointer(), xmax.GetPointer());

        opt->design->ReadDesignFromExtern(xnew, false); // don't write this test design to density.xml
        fnew = EvalObjective(n, xnew.GetPointer(), true);
        EvalConstraints(n, xnew.GetPointer(), m, true, gnew.GetPointer(), true);

        conserv = gcmma->ConCheck(fnew, gnew.GetPointer());
      }
      xval = xnew;
    }
    optimizer_timer_->Stop();

    opt->design->ReadDesignFromExtern(xval, true);
    CommitIteration();
    optimizer_timer_->Start();
  }
  while(!optimization->DoStopOptimization() && optimization->GetCurrentIteration() <= optimization->GetMaxIterations());
  optimizer_timer_->Stop();

  PtrParamNode in = optimization->optInfoNode->Get(ParamNode::SUMMARY)->Get("break");

  if(optimization->GetCurrentIteration() >= optimization->GetMaxIterations()-1)
  {
    optimization->DoStopOptimizationHelper(false, "Maximum iterations exceeded");
    std::cout << "++ max iterations reached" << std::endl;
  }
}
