#include <iostream>

#include "Domain/Domain.hh"
#include "Driver/BaseDriver.hh"
#include "Driver/HarmonicDriver.hh"
#include "General/Enum.hh"
#include "General/Environment.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Optimizer/EvaluateOnly.hh"
#include "Utils/StdVector.hh"
#include "Utils/Timer.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"


DECLARE_LOG(eval)
DEFINE_LOG(eval, "eval")

using namespace CoupledField;

EvaluateOnly::EvaluateOnly(Optimization* optimization, PtrParamNode pn)
 : BaseOptimizer(optimization, pn, Optimization::EVALUATE_INITIAL_DESIGN)
{
  // reduce to our actual ParamNode
  pn = pn->Get(Optimization::optimizer.ToString(Optimization::EVALUATE_INITIAL_DESIGN), ParamNode::PASS);
  
  eval_grad = true;
  if(pn != NULL){ // tag can be omitted
    eval_grad = pn->Get("objective_gradient")->As<bool>();
  }

  PostInitScale(1.0, true);
}

void EvaluateOnly::SolveProblem()
{
  // solve the state problem with the initial guess.
  std::cout << "Evaluate state problem for initial guess ..." << std::endl;

  // there is no time within this optimizer spent
  optimizer_timer_->Stop();

  // in the harmonic case we sweep over multiple frequencies if we have not "multipleExcitation"
  HarmonicDriver* hd = Optimization::context->GetHarmonicDriver();
  int end = optimization->context->IsHarmonic() && !optimization->GetMultipleExcitation()->IsEnabled() ? hd->freqs.GetSize() : 1;

  // space to store the gradient values, we need it to evaluate density filtering.
  StdVector<double> grad(optimization->GetDesign()->GetNumberOfVariables());
  grad.Init(0.0);
  // scale the window to the whole data domain
  grad.window.Set(grad);

  // our initial design
  StdVector<double> x;
  optimization->GetDesign()->WriteDesignToExtern(x);

  for(int i = 0; i < end; i++)
  {
    // we do multiple excitations only when end > 1. Otherwise it is unused
    if(end > 1)
    {
      Excitation& excite = optimization->GetMultipleExcitation()->excitations[0];
      excite.index = 0; // we solve always at the same position
      excite.f_link = &hd->freqs[i];
      excite.frequency = excite.f_link->freq;
    }

    // special case only in harmonic case with more frequencies but not multiple_loads optimization
    optimization->SolveStateProblem();

    double v = optimization->CalcObjective();
    LOG_DBG(eval) << "SP: obj=" << v;
    // calc gradients, they might be stored in store results!
    // gradients might need adjoints
    if(eval_grad){
      optimization->SolveAdjointProblems();
      optimization->CalcObjectiveGradient(&grad);
      for(unsigned int i = 0; i < grad.GetSize(); i++) {
        BaseDesignElement* de = optimization->GetDesign()->GetDesignElement(i);
        LOG_DBG2(eval) << "SP: obj grad i=" << i << " (" << (i+1) <<  ") de=\"" << de->ToString() << "\" -> " << grad[i];
      }
    }
    
    for(int c = 0; c < optimization->constraints.view->GetNumberOfTotalConstraints(); c++)
    {
      Condition* g = optimization->constraints.view->Get(c);
      v = EvalConstraint(g, false, false);

      double scaling = g->DoObjectiveScaling() ? objective->scaling.value : g->manual_scaling_value;

      LOG_DBG(eval) << "SP: g[" << c << " (" << (c+2) << ")]=" << g->ToString() << " -> " << v * scaling; // snopt index in brackets

      if(!g->IsObservation()) // not for observation stuff
      {
        StdVector<unsigned int>& pattern = g->GetSparsityPattern();
        grad.window.Set(0, pattern.GetSize()); // necessary for a local condition assert
        optimization->CalcConstraintGradient(g, &grad);
        for(unsigned int i = 0; i < pattern.GetSize(); i++) {
          BaseDesignElement* de = optimization->GetDesign()->GetDesignElement(pattern[i]);
          LOG_DBG2(eval) << "SP: grad g[" << c << " (" << (c+2) << ")]=" << g->ToString() << " i=" << i
                         << "(" << (i+1) << ") pi=" << pattern[i] << "(" << (pattern[i]+1) <<  ") de=\"" << de->ToString() << "\" -> " << grad[i] * scaling;
        }
      }
    }
    optimization->constraints.view->Done(); // reset the slope constraints to global

    // multiple excitations in evaluate only are identified as increasing "iterations".
    optimization->CommitIteration();
  }
}
