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

  opt_timer->Stop();
  PostInitScale(1.0, true);
}

void EvaluateOnly::SolveProblem()
{
  opt_timer->Stop(); // we don't need this time

  // solve the state problem with the initial guess.
  std::cout << "Evaluate state problem for initial guess ..." << std::endl;

  StdVector<double> xl(optimization->GetDesign()->GetNumberOfVariables());
  StdVector<double> xu(xl.GetSize());
  StdVector<double> gl(optimization->constraints.view->GetNumberOfActiveConstraints());
  StdVector<double> gu(gl.GetSize());
  GetBounds(xl.GetSize(), xl.GetPointer(), xu.GetPointer(), gl.GetSize(), gl.GetPointer(), gu.GetPointer());

#ifndef NDEBUG
  for(int i = 0; i < optimization->constraints.view->GetNumberOfActiveConstraints(); i++)
  {
    Condition* g = optimization->constraints.view->Get(i);
    LOG_DBG(eval) << "SP: bnds g[" << i << " (" << (g->GetIndex()+1) << ")]=" << g->ToString() << " -> " << gl[i] << " ... " << gu[i];
  }
  optimization->constraints.view->Done();
#endif

  // in the harmonic case we sweep over multiple frequencies if we have not "multipleExcitation"
  HarmonicDriver* hd = Optimization::context->GetHarmonicDriver();
  // end is > 1 for "multiple_excitations" set to false in order to evaluate the functions separately for each frequency
  int end = optimization->context->IsHarmonic() && !optimization->GetMultipleExcitation()->IsEnabled() ? hd->freqs.GetSize() : 1;

  // space to store the gradient values, we need it to evaluate sensitivity filtering, special results and debugging.
  StdVector<double> grad(optimization->GetDesign()->GetNumberOfVariables());
  grad.Init(0.0);
  // scale the window to the whole data domain
  grad.window.Set(grad);

  // our initial design
  StdVector<double> x;
  optimization->GetDesign()->WriteDesignToExtern(x);

  // if we really loop here in evaluate we want to show the function values separately for each frequency
  for(int i = 0; i < end; i++)
  {
    Excitation* excite = NULL;
    // we do multiple excitations only when end > 1. Otherwise it is unused
    // end = 1 means that the functions are evaluated as weighted sums ("multiple_excitation=true")
    if(end > 1)
    {
      excite = &(optimization->GetMultipleExcitation()->excitations[0]);
      excite->index = 0; // we solve always at the same position
      excite->f_link = &hd->freqs[i];
      excite->frequency = excite->f_link->freq;
    }

    // special case only in harmonic case with more frequencies but not multiple_loads optimization
    optimization->SolveStateProblem(excite);

    eval_obj_timer_->Start();
    design_.value = optimization->CalcObjective(excite);
    eval_obj_timer_->Stop();
    LOG_DBG(eval) << "SP: obj=" << design_.value;
    // calc gradients, they might be stored in store results!
    // gradients might need adjoints
    optimization->SolveAdjointProblems(excite);

    if(eval_grad){
      eval_grad_obj_timer_->Start();
      optimization->CalcObjectiveGradient(&grad, excite);
      eval_grad_obj_timer_->Stop();
      for(unsigned int i = 0; i < grad.GetSize(); i++) {
        BaseDesignElement* de = optimization->GetDesign()->GetDesignElement(i);
        LOG_DBG2(eval) << "SP: obj grad i=" << i << " (" << (i+1) <<  ") de=\"" << de->ToString() << "\" -> " << std::setprecision(10) << grad[i];
      }
    }

    for(int c = 0; c < optimization->constraints.view->GetNumberOfTotalConstraints(); c++)
    {
      Condition* g = optimization->constraints.view->Get(c);
      opt_timer->Start(); // only for the assert
      double v = EvalConstraint(g, false, false, true, excite); // sets the timer itself
      opt_timer->Stop();

      double scaling = g->DoObjectiveScaling() ? objective->scaling.value : g->manual_scaling_value;

      LOG_DBG(eval) << "SP: g[" << c << " (" << (c+2) << ")]=" << g->ToString() << " -> " << std::setprecision(10) << v * scaling; // snopt index in brackets

      if(!g->IsObservation()) // not for observation stuff
      {
        StdVector<unsigned int>& pattern = g->GetSparsityPattern();
        grad.window.Set(0, pattern.GetSize()); // necessary for a local condition assert
        eval_grad_const_timer_->Start();
        optimization->CalcConstraintGradient(g, &grad, excite);
        eval_grad_const_timer_->Stop();
#ifndef NDEBUG
        for(unsigned int i = 0; i < pattern.GetSize(); i++) {
          BaseDesignElement* de = optimization->GetDesign()->GetDesignElement(pattern[i]);
          LOG_DBG2(eval) << "SP: grad g[" << c << " (" << (c+2) << ")]=" << g->ToString() << " i=" << i
                         << "(" << (i+1) << ") pi=" << pattern[i] << "(" << (pattern[i]+1) <<  ") de=\"" << de->ToString()
                         << "\" -> " << std::setprecision(10) << grad[i] * scaling;
        }
#endif
      }
    }
    optimization->constraints.view->Done(); // reset the slope constraints to global


    // multiple excitations in evaluate only are identified as increasing "iterations".
    CommitIteration();
  }
}
