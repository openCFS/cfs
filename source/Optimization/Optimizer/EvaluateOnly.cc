#include "Optimization/Optimizer/EvaluateOnly.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Optimization.hh"
#include "Driver/harmonicDriver.hh"
#include "Domain/domain.hh"
#include "Utils/Timer.hh"

using namespace CoupledField;

EvaluateOnly::EvaluateOnly(Optimization* optimization, PtrParamNode pn)
 : BaseOptimizer(optimization, pn)
{
  // reduce to our actual ParamNode
  pn = pn->Get(Optimization::optimizer.ToString(Optimization::EVALUATE_INITIAL_DESIGN), ParamNode::PASS);

  PostInit(1.0, true);
}

void EvaluateOnly::SolveProblem()
{
  // solve the state problem with the initial guess.
  std::cout << "Evaluate state problem for initial guess ..." << std::endl;

  // there is no time within this optimizer spent
  timer_->Stop();

  // in the harmonic case we sweep over multiple frequencies if we have not "multipleExcitation"
  HarmonicDriver* hd = dynamic_cast<HarmonicDriver*>(domain->GetDriver());
  int end = optimization->IsHarmonic() && !optimization->GetMultipleExcitation()->IsEnabled() ? hd->freqs.GetSize() : 1;

  // to store the gradient values, we need it to evaluate density filtering
  StdVector<double> grad(optimization->GetDesign()->GetNumberOfVariables());
  grad.Init(0.0);
  // scale the window to the whole data domain
  grad.window.Set(grad);

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

    optimization->CalcObjective();
    // calc gradients, they might be stored in store results!
    // gradients might need adjoints
    optimization->SolveAdjointProblems();
    optimization->CalcObjectiveGradient(&grad);
    
    for(int c = 0; c < optimization->constraints.view->GetNumberOfTotalConstraints(); c++)
    {
      Condition* g = optimization->constraints.view->Get(c);
      optimization->CalcConstraint(g);
      if(!g->IsObservation()) // not for observation stuff
        optimization->CalcConstraintGradient(g, &grad);
    }
    optimization->constraints.view->Done(); // reset the slope constraints to global

    // multiple excitations in evaluate only are identified as increasing "iterations".
    optimization->CommitIteration(false);
  }
}
