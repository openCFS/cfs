#include "Optimization/EvaluateOnly.hh"
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

  for(int i = 0; i < end; i++)
  {
    // we do multiple excitations only when end > 1. Otherwise it is unused
    if(end > 1)
    {
      Excitation& excite = dynamic_cast<ErsatzMaterial*>(optimization)->excitations[0];
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
    optimization->CalcObjectiveGradient(NULL);
    
    for(int c = 0; c < optimization->constraints.view->GetNumberOfTotalConstraints(); c++)
    {
      Condition* g = optimization->constraints.view->Get(c);
      optimization->CalcConstraint(g);
      if(g->IsActive()) // not for observation stuff
        optimization->CalcConstraintGradient(g);
    }
    optimization->constraints.view->Done(); // reset the slope constraints to global

    // multiple excitations in evaluate only are identified as increasing "iterations".
    optimization->CommitIteration(false);
  }
}
