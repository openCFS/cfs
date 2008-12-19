#include "Optimization/EvaluateOnly.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/SIMP.hh"
#include "Optimization/PiezoSIMP.hh"
#include "Driver/harmonicDriver.hh"
#include "Domain/domain.hh"

using namespace CoupledField;

EvaluateOnly::EvaluateOnly(Optimization* optimization, ParamNode* pn)
 : BaseOptimizer(optimization, pn)
{
  // reduce to our actual ParamNode
  pn = pn->Get(Optimization::optimizer.ToString(Optimization::EVALUATE_INITIAL_DESIGN), false);

  PostInit(1.0, true);
}

void EvaluateOnly::SolveProblem()
{
  // solve the state problem with the initial guess.
  std::cout << "Evaluate state problem for initial guess ..." << std::endl;

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
    optimization->CalcObjectiveGradient(NULL);

    StdVector<Condition>& cns = optimization->constraints;
    StdVector<Condition>& ops = optimization->outputs; // The "inactive" constraits with output_only mode in xml

    for(unsigned int c = 0; c < cns.GetSize(); c++)
      optimization->CalcConstraintGradient(&cns[c], NULL);

    for(unsigned int c = 0; c < ops.GetSize(); c++)
      optimization->CalcConstraintGradient(&ops[c], NULL);

    // multiple excitations in evaluate only are identified as increasing "iterations".
    optimization->CommitIteration(false);
  }
}
