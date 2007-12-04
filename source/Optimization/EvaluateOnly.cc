#include "Optimization/EvaluateOnly.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/SIMP.hh"
#include "Optimization/PiezoSIMP.hh"

using namespace CoupledField;

EvaluateOnly::EvaluateOnly(Optimization* optimization, ParamNode* pn)
 : BaseOptimizer(optimization, pn)
{
  // reduce to our actual ParamNode
  pn = pn->Get(Optimization::optimizer.ToString(Optimization::EVALUATE_INITIAL_DESIGN), false);
  
  PostInit(false);  
}

void EvaluateOnly::SolveProblem()
{
  // solve the state problem with the initial guess.
  std::cout << "Evaluate state problem for initial guess ..." << std::endl;
  // note that when PiezoSIMP and storage is not PiezoSIMP::COMMIT we might have wrong 
  // special results before CalcObjective()
  optimization->SolveStateProblem();
  std::cout << "objective: " << optimization->CalcObjective() << std::endl;
  // see note above!
  if(optimization->GetOptimizationType() == Optimization::SIMP_TYPE 
     && (dynamic_cast<SIMP*>(optimization))->GetSystem() == SIMP::PIEZO
     && (dynamic_cast<PiezoSIMP*>(optimization))->GetStorage() != PiezoSIMP::COMMIT)
    optimization->SolveStateProblem();
  
  // calc gradients, they might be stored in store results!
  optimization->CalcObjectiveGradient(NULL);

  StdVector<Condition>& cns = optimization->constraints;   
  /** The "inactive" constraits with output_only mode in xml */
  StdVector<Condition>& ops = optimization->outputs;   
  
  for(unsigned int i = 0; i < cns.GetSize(); i++)
  { 
     std::cout << "constraint " << cns[i].ToString() << ": " 
               << optimization->CalcConstraint(&cns[i]) << std::endl;
     optimization->CalcConstraintGradient(&cns[i], NULL);          
  }

  for(unsigned int i = 0; i < ops.GetSize(); i++)
  { 
     std::cout << "observation " << ops[i].ToString() << ": "
               << optimization->CalcConstraint(&ops[i]) << std::endl;
  }
  optimization->CommitIteration();
}
