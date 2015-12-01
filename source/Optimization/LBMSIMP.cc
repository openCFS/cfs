#include "Optimization/LBMSIMP.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Domain/Domain.hh"
#include "PDE/LatticeBoltzmannPDE.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"

DECLARE_LOG(lbmsimp)
DEFINE_LOG(lbmsimp, "lbmsimp")

LBMSIMP::LBMSIMP()
{
	lbm = NULL;
}

LBMSIMP::~LBMSIMP()
{

}

void LBMSIMP::SolveStateProblem(Excitation* ev_only_excite)
{
  assert(context->pde != NULL);
  lbm = dynamic_cast<LatticeBoltzmannPDE*>(context->pde);
  assert(lbm != NULL);
  LOG_DBG(lbmsimp) << "SSP -> solve";
  lbm->Solve();
}

/** overloads SIMP::CalcFunction()
 * @see ErsatzMaterial::CalcFunction */
double LBMSIMP::CalcFunction(Excitation& excite, Function* f, bool derivative)
{
  LOG_DBG(lbmsimp) << "CF -> f=" << f->ToString() << " d=" << derivative;

  // assume the problem is solved for the current design by DesignChanged()

  switch(f->GetType())
  {
  case Function::PRESSURE_DROP:
  {
    if(!derivative)
      return lbm->CalcPressureDrop();
    else
    {
      switch(lbm->GetIface())
      {
      case LatticeBoltzmannPDE::EXTERNAL:
      case LatticeBoltzmannPDE::INTERNAL:
//    	std::cout << "size of f->elements:" <<  f->elements.GetSize() << std::endl;
        lbm->SensitivityAnalysis(design->GetTransferFunction(f->elements[0]), f, design);
        break;
      }
      return 0.0;
    }
  }
  break;

  default: // return below as we don't implement
    break;
  }

  return SIMP::CalcFunction(excite, f, derivative);
}

