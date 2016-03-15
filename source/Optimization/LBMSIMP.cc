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

}

LBMSIMP::~LBMSIMP()
{

}

void LBMSIMP::SolveLBMState(Excitation* ev_only_excite)
{
  assert(context->pde != NULL);
  std::cout << "num pdes: " << context->pdes.size() << std::endl;
  std::cout << context->pdes[App::LBM]->GetName() << std::endl;
  LatticeBoltzmannPDE* lbm = dynamic_cast<LatticeBoltzmannPDE*>(context->pde);
  assert(lbm != NULL);
  LOG_DBG(lbmsimp) << "SSP -> solve";
  lbm->Solve();
}

/** overloads SIMP::CalcFunction()
 * @see ErsatzMaterial::CalcFunction */
double LBMSIMP::CalcFunction(Excitation& excite, Function* f, bool derivative)
{
  LOG_DBG(lbmsimp) << "CF -> f=" << f->ToString() << " d=" << derivative;
  LatticeBoltzmannPDE* lbm = dynamic_cast<LatticeBoltzmannPDE*>(context->pde);

  // assume the problem is solved for the current design by DesignChanged()
  if(!derivative)
  {
    switch(f->GetType())
      {
        case Function::PRESSURE_DROP:
          return lbm->CalcPressureDrop();
        default: // return below as we don't implement
          return SIMP::CalcFunction(excite, f, derivative);
      }

  }
  else
  {
    std::cout << "function: " << Function::type.ToString(f->GetType()) << std::endl;
    if(f->GetType() == Function::PRESSURE_DROP)
    {
      assert(lbm->GetIface() != LatticeBoltzmannPDE::EXTERNAL);
      lbm->SensitivityAnalysis(design->GetTransferFunction(f->elements[0]), f, design);
      return 0.0;
    }
    else
      return SIMP::CalcFunction(excite, f, derivative);
  }
}
