#include "Optimization/LBMSIMP.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Domain/Domain.hh"
// FIXME #include "PDE/LatticeBoltzmannPDE.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"

DECLARE_LOG(simp)

LBMSIMP::LBMSIMP()
{
  /* FIXME
  lbm = dynamic_cast<LatticeBoltzmannPDE*>(pde);

  for(unsigned int r = 0; r < design->GetRegionIds().GetSize(); r++)
    GetForm(design->GetRegionIds()[r], lbm, lbm, "LatticeBoltzmannInt")->SetSolDependent(true);
  */
}

LBMSIMP::~LBMSIMP()
{

}

void LBMSIMP::SolveStateProblem(Excitation* ev_only_excite)
{
  LOG_DBG(simp) << "SSP -> solve";
  // FIXME lbm->Solve();

}

/** overloads SIMP::CalcFunction()
 * @see ErsatzMaterial::CalcFunction */
double LBMSIMP::CalcFunction(Excitation& excite, Function* f, bool derivative)
{
  LOG_DBG(simp) << "CF -> f=" << f->ToString() << " d=" << derivative;

  // assume the problem is solved for the current design by DesignChanged()

  switch(f->GetType())
  {
  case Function::PRESSURE_DROP:
  {
    /* FIXME
    if(!derivative)
      return lbm->CalcPressureDrop();
    else
    {
      switch(lbm->GetIface())
      {
      case LatticeBoltzmannPDE::EXT_MATLAB:
        lbm->SetPrecalculatedGradient(f->elements, f);
        break;

      case LatticeBoltzmannPDE::EXT_CFSxLBM:
      case LatticeBoltzmannPDE::INTERNAL:
        lbm->SensitivityAnalysis(design->GetTransferFunction(f->elements[0]), f, design);
        break;
      }
      return 0.0;
    }
    */
  }
  break;

  default: // return below as we don't implement
    break;
  }

  return SIMP::CalcFunction(excite, f, derivative);
}

