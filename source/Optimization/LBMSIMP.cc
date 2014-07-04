#include "Optimization/LBMSIMP.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Domain/domain.hh"
#include "PDE/LatticeBoltzmannPDE.hh"



LBMSIMP::LBMSIMP()
{
  lbm = dynamic_cast<LatticeBoltzmannPDE*>(pde);

  for(unsigned int r = 0; r < design->GetRegionIds().GetSize(); r++)
    GetForm(design->GetRegionIds()[r], lbm, lbm, "LatticeBoltzmannInt")->SetSolDependent(true);

}

LBMSIMP::~LBMSIMP()
{

}

void LBMSIMP::SolveStateProblem(Excitation* ev_only_excite)
{
  lbm->Solve();
}

/** overloads SIMP::CalcFunction()
 * @see ErsatzMaterial::CalcFunction */
double LBMSIMP::CalcFunction(Excitation& excite, Function* f, bool derivative)
{
  switch(f->GetType())
  {
  case Function::PRESSURE_DROP:
  {
    if(!derivative)
      return lbm->CalcPressureDrop();
    else
    {
      CalcPressureDropDerivative(f);
      return 0.0;
    }
  }
  break;

  default: // return below as we don't implement
    break;
  }

  return SIMP::CalcFunction(excite, f, derivative);
}

void LBMSIMP::CalcPressureDropDerivative(Function* f)
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
}


void LBMSIMP::ConstructAdjointRHS(Excitation& excite, Function* f)
{
}


/** This is our part for CalcU1KU2() -> This set the matrix derivatives form ELEC and
 * PIEZO_COUPLING ( + transposed) */
void LBMSIMP::SetElementK(DesignElement* de, const TransferFunction* tf, Application app, DenseMatrix* out, CalcMode calcMode, bool derivative)
{
}

