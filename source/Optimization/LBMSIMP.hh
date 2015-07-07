#ifndef LBMSIMP_HH_
#define LBMSIMP_HH_

#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Function.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/SIMP.hh"

namespace CoupledField
{
class DenseMatrix;
class DesignElement;
class Excitation;
class TransferFunction;
}  // namespace CoupledField

namespace CoupledField
{
class LatticeBoltzmannPDE;

/** Extension for LBM optimization based on an non-FE solver */
class LBMSIMP : public SIMP
{
public:
  LBMSIMP();

  ~LBMSIMP();

protected:

  /** overloads SIMP::CalcFunction()
   * @see ErsatzMaterial::CalcFunction */
  double CalcFunction(Excitation& excite, Function* f, bool derivative);

  /** overwrite. We have no FE and don't store the results, all done by LatticeBoltzmannPDE. */
  void SolveStateProblem(Excitation* ev_only_excite = NULL);

private:

  /** shortcut to our pde, is also in ErsatzMaterial::pdes */
  LatticeBoltzmannPDE* lbm;
};

} // end of namespace




#endif /* LBMSIMP_HH_ */
