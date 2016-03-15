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

/** Extension for App::LBM optimization based on an non-FE solver */
class LBMSIMP : public SIMP
{
public:
  LBMSIMP();

  ~LBMSIMP();

  void SolveLBMState(Excitation* ev_only_excite = NULL);

protected:

  /** overloads SIMP::CalcFunction()
   * @see ErsatzMaterial::CalcFunction */
  virtual double CalcFunction(Excitation& excite, Function* f, bool derivative);

};

} // end of namespace




#endif /* LBMSIMP_HH_ */
