#ifndef PIEZOPARAMMAT_HH_
#define PIEZOPARAMMAT_HH_

/** The clou with PiezoParamMat ist, that it derives from PiezoSIMP.hh as there is not much to be used from the
 * lightweight ParamMat class.
 * Used to do Piezo FMO */

#include "Optimization/PiezoSIMP.hh"

namespace CoupledField
{

class PiezoParamMat : public PiezoSIMP
{
public:
  PiezoParamMat();

  ~PiezoParamMat();

private:

  /** @see PiezoSIMO::SetElementK() */
  virtual void SetElementK(DesignElement* de, const TransferFunction* tf, Application app, DenseMatrix* out, CalcMode calcMode, bool derivative = true);
};

} // end of namespace

#endif /* PIEZOPARAMMAT_HH_ */
