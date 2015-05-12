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
  virtual void SetElementK(DesignElement* de, const TransferFunction* tf, Application app, DenseMatrix* mat_out, bool derivative = true, CalcMode calcMode = STANDARD, double ev = -1.0);

  /** @see PiezoSIMO::SetElementK() */
   virtual void SetElementKMapping(DesignElement* de, BaseDesignElement::Type type, const TransferFunction* tf, Application app, DenseMatrix* mat_out, CalcMode calcMode, bool derivative = true);

};

} // end of namespace

#endif /* PIEZOPARAMMAT_HH_ */
