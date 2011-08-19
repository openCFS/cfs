#ifndef PARAMMAT_HH_
#define PARAMMAT_HH_

#include "Optimization/ErsatzMaterial.hh"

namespace CoupledField {

class MechMat;

  class ParamMat : public ErsatzMaterial
  {
  public:
    /** constructor for parametric material optimization, most is already done in ErsatzMaterial */
    ParamMat();

    /** virtual second phase initializer */
    virtual void PostInit();
    
  protected:
    
    /** provides the derivative of the material matrix/tensor at the current element, 
     * is called by CalcU1KU2 in ErsatzMaterial
     * @param de the current DesignElement (this provides the element as well as the direction)
     * @param app is ignored
     * @param outn pointer where the matrix should be stored */
    virtual void SetElementK(DesignElement* de, const TransferFunction* tf, Application app, DenseMatrix* mat_out, CalcMode calcMode, bool derivative);
    
    /** this is a shortcut to the material class */
    MechMat* mech_mat_;
  };

}

#endif /*PARAMMAT_HH_*/
