#ifndef PARAMMAT_HH_
#define PARAMMAT_HH_

#include "Optimization/ErsatzMaterial.hh"

namespace CoupledField {

  class ParamMat : public ErsatzMaterial
  {
  public:
    /** constructor for parametric material optimization, most is already done in ErsatzMaterial */
    ParamMat();
    
  protected:
    
    /** provides the derivative of the material matrix/tensor at the current element, 
     * is called by CalcU1KU2 in ErsatzMaterial
     * @param de the current DesignElement (this provides the element as well as the direction)
     * @param app is ignored
     * @param outn pointer where the matrix should be stored */
    virtual void SetElementK(DesignElement* de, Application app, DenseMatrix* mat_out);
    
  };

}

#endif /*PARAMMAT_HH_*/
