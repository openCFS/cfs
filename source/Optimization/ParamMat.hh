#ifndef PARAMMAT_HH_
#define PARAMMAT_HH_

#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Optimization.hh"
#include "Domain/elem.hh"

namespace CoupledField {

class MechMat;
class DenseMatrix;
class DesignElement;
class TransferFunction;

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
    virtual void SetElementK(DesignElement* de, const TransferFunction* tf, Application app, DenseMatrix* out, CalcMode calcMode, bool derivative = true)
    {
      if(harmonic){
        if(pde->HasComplexMatData(de->elem->regionId))
          SetElementK<Complex, Complex >(de, tf, app, out, calcMode, derivative);
        else
          SetElementK<Complex, double >(de, tf, app, out, calcMode, derivative);
      }
      else SetElementK<double,double>(de, tf, app, out, calcMode, derivative);
    }
    
    /** this is a shortcut to the material class */
    MechMat* mech_mat_;

  private:


    template <class T1, class T2>
    void SetElementK(DesignElement* de, const TransferFunction* tf, Application app, DenseMatrix* out, CalcMode calcMode, bool derivative = true);


  };

}

#endif /*PARAMMAT_HH_*/
