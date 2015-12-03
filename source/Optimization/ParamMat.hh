#ifndef PARAMMAT_HH_
#define PARAMMAT_HH_

#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Optimization.hh"

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
    virtual void SetElementK(DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* mat_out, bool derivative, CalcMode calcMode = STANDARD, double ev = -1.0)
    {
      if(context->IsComplex())
      {
        if(context->mat->ComplexElementMatrix(de->elem->regionId)) // handles also bloch which real material but complex BOp
          SetElementK<Complex, Complex >(de, tf, app, mat_out, derivative, calcMode, ev);
        else
          SetElementK<Complex, double >(de, tf, app, mat_out, derivative, calcMode, ev);
      }
      else
        SetElementK<double,double>(de, tf, app, mat_out, derivative, calcMode, ev);
    }
    
    virtual void SetElementKMapping(DesignElement* de, BaseDesignElement::Type type, const TransferFunction* tf, App::Type app, DenseMatrix* mat_out, CalcMode calcMode, bool derivative);

    /** see the SIMP case
     * T1 is the out matrix type
     * T2 is the element matrix type */
    template <class T1, class T2>
    void SetElementK(DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* out, bool derivative = true, CalcMode mode = STANDARD, double ev = -1.0);

  };

}

#endif /*PARAMMAT_HH_*/
