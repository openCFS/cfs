#ifndef PARAMMAT_HH_
#define PARAMMAT_HH_

#include "Optimization/SIMP.hh"

namespace CoupledField {

class MechMat;
class DenseMatrix;
class DesignElement;
class TransferFunction;

  class ParamMat : public SIMP
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
    virtual void SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* mat_out, bool derivative, CalcMode calcMode = STANDARD, double ev = -1.0) override
    {
      if(f->ctxt->IsComplex())
      {
        if(f->ctxt->mat->ComplexElementMatrix(de->elem->regionId)) // handles also bloch which real material but complex BOp
          SetElementK<Complex, Complex >(f, de, tf, app, mat_out, derivative, calcMode, ev);
        else
          SetElementK<Complex, double >(f, de, tf, app, mat_out, derivative, calcMode, ev);
      }
      else
        SetElementK<double,double>(f, de, tf, app, mat_out, derivative, calcMode, ev);
    }

    /** see the SIMP case
     * T1 is the out matrix type
     * T2 is the element matrix type */
    template <class T1, class T2>
    void SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* out, bool derivative = true, CalcMode mode = STANDARD, double ev = -1.0);

  };

}

#endif /*PARAMMAT_HH_*/
