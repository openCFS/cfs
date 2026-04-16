#ifndef PIEZOSIMP_HH_
#define PIEZOSIMP_HH_

#include <complex>

#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Function.hh"
#include "Optimization/Objective.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/SIMP.hh"
#include "Utils/ToolsFull.hh"

namespace CoupledField {
class DenseMatrix;
class DesignElement;
class Excitation;
class TransferFunction;
}  // namespace CoupledField

namespace CoupledField
{
class ElecPDE;
class PiezoElecMat;

/** Extension from lin elast SIMP to the piezoelectric case */
class PiezoSIMP : public SIMP
{
public:
  PiezoSIMP();

  ~PiezoSIMP();

  /** e.g. stuff that needs a this pointer of the optimization problem */
  void PostInit();

protected:

  /** overloads SIMP::CalcFunction()
   * @see ErsatzMaterial::CalcFunction */
  double CalcFunction(Excitation& excite, Function* f, bool derivative);

  virtual void ConstructAdjointRHS(Excitation& excite, Function* f)
  {
    if(f->GetType() == Objective::ELEC_ENERGY)
    {
      if(context->IsComplex()) ConstructAdjointRHS<std::complex<double> >(excite, f);
                         else ConstructAdjointRHS<double>(excite, f);
    }
    else
    {
      SIMP::ConstructAdjointRHS(excite, f); // EM
    }
  }
  
  
  /** This is our part for CalcU1KU2() -> This set the matrix derivatives form App::ELEC and
   * App::PIEZO_COUPLING ( + transposed) */
  virtual void SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* out, bool derivative = true, CalcMode calcMode = STANDARD, double ev = -1.0)
  {
    if(f->ctxt->IsComplex())
      SetElementK<std::complex<double>, double >( f, de, tf, app, out, derivative, calcMode, ev);
    else
      SetElementK<double, double>( f, de, tf, app, out, derivative, calcMode, ev);
  }
private:

  template <class T1, class T2>
  void SetElementK(Function* f, DesignElement* de, const TransferFunction* tf, App::Type app, DenseMatrix* out, bool derivative = true, CalcMode calcMode = STANDARD, double ev = -1.0);


  /** Calculate the electrix enegy p^T K_pp p resp. p^T K_pp p^* */
  template <class T>
  double CalcElecEnergy(Excitation& excite);

  /** Sets -K_pp p or -K_pp p^* */
  template <class T>
  void ConstructAdjointRHS(Excitation& excite, Function* cost);


  /** The electric rhs, real or complex */
  DesignDependentRHS elecRHS;
};

} // end of namespace


#endif /*PIEZOSIMP_HH_*/
