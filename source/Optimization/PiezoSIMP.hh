#ifndef PIEZOSIMP_HH_
#define PIEZOSIMP_HH_

#include <complex>

#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Function.hh"
#include "Optimization/Objective.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/SIMP.hh"
#include "Utils/tools.hh"

namespace CoupledField {
class DenseMatrix;
class DesignElement;
class Excitation;
class TransferFunction;
}  // namespace CoupledField

namespace CoupledField
{
class ElecPDE;
class PiezoelecMat;

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
      if(harmonic) ConstructAdjointRHS<std::complex<double> >(excite, f);
              else ConstructAdjointRHS<double>(excite, f);
    }
    else
    {
      SIMP::ConstructAdjointRHS(excite, f); // EM
    }
  }
  
private:

  /** Calculate the electrix enegy p^T K_pp p resp. p^T K_pp p^* */  
  template <class T>
  double CalcElecEnergy(Excitation& excite);

  /** Sets -K_pp p or -K_pp p^* */
  template <class T>
  void ConstructAdjointRHS(Excitation& excite, Function* cost);
  
  
  /** This is our part for CalcU1KU2() -> This set the matrix derivatives form ELEC and
   * PIEZO_COUPLING ( + transposed) */
  virtual void SetElementK(DesignElement* de, const TransferFunction* tf, Application app, DenseMatrix* out, CalcMode calcMode, bool derivative = true)
  {
    if(harmonic) SetElementK<std::complex<double> >(de, tf, app, out, calcMode, derivative);
            else SetElementK<double>(de, tf, app, out, calcMode, derivative);
  }

  template <class T>
  void SetElementK(DesignElement* de, const TransferFunction* tf, Application app, DenseMatrix* out, CalcMode calcMode, bool derivative = true);

  /** The electric rhs, real or complex */
  DesignDependentRHS elecRHS;
  
  /** shortcut to our pde, is also in ErsatzMaterial::pdes */
  ElecPDE* elec;
  
  /** This are logging variables for LogFileLine */
  double log_coupled_;
  double log_coupled_simp_;
  double log_elec_;
  double log_elec_simp_;

  /** is a cast of the ErsatzMaterial::material attrbiute. Set in PostInit() */
  PiezoelecMat* piezo_mat_;
  
};

} // end of namespace


#endif /*PIEZOSIMP_HH_*/
