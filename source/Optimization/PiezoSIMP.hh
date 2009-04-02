#ifndef PIEZOSIMP_HH_
#define PIEZOSIMP_HH_

#include "Optimization/SIMP.hh"

namespace CoupledField
{
class ElecPDE;
class OptPiezoMat;

/** Extension from lin elast SIMP to the piezoelectric case */
class PiezoSIMP : public SIMP
{
public:
  PiezoSIMP();

  ~PiezoSIMP();

  /** e.g. stuff that needs a this pointer of the optimization problem */
  void PostInit();

protected:

  /** overwrite this method for own objectives. */
  virtual double CalcObjective(Excitation& excite)
  {
    if(cost->type == ELEC_ENERGY)
    {
      return harmonic ? CalcElecEnergy<std::complex<double> >(excite) : CalcElecEnergy<double>(excite); 
    }
    else
      return SIMP::CalcObjective(excite);
  }
  
  virtual void ConstructAdjointRHS(Excitation& excite)
  {
    if(cost->type == ELEC_ENERGY)
    {
      if(harmonic) ConstructAdjointRHS<std::complex<double> >(excite);
              else ConstructAdjointRHS<double>(excite);
    }
    // else
    SIMP::ConstructAdjointRHS(excite);
  }
  
  /** Calculates gradients in the form <l, Ku> */
  void CalcObjectiveGradient(Excitation& excite);

private:

  /** Calculate the electrix enegy p^T K_pp p resp. p^T K_pp p^* */  
  template <class T>
  double CalcElecEnergy(Excitation& excite);

  /** Sets -K_pp p or -K_pp p^* */
  template <class T>
  void ConstructAdjointRHS(Excitation& excite);  
  
  
  /** This is our part for CalcU1KU2() -> This set the matrix derivatives form ELEC and
   * PIEZO_COUPLING ( + transposed) */
  virtual void SetElementK(DesignElement* de, Application app, DenseMatrix* out, CalcMode calcMode)
  {
    if(harmonic) SetElementK<std::complex<double> >(de, app, out, calcMode);
            else SetElementK<double>(de, app, out, calcMode);
  }

  template <class T>
  void SetElementK(DesignElement* de, Application app, DenseMatrix* out, CalcMode calcMode);

  /** The electric rhs, real or complex */
  SurfaceRef elecRHS;
  
  /** shortcut to our pde, is also in ErsatzMaterial::pdes */
  ElecPDE* elec;
  
  /** This are logging variables for LogFileLine */
  double log_coupled_;
  double log_coupled_simp_;
  double log_elec_;
  double log_elec_simp_;

  /** is a cast of the ErsatzMaterial::material attrbiute. Set in PostInit() */
  OptPiezoMat* piezo_mat_;
  
};

} // end of namespace


#endif /*PIEZOSIMP_HH_*/
