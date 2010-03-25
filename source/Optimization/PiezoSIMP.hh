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
  virtual double CalcObjective(Excitation& excite, Objective* cost)
  {
    if(cost->GetType() == Objective::ELEC_ENERGY)
    {
      return harmonic ? CalcElecEnergy<std::complex<double> >(excite) : CalcElecEnergy<double>(excite); 
    }
    else
      return SIMP::CalcObjective(excite, cost);
  }
  
  virtual void ConstructAdjointRHS(Excitation& excite, Objective* cost)
  {
    if(cost->GetType() == Objective::ELEC_ENERGY)
    {
      if(harmonic) ConstructAdjointRHS<std::complex<double> >(excite, cost);
              else ConstructAdjointRHS<double>(excite, cost);
    }
    else
    {
      SIMP::ConstructAdjointRHS(excite, cost); // EM
    }
  }
  
  /** Calculates gradients in the form <l, Ku> */
  void CalcObjectiveGradient(Excitation& excite, Objective* cost);


private:

  /** Calculate the electrix enegy p^T K_pp p resp. p^T K_pp p^* */  
  template <class T>
  double CalcElecEnergy(Excitation& excite);

  /** Sets -K_pp p or -K_pp p^* */
  template <class T>
  void ConstructAdjointRHS(Excitation& excite, Objective* cost);
  
  
  /** This is our part for CalcU1KU2() -> This set the matrix derivatives form ELEC and
   * PIEZO_COUPLING ( + transposed) */
  virtual void SetElementK(DesignElement* de, Application app, DenseMatrix* out, CalcMode calcMode, bool derivative = true)
  {
    if(harmonic) SetElementK<std::complex<double> >(de, app, out, calcMode, derivative);
            else SetElementK<double>(de, app, out, calcMode, derivative);
  }

  template <class T>
  void SetElementK(DesignElement* de, Application app, DenseMatrix* out, CalcMode calcMode, bool derivative = true);

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
