#ifndef PIEZOSIMP_HH_
#define PIEZOSIMP_HH_

#include "Optimization/SIMP.hh"

namespace CoupledField
{
class ElecPDE;

/** Extension from lin elast SIMP to the piezoelectric case */
class PiezoSIMP : public SIMP
{
public:
  PiezoSIMP();

  ~PiezoSIMP();

  /** e.g. stuff that needs a this pointer of the optimization problem */
  void PostInit();
  
  
  /** Evaluates the gradient of the cost function. Saves the result to data.objective_gradient  
   * @param grad_out if not NULL copy write the (design space size) results. 
   * If filtering is enabled, this is automatically filtered. */
  void CalcObjectiveGradient(double* grad_out);

private:

  /** This is our part for CalcU1KU2() -> This set the matrix derivatives form ELEC and
   * PIEZO_COUPLING ( + transposed) */
  virtual void SetElementK(DesignElement* de, Application app, CFSMatrix* out)
  {
    if(harmonic) SetElementK<std::complex<double> >(de, app, out);
            else SetElementK<double>(de, app, out);  
  }
  
  template <class T>
  void SetElementK(DesignElement* de, Application app, CFSMatrix* out);
  

  /** We have for the direct coupled piezo the corresponding ElecPDE, is also in SIMP::pde */
  ElecPDE* elec;
  

  /** The elec stiffness matrix $K_{\phi \phi}$. $K_{uu}$ is already as mechStiffness in SIMP.hh */
  Matrix<double> elecStiffness;

  /** The coupling stiffness matrix $K_{u \phi}$ */
  Matrix<double> coupledStiffness;

  /** The transposed coupling stiffness matrix $K_{u \phi}^T$ */         
  Matrix<double> coupledStiffnessTransposed;

  /** The electric rhs, real or complex */
  SurfaceRef elecRHS;

  /** This are logging variables for LogFileLine */
  double log_coupled_;
  double log_coupled_simp_;
  double log_elec_;
  double log_elec_simp_;

};

} // end of namespace


#endif /*PIEZOSIMP_HH_*/
