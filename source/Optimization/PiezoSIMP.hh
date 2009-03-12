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

protected:

  /** Calculates gradients in the form <l, Ku> */
  void CalcObjectiveGradient(Excitation& excite);

private:

  /** This is our part for CalcU1KU2() -> This set the matrix derivatives form ELEC and
   * PIEZO_COUPLING ( + transposed) */
  virtual void SetElementK(DesignElement* de, Application app, DenseMatrix* out)
  {
    if(harmonic) SetElementK<std::complex<double> >(de, app, out);
            else SetElementK<double>(de, app, out);
  }

  template <class T>
  void SetElementK(DesignElement* de, Application app, DenseMatrix* out);


  /** We have for the direct coupled piezo the corresponding ElecPDE, is also in SIMP::pde */
  ElecPDE* elec;


  /** The elec stiffness matrix $K_{\phi \phi}$. $K_{uu}$ is already as mechStiffness in SIMP.hh */
  std::map<RegionIdType, Matrix<double> > elecStiffness_map;

  /** The coupling stiffness matrix $K_{u \phi}$ */
  std::map<RegionIdType, Matrix<double> > coupledStiffness_map;

  /** The transposed coupling stiffness matrix $K_{u \phi}^T$ */
  std::map<RegionIdType, Matrix<double> > coupledStiffnessTransposed_map;

  /** The electric rhs, real or complex */
  SurfaceRef elecRHS;

  /** Get the elec stiffness matrix $K_{\phi \phi}$ for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @return a pointer to the Element Mass Matrix*/
  const Matrix<double>& ElecStiffness(Elem* elem);

  /** Get the coupling stiffness matrix $K_{u \phi}$ for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @return a pointer to the Element Mass Matrix*/
  const Matrix<double>& CoupledStiffness(Elem* elem);

  /** Get the elec stiffness matrix $K_{\phi \phi}$ for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @return a pointer to the Element Mass Matrix*/
  const Matrix<double>& CoupledStiffnessTransposed(Elem* elem);

  /** This are logging variables for LogFileLine */
  double log_coupled_;
  double log_coupled_simp_;
  double log_elec_;
  double log_elec_simp_;

};

} // end of namespace


#endif /*PIEZOSIMP_HH_*/
