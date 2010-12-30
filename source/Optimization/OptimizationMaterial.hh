#ifndef OPTMATERIAL_HH_
#define OPTMATERIAL_HH_

#include "General/Enum.hh"
#include "General/environment.hh"
#include "Optimization/Design/DesignElement.hh"

namespace CoupledField 
{

class ErsatzMaterial;
class ElecPDE;
class MechPDE;
class HeatCondPDE;
class BaseForm;
class LinearForm;

/** For Optimization problems does this class provide an interface to the actual physics.
 * While ErsatzMaterial itself contains a vector of pdes and the solutions for these
 * pdes, does the OptimizationMaterial and it's concrete child classes contain the matrices. 
 * 
 * This is not much content but by this extraction we can optimize mechanics, acoustics,
 * piezoelectric materials, ... and furthermore do it by SIMP and LevelSet without code 
 * doubling. 
 */
class OptimizationMaterial
{
public:
  
  OptimizationMaterial(ErsatzMaterial* em);
  virtual ~OptimizationMaterial();
  
  /** Id of our material class */
  typedef enum { PIEZOCOUPLING, MECH, HEAT } System;

  /** Here we store the system enum */
  static Enum<System> system;

  System GetSystem() const { return system_; }
  
  /** <p>Get the original element matrix (stiffness, mass, ...)
   * which is constant for all isotropic elements.
   * This method is not only for mechanical SIMP but is also used by PiezoSIMP,
   * therefore it is generic.</p>
   * <p>If no element is given, the one from the first design element is used.</p>
   * <p>All transfer functions are disabled during this method. Call only for
   * enabled transfer functions (default)</p>
   * @param form to be extracted via GetForm()
   * @param out here the element stiffness matrix written. e.h. K_uu which is \int B E B
   * @param elem if not given the first design element is used, otherwise the provided one
   * @param factor in piezoelectricity K_pp is -1* BDBInt */
  void GetElementMatrix(BaseForm* form, Matrix<double>& out, const Elem* elem = NULL,
                        BaseMaterial* bimaterial = NULL,
                        const DesignElement::Type direction = DesignElement::NO_DERIVATIVE, double factor = 1.0);
  
  /** Very similar to GetElementMatrix() but for the vector, e.g. for rhs linear forms */
  void GetElementVector(LinearForm* form, Vector<double>& out, const Elem* elem = NULL,
                        BaseMaterial* bimaterial = NULL, const Vector<double>* ts = NULL);

protected:
  StdVector<RegionIdType> regionIds;

  ErsatzMaterial* opt;
  
  // what we are;
  System system_;

private:

  /** This is the common implementation for GetElementMatrix() and GetElementVector() */
  void GetElementEntity(BaseForm* form, Matrix<double>* mat_out, Vector<double>* vec_out, const Elem* elem = NULL,
                        BaseMaterial* bimaterial = NULL,
                        const DesignElement::Type direction = DesignElement::NO_DERIVATIVE, const Vector<double>* ts = NULL);
  
};

class OptMechMat : public OptimizationMaterial
{
public:
  
  OptMechMat(ErsatzMaterial* em);
  
  /** Get the ElementStiffness Matrix for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @param bimaterial if true gets the material from the design space by the element's region
   * @param direction if given, calculate derivative of Stiffness Matrix instead
   * @return a pointer to the Element Stiffness Matrix*/
  const Matrix<double>& MechStiffness(const Elem* elem, bool bimaterial = false, const DesignElement::Type direction = DesignElement::NO_DERIVATIVE);

  /** Get the ElementMass Matrix for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @param direction if given, calculate derivative of mass Matrix instead
   * @return a pointer to the Element Mass Matrix*/
  const Matrix<double>& MechMass(const Elem* elem,  bool bimaterial = false, const DesignElement::Type direction = DesignElement::NO_DERIVATIVE);
  
  /** The the rhs-contribution for full material for the current test strain. There is no caching!
   * @param testStrain optional value, otherwise the current set excitation set. You need it for homogenization! */
  const Vector<double>& MechStrainRHS(const Elem* elem, MechPDE::TestStrain testStrain = MechPDE::NOT_SET);

protected:  
  /** The mechanical element stiffness matrix is constant.
   * We store the results for standard (first) and bimaterial (second)  */
  std::map<RegionIdType, std::pair<Matrix<double>, Matrix<double> > > mechStiffness_map;

  /** The mechanical element mass matrix is also constant. Only for harmonic!
   * @see mechStiffness_map*/
  std::map<RegionIdType, std::pair<Matrix<double>, Matrix<double> > > mechMass_map;
  
  /** We do not cache the vectors but always precalculate them */
  Vector<double> mechStrainRHS;

  MechPDE* mech;
};

class OptPiezoMat : public OptMechMat
{
public:
  OptPiezoMat(ErsatzMaterial* em);
  
  /** Get the elec stiffness matrix $K_{\phi \phi}$ for this element, this is the region constant version.
   * Note, that in piezoelectriciy one usually requires -K_pp, use the the factor -1! 
   * @param elem the Element for which the Matrix should be returned
   * @param factor has to be either 1 or -1 
   * @return a pointer to the Element Mass Matrix*/
  const Matrix<double>& ElecStiffness(Elem* elem, int factor);

  /** Get the coupling stiffness matrix $K_{u \phi}$ for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @return a pointer to the Element Mass Matrix*/
  const Matrix<double>& CoupledStiffness(Elem* elem);

  /** Get the elec stiffness matrix $K_{\phi \phi}$ for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @return a pointer to the Element Mass Matrix*/
  const Matrix<double>& CoupledStiffnessTransposed(Elem* elem);
  
private:
  /** The elec stiffness matrix $K_{\phi \phi}$. */
  std::map<RegionIdType, Matrix<double> > elecStiffness_map;

  /** The negative elec stiffness matrix $-K_{\phi \phi}$. */
  std::map<RegionIdType, Matrix<double> > elecStiffness_neg_map;
  
  /** The coupling stiffness matrix $K_{u \phi}$ */
  std::map<RegionIdType, Matrix<double> > coupledStiffness_map;

  /** The transposed coupling stiffness matrix $K_{u \phi}^T$ */
  std::map<RegionIdType, Matrix<double> > coupledStiffnessTransposed_map;
  
  ElecPDE* elec;
};

class HeatMat : public OptimizationMaterial
{
public:
  HeatMat(ErsatzMaterial* em);
  
protected:
  HeatCondPDE* heat;
};
}


#endif /* OPTMATERIAL_HH_ */
