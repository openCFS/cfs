#ifndef OPTMATERIAL_HH_
#define OPTMATERIAL_HH_

#include <stddef.h>
#include <complex>
#include <map>
#include <utility>

#include "General/Enum.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Optimization/Design/DesignElement.hh"
#include "PDE/MechPDE.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {
class BaseMaterial;
struct Elem;
}  // namespace CoupledField

namespace CoupledField 
{

class AcousticPDE;
class BiLinForm;
class BiLinFormContext;
class DesignSpace;
class ElecPDE;
class ErsatzMaterial;
class HeatCondPDE;
class LinearForm;
class LatticeBoltzmannPDE;
class CoefFunctionOpt;

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
  
  virtual ~OptimizationMaterial();
  
  /** Id of our material class */
  typedef enum { PIEZOCOUPLING, MECH, ELEC, HEAT, ACOUSTIC, LBM } System;

  /** calls the proper constructor */
  static OptimizationMaterial* CreateInstance(System sys, ErsatzMaterial* em);

  /** Here we store the system enum */
  static Enum<System> system;

  System GetSystem() const { return system_; }

  /** works fine for standard single pde SIMP stuff */
  virtual DenseMatrix& Stiffness(const Elem* elem, bool bimaterial = false, int multimaterial = -1) {
    EXCEPTION("overload!");
  }

  virtual const Matrix<double>& Mass(const Elem* elem, bool bimaterial = false, int multimaterial = -1) {
    EXCEPTION("overload!");
  }

protected:


  /** @param em for the standard constructor to be used in ErsatzMaterial.
   * @param space  This allows to gain the MassMatrix for pamping where we might have no ErsatzMaterial  */
  OptimizationMaterial(ErsatzMaterial* em, DesignSpace* space = NULL);

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
  void GetElementMatrix(BiLinearForm* form, Matrix<double>& out, const Elem* elem = NULL,
                        BaseMaterial* bimaterial = NULL,
                        const DesignElement::Type direction = DesignElement::NO_DERIVATIVE, double factor = 1.0);

  void GetElementMatrix(Matrix<double>& out, const std::string& integrator, const Elem* elem = NULL,
                        BaseMaterial* bimaterial = NULL, DesignElement::Type direction = DesignElement::NO_DERIVATIVE,  Global::ComplexPart entryType =  (Global::ComplexPart) 4711);

  /** Very similar to GetElementMatrix() but for the vector, e.g. for rhs linear forms */
  void GetElementVector(LinearForm* form, Vector<double>& out, const Elem* elem = NULL,
                        BaseMaterial* bimaterial = NULL, const Vector<double>* ts = NULL);

  /** Helper which extracts the FormContext from assemble using the optimization region
          * @param regionId the corresponding region
          * @param pde1 the first pde (e.g. mech)
          * @param pde2 this is either the same as pde1 or the coupling partner
          * @param integrator there is no nice enum yet :( e.g. linElastInt, MechInt, ... */
  BiLinearForm* GetForm(RegionIdType regionId, const std::string& integrator, Global::ComplexPart entryType = (Global::ComplexPart) 4711);


  /** service function handling multimaterial
   * @param de we need elem->region, multimaterial might be NULL otherwise we check for index */
  Matrix<double>& GeneralStiffness(std::map<RegionIdType, StdVector<Matrix<double> > >& map,
      const DesignElement* de, MaterialClass mc, StdPDE* pde1, StdPDE* pde2, DesignElement::Type direction,
      double factor, bool transposed);

  virtual shared_ptr<CoefFunctionOpt> GetMatCoef(const std::string& integrator, BiLinFormContext* context = NULL, RegionIdType reg_id = NO_REGION_ID);

  virtual SinglePDE* GetPDE() = 0;

  StdVector<RegionIdType> regionIds;

  /** might be NULL when using other constructor */
  ErsatzMaterial* opt;

  DesignSpace* space;
  
  // short cuts
  bool harmonic_;
  bool transient_;
  bool structured_;

  // what we are;
  System system_;



private:

  /** This is the common implementation for GetElementMatrix() and GetElementVector() */
  void GetElementEntity(BiLinearForm* form, Matrix<double>* mat_out, Vector<double>* vec_out, const Elem* elem = NULL,
                        BaseMaterial* bimaterial = NULL,
                        DesignElement::Type direction = DesignElement::NO_DERIVATIVE, const Vector<double>* ts = NULL);

};

class MechMat : public OptimizationMaterial
{
public:
  
  MechMat(ErsatzMaterial* em);
  MechMat(DesignSpace* space);
  
  /** Get the ElementStiffness Matrix for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @param bimaterial if true gets the material from the design space by the element's region
   * @param multimaterial index or negative
   * @param direction if given, calculate derivative of Stiffness Matrix instead
   * @return a pointer to the Element Stiffness Matrix*/
  DenseMatrix& MechStiffness(const Elem* elem, bool bimaterial = false, int multimaterial = -1, DesignElement::Type direction = DesignElement::NO_DERIVATIVE);

  /** overwrites OptimizationMaterial::Stiffness */
  DenseMatrix& Stiffness(const Elem* elem, bool bimaterial = false, int multimaterial = -1 ) {
    return MechStiffness(elem, bimaterial, multimaterial, DesignElement::NO_DERIVATIVE);
  }

  /** Get the ElementMass Matrix for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @param direction if given, calculate derivative of mass Matrix instead
   * @return a pointer to the Element Mass Matrix*/
  const Matrix<double>& MechMass(const Elem* elem,  bool bimaterial = false, int multimaterial = -1, DesignElement::Type direction = DesignElement::NO_DERIVATIVE);

  /** overwrites OptimizationMaterial::Mass */
  const Matrix<double>& Mass(const Elem* elem, bool bimaterial = false, int multimaterial = -1) {
    return MechMass(elem, bimaterial, multimaterial, DesignElement::NO_DERIVATIVE);
  }

  
  /** The the rhs-contribution for full material for the current test strain. There is no caching!
   * @param testStrain optional value, otherwise the current set excitation set. You need it for homogenization! */
  const Vector<double>& MechStrainRHS(const Elem* elem, MechPDE::TestStrain testStrain = MechPDE::NOT_SET);


protected:

  SinglePDE* GetPDE();

  /** The mechanical element stiffness matrix is constant.
   * We store multimaterial as a vector. No material one entry and legacy bimaterial two entries */
  std::map<RegionIdType, StdVector<Matrix<double> > > mechStiffness_map;
  std::map<RegionIdType, StdVector<Matrix<Complex> > > mechStiffness_mapC;

  /** The mechanical element mass matrix is also constant. Only for harmonic!
   * @see mechStiffness_map*/
  std::map<RegionIdType, StdVector<Matrix<double> > > mechMass_map;
  
  /** We do not cache the vectors but always precalculate them */
  Vector<double> mechStrainRHS;

  MechPDE* mech;

private:

  /** actual Constructor */
  void Init();
};


class AcouMat : public OptimizationMaterial
{
public:
  AcouMat(ErsatzMaterial* em);

  Matrix<double>& AcouStiffness(const Elem* elem, bool bimaterial);

  /** overwrites OptimizationMaterial::Stiffness */
  DenseMatrix& Stiffness(const Elem* elem, bool bimaterial = false) {
    return AcouStiffness(elem, bimaterial);
  }

  const Matrix<double>& AcouMass(const Elem* elem,  bool bimaterial);

  /** overwrites OptimizationMaterial::Mass */
  const Matrix<double>& Mass(const Elem* elem, bool bimaterial = false) {
    return AcouMass(elem, bimaterial);
  }

protected:

  SinglePDE* GetPDE();

  std::map<RegionIdType, std::pair<Matrix<double>, Matrix<double> > > acouStiffness_map;
  std::map<RegionIdType, std::pair<Matrix<double>, Matrix<double> > > acouMass_map;

  AcousticPDE* acou;
};


class PiezoElecMat : public MechMat
{
public:
  PiezoElecMat(ErsatzMaterial* em);
  
  /** Get the elec stiffness matrix $K_{\phi \phi}$ for this element, this is the region constant version.
   * Note, that in piezoelectriciy one usually requires -K_pp, use the the factor -1! 
   * @param de only elem for number and region is used and multimaterial is checked
   * @return a pointer to the Element Mass Matrix*/
  const Matrix<double>& ElecStiffnessPos(const DesignElement* de, DesignElement::Type direction = DesignElement::NO_DERIVATIVE);

  const Matrix<double>& ElecStiffnessNeg(const DesignElement* de, DesignElement::Type direction = DesignElement::NO_DERIVATIVE);

  /** Get the coupling stiffness matrix $K_{u \phi}$ for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @return a pointer to the Element Mass Matrix*/
  const Matrix<double>& CoupledStiffness(const DesignElement* de, DesignElement::Type direction = DesignElement::NO_DERIVATIVE);

  /** Get the elec stiffness matrix $K_{\phi \phi}$ for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @return a pointer to the Element Mass Matrix*/
  const Matrix<double>& CoupledStiffnessTransposed(const DesignElement* de, DesignElement::Type direction = DesignElement::NO_DERIVATIVE);
  
protected:
  SinglePDE* GetPDE();

private:
  /** The elec stiffness matrix $K_{\phi \phi}$. */
  std::map<RegionIdType, StdVector<Matrix<double> > > elecStiffness_map;

  /** The negative elec stiffness matrix $-K_{\phi \phi}$. */
  std::map<RegionIdType, StdVector<Matrix<double> > > elecStiffness_neg_map;
  
  /** The coupling stiffness matrix $K_{u \phi}$ */
  std::map<RegionIdType, StdVector<Matrix<double> > > coupledStiffness_map;

  /** The transposed coupling stiffness matrix $K_{u \phi}^T$ */
  std::map<RegionIdType, StdVector<Matrix<double> > > coupledStiffnessTransposed_map;
  
  ElecPDE* elec;
};




/** For Jannis' Maxwell homogenization. The PiezoElecMat has elec for piezo */
class ElecMat : public OptimizationMaterial
{
public:

  ElecMat(ErsatzMaterial* em);

  /** Initialize elecStiffness_map again (needed for BITENSOR objective) */
  void ReInit();

  /** Get the ElementStiffness Matrix for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @param bimaterial if true gets the material from the design space by the element's region
   * @param direction if given, calculate derivative of Stiffness Matrix instead
   * @return a pointer to the Element Stiffness Matrix*/
  Matrix<std::complex<double> >& ElecStiffness(const Elem* elem, bool bimaterial = false, const DesignElement::Type direction = DesignElement::NO_DERIVATIVE);


//  /** overwrites OptimizationMaterial::Stiffness */
//  Matrix<std::complex<double> >& Stiffness(const Elem* elem, bool bimaterial = false) {
//    return ElecStiffness(elem, bimaterial, DesignElement::NO_DERIVATIVE);
//  }

protected:

  SinglePDE* GetPDE();

  /** The electrostatic element stiffness matrix is constant.
   * We store the results for standard (first) and bimaterial (second)  */
  std::map<RegionIdType, std::pair<Matrix<std::complex <double> >, Matrix<std::complex <double> > > > elecStiffness_map;

  ElecPDE* elec;
};


class HeatMat : public OptimizationMaterial
{
public:
  HeatMat(ErsatzMaterial* em);

protected:
  SinglePDE* GetPDE();

  HeatCondPDE* heat;
};



class LBMMat : public OptimizationMaterial
{
public:
  LBMMat(ErsatzMaterial* em);

protected:

  SinglePDE* GetPDE() { assert(false); return NULL; } // FIXME

  LatticeBoltzmannPDE* lbm;
};


}

#endif /* OPTMATERIAL_HH_ */
