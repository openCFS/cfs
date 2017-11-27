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
class HeatPDE;
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
  typedef enum { NO_SYSTEM = -1, PIEZOCOUPLING, MECH, ELEC, HEAT, ACOUSTIC, LBM } System;

  /** calls the proper constructor */
  static OptimizationMaterial* CreateInstance(System sys, ErsatzMaterial* em, Context* ctxt);

  /** Here we store the system enum */
  static Enum<System> system;

  System GetSystem() const { return system_; }

  /** works fine for standard single pde SIMP stuff */
  virtual const DenseMatrix& Stiffness(const Elem* elem, bool bimaterial = false, int multimaterial = -1) {
    EXCEPTION("overload!");
  }

  virtual const DenseMatrix& Mass(const Elem* elem, bool bimaterial = false, int multimaterial = -1, DesignElement::Type direction = DesignElement::NO_DERIVATIVE) {
    EXCEPTION("overload!");
  }

  /** determines if we have a complex element matrix. This is the case for damped material or Bloch mode analysis with complex BOp*/
  bool ComplexElementMatrix(RegionIdType reg = NO_REGION_ID) const;

  virtual shared_ptr<CoefFunctionOpt> GetMatCoef(const std::string& integrator, BiLinFormContext* context = NULL, RegionIdType reg_id = NO_REGION_ID);

protected:


  /** @param em for the standard constructor to be used in ErsatzMaterial.
   * @param space  This allows to gain the MassMatrix for pamping where we might have no ErsatzMaterial  */
  OptimizationMaterial(ErsatzMaterial* em, Context* ctxt, DesignSpace* space = NULL);


  /** This is actually legacy code - the data shall always come from DesignSpace::localElementCace */
  template <class T>
  void GetElementMatrix(Matrix<T>& out, const std::string& integrator, const Elem* elem, bool lower_bimat = false,
                        DesignElement::Type direction = DesignElement::NO_DERIVATIVE, Global::ComplexPart entryType =  (Global::ComplexPart) 4711);


  /** retrieves the data from LocalElementCache */
  template <class T>
  const Matrix<T>& GetCachedElementMatrix(const std::string& integrator, const Elem* elem, bool lower_bimat = false, DesignElement::Type direction = DesignElement::NO_DERIVATIVE);

  const DenseMatrix& GetCachedElementMatrixDM(const std::string& integrator, const Elem* elem, bool lower_bimat = false, DesignElement::Type direction = DesignElement::NO_DERIVATIVE);

  /** Very similar to GetElementMatrix() but for the vector, e.g. for rhs linear forms */
  void GetElementVector(LinearForm* form, Vector<double>& out, const Elem* elem = NULL,
                          BaseMaterial* bimaterial = NULL, const Vector<double>* ts = NULL);


  /** service function handling multimaterial
   * @param de we need elem->region, multimaterial might be NULL otherwise we check for index */
  /*Matrix<double>& GeneralStiffness(std::map<RegionIdType, StdVector<Matrix<double> > >& map,
      const DesignElement* de, MaterialClass mc, DesignElement::Type direction, double factor, bool transposed); */

  DesignElement* FirstDesignByRegion(RegionIdType reg);

  StdVector<RegionIdType> regionIds;

  /** might be NULL when using other constructor */
  ErsatzMaterial* opt;

  DesignSpace* space;
  
  /** our context, we are assigned to context */
  Context* ctxt_;

  // short cuts
  bool transient_; // shall become ctxt->transient
  // harmonic, eigenvalue, transient
  bool needs_mass_;

  // what we are;
  System system_;
};

class MechMat : public OptimizationMaterial
{
public:
  
  MechMat(ErsatzMaterial* em, Context* ctxt);
  MechMat(DesignSpace* space);
  
  /** Get the ElementStiffness Matrix for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @param bimaterial if true gets the material from the design space by the element's region
   * @param multimaterial index or negative
   * @param direction if given, calculate derivative of Stiffness Matrix instead
   * @return a pointer to the Element Stiffness Matrix*/
  const DenseMatrix& MechStiffness(const Elem* elem, bool bimaterial = false, int multimaterial = -1, DesignElement::Type direction = DesignElement::NO_DERIVATIVE);

 /** overwrites OptimizationMaterial::Stiffness */
  const DenseMatrix& Stiffness(const Elem* elem, bool bimaterial = false, int multimaterial = -1 ) {
    return MechStiffness(elem, bimaterial, multimaterial, DesignElement::NO_DERIVATIVE);
  }

  /** Get the ElementMass Matrix for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @param direction if given, calculate derivative of mass Matrix instead
   * @return a pointer to the Element Mass Matrix*/
  const DenseMatrix& MechMass(const Elem* elem,  bool bimaterial = false, int multimaterial = -1, DesignElement::Type direction = DesignElement::NO_DERIVATIVE);

  /** overwrites OptimizationMaterial::Mass */
  const DenseMatrix& Mass(const Elem* elem, bool bimaterial = false, int multimaterial = -1, DesignElement::Type direction = DesignElement::NO_DERIVATIVE) {
    return MechMass(elem, bimaterial, multimaterial, direction);
  }

  /** The the rhs-contribution for full material for the current test strain. There is no caching!
   * @param testStrain optional value, otherwise the current set excitation set. You need it for homogenization! */
  const Vector<double>& MechStrainRHS(const Elem* elem, MechPDE::TestStrain testStrain = MechPDE::NOT_SET);

protected:
  
  /** We do not cache the vectors but always precalculate them */
  Vector<double> mechStrainRHS;

private:

  /** actual Constructor */
  void Init();
};


class AcouMat : public OptimizationMaterial
{
public:
  AcouMat(ErsatzMaterial* em, Context* ctxt);

  const Matrix<double>& AcouStiffness(const Elem* elem, bool bimaterial);

  /** overwrites OptimizationMaterial::Stiffness */
  const DenseMatrix& Stiffness(const Elem* elem, bool bimaterial = false) {
    return AcouStiffness(elem, bimaterial);
  }

  const Matrix<double>& AcouMass(const Elem* elem,  bool bimaterial);

  /** overwrites OptimizationMaterial::Mass */
  const Matrix<double>& Mass(const Elem* elem, bool bimaterial = false) {
    return AcouMass(elem, bimaterial);
  }

};


class PiezoElecMat : public MechMat
{
public:
  PiezoElecMat(ErsatzMaterial* em, Context* ctxt);
  
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
  
private:
  /** The elec stiffness matrix $K_{\phi \phi}$. */
  //std::map<RegionIdType, StdVector<Matrix<double> > > elecStiffness_map;

  /** The negative elec stiffness matrix $-K_{\phi \phi}$. */
  //std::map<RegionIdType, StdVector<Matrix<double> > > elecStiffness_neg_map;
  
  /** The coupling stiffness matrix $K_{u \phi}$ */
  //std::map<RegionIdType, StdVector<Matrix<double> > > coupledStiffness_map;

  /** The transposed coupling stiffness matrix $K_{u \phi}^T$ */
  //std::map<RegionIdType, StdVector<Matrix<double> > > coupledStiffnessTransposed_map;
};




/** For Jannis' Maxwell homogenization. The PiezoElecMat has elec for piezo */
class ElecMat : public OptimizationMaterial
{
public:

  ElecMat(ErsatzMaterial* em, Context* ctxt);

  /** Initialize elecStiffness_map again (needed for BITENSOR objective) */
  void ReInit();

  /** Get the ElementStiffness Matrix for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @param bimaterial if true gets the material from the design space by the element's region
   * @param direction if given, calculate derivative of Stiffness Matrix instead
   * @return a pointer to the Element Stiffness Matrix*/
  const Matrix<std::complex<double> >& ElecStiffness(const Elem* elem, bool bimaterial = false, const DesignElement::Type direction = DesignElement::NO_DERIVATIVE);


//  /** overwrites OptimizationMaterial::Stiffness */
//  Matrix<std::complex<double> >& Stiffness(const Elem* elem, bool bimaterial = false) {
//    return ElecStiffness(elem, bimaterial, DesignElement::NO_DERIVATIVE);
//  }
};


class HeatMat : public OptimizationMaterial
{
public:
  HeatMat(ErsatzMaterial* em, Context* ctxt);

  /** Get the ElementStiffness Matrix for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @param bimaterial if true gets the material from the design space by the element's region
   * @param multimaterial index or negative
   * @param direction if given, calculate derivative of Stiffness Matrix instead
   * @return a pointer to the Element Stiffness Matrix*/
  const DenseMatrix& HeatStiffness(const Elem* elem, bool bimaterial = false, int multimaterial = -1);

  /** overwrites OptimizationMaterial::Stiffness */
  const DenseMatrix& Stiffness(const Elem* elem, bool bimaterial = false, int multimaterial = -1 ) {
    return HeatStiffness(elem, bimaterial, multimaterial);
  }

  /** Get the ElementMass Matrix for this element, this is the region constant version
   * @param elem the Element for which the Matrix should be returned
   * @param direction if given, calculate derivative of mass Matrix instead
   * @return a pointer to the Element Mass Matrix*/
  const DenseMatrix& HeatMass(const Elem* elem,  bool bimaterial = false, int multimaterial = -1);

protected:
  SinglePDE* GetPDE();
};



class LBMMat : public OptimizationMaterial
{
public:
  LBMMat(ErsatzMaterial* em, Context* ctxt);
};


}

#endif /* OPTMATERIAL_HH_ */
