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
#include "PDE/HeatPDE.hh"
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

  /** this structure describes a material form. Most systems have a stiffness and a mass variant. */
  struct FormID
  {
    /** form name */
    std::string   integrator;
    /** not same for stiff and mass but to have it all together */
    MaterialClass mc = NO_CLASS;
    /** specific material type */
    MaterialType  mt = NO_MATERIAL;
    /** here store the matrix we return by Stiffness() or Mass when we need to compute it and cannot return a cached one.
     * ThreadLocalStorage to allow parallelization */
    CfsTLS<Matrix<double> >  calc_real;
    CfsTLS<Matrix<Complex> > calc_cplx;
  };
  
  /** Id of our material class */
  typedef enum { NO_SYSTEM = -1, PIEZOCOUPLING, MECH, ELEC, HEAT, MAG, ACOUSTIC, LBM } System;

  /** calls the proper constructor */
  static OptimizationMaterial* CreateInstance(System sys, ErsatzMaterial* em, Context* ctxt);

  /** Here we store the system enum */
  static Enum<System> system;

  System GetSystem() const { return system_; }

  /** works fine for standard single pde SIMP stuff. We have a complex stiffness matrix for bloch.
   * The direction (param mat derivatives are only for MechMat::MechStiffness() */
  const DenseMatrix& Stiffness(const Elem* elem, bool bimaterial = false, int multimaterial = -1, DesignElement::Type mat_deriv = DesignElement::NO_DERIVATIVE) {
    return GetElementMatrix(stiff, elem, bimaterial, multimaterial, mat_deriv);
  }

  /** we have a dense mass matrix for bloch */
  const DenseMatrix& Mass(const Elem* elem, bool bimaterial = false, int multimaterial = -1, DesignElement::Type mat_deriv = DesignElement::NO_DERIVATIVE) {
    return GetElementMatrix(mass, elem, bimaterial, multimaterial, mat_deriv);
  }

  const DenseMatrix& GeometricStiffness(const Elem* elem, bool bimaterial = false, int multimaterial = -1, DesignElement::Type mat_deriv = DesignElement::NO_DERIVATIVE) {
    return GetElementMatrix(gstiff, elem, bimaterial, multimaterial, mat_deriv);
  }

  /** service function which indeed computes the matrix.
   * The underlaying CoefFunctionOpt needs to be org or opt before. State is resetted!
   * @param use stiff/mass from your OptimizationMaterial
   * @param this material is used for the BDBInt */
  template <class T>
  const Matrix<T>& ComputeElementMatrix(Matrix<T>& out, const FormID& form_id, const Elem* elem, shared_ptr<CoefFunction> shadow);

  /** If we don't have cached element matrices we compute them here. This is necessary for many param mat derivatives and also fallback.
   * @return this is out, just a convenience function. */
  template <class T>
  const Matrix<T>& ComputeElementMatrix(Matrix<T>& out, const std::string& integrator, const Elem* elem, bool lower_bimat = false,
                        DesignElement::Type direction = DesignElement::NO_DERIVATIVE, Global::ComplexPart entryType =  (Global::ComplexPart) 4711);

  /** determines if we have a complex element matrix. This is the case for damped material or Bloch mode analysis with complex BOp*/
  bool ComplexElementMatrix(RegionIdType reg = NO_REGION_ID) const;


  /** in contrast to CoefFunctionOpt* GetMatCoef() this is always the original material. Works for optimization domains and shall be extended for non-optimization domains!
   * @return don't store the coef!! */
  CoefFunction* GetOrgMatCoef(BaseBDBInt* bdb);

  CoefFunctionOpt* GetMatCoef(BiLinFormContext* context);

  CoefFunctionOpt* GetMatCoef(LinearFormContext* context);

  CoefFunctionOpt* GetMatCoef(const string& integrator, RegionIdType reg_id);

  CoefFunctionOpt* GetMatCoef(const FormID& form_id, RegionIdType reg_id) {
    return GetMatCoef(form_id.integrator, reg_id);
  }


  FormID stiff;
  FormID mass;
  FormID gstiff;

protected:

  /** @param em for the standard constructor to be used in ErsatzMaterial.
   * @param space  This allows to gain the MassMatrix for pamping where we might have no ErsatzMaterial  */
  OptimizationMaterial(ErsatzMaterial* em, Context* ctxt, DesignSpace* space = NULL);


  /** This is the common function for Stiffness() and Mass(). It either calls ComputeElementMatrix() or gets the stuff from LocalElementCache */
  const DenseMatrix& GetElementMatrix(FormID& id, const Elem* elem, bool bimaterial = false, int multimaterial = -1, DesignElement::Type mat_deriv = DesignElement::NO_DERIVATIVE);

  /** Very similar to GetElementMatrix() but for the vector, e.g. for rhs linear forms */
  void GetElementVector(LinearForm* form, Vector<double>& out, const Elem* elem = NULL,
                          BaseMaterial* bimaterial = NULL, const Vector<double>* ts = NULL);


  void GetOrgElementVector(LinearFormContext* lc, const Elem* elem, Vector<double>& elemVec);

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

  /** needs to be replaced by derived classes */
  System system_ = NO_SYSTEM;

};

class MechMat : public OptimizationMaterial
{
public:
  
  MechMat(ErsatzMaterial* em, Context* ctxt);
  MechMat(DesignSpace* space);
  

  /** specific version for piezo */
  const Matrix<double>& MechStiffness(const DesignElement* de, bool lower_bimat = false, int multi_material = -1, DesignElement::Type direction = DesignElement::NO_DERIVATIVE) {
     return dynamic_cast<const Matrix<double>&>(Stiffness(de->elem, lower_bimat, multi_material, direction));
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

  // for numerical homogenization:
  // we compute chi_e^0 by solving element-wise system K_e chi_e^0 = f_e
  // with respect to scalar temperature the dimensions are: K-e (4 x 4), chi_e^0 (4 x 1), f_e (4 x 1)
  // if mesh is regular, we only need to do this once
  // Ke is the element stiffness matrix
  Vector<double> CalcElementTemperature(Context* ctxt, const Elem* elem, HeatPDE::TestStrain testStrain);

private:
  // temperature on each node of finite element
  // assume a regular FE mesh
  // for each region and for each excitation, store chi_e^0
  std::map<RegionIdType, StdVector<Vector<double> > > elem_temp;
};


class MagMat : public OptimizationMaterial
{
public:
  MagMat(ErsatzMaterial* em, Context* ctxt);

  /** To get the rhs in magTopOpt, it is necessary to switch to the ORG Material.
   * @return the elem rhs contribution for the original case, not design dependent. */
  const Vector<double>& MagExcitationRHS(const std::string& integrator, const Elem* elem);

protected:

  /** We do not cache the vectors but always precalculate them. NOT thread safe */
  Vector<double> mag_excitation;
};



class LBMMat : public OptimizationMaterial
{
public:
  LBMMat(ErsatzMaterial* em, Context* ctxt);
};


}

#endif /* OPTMATERIAL_HH_ */
