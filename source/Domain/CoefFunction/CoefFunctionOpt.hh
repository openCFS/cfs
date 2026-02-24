#ifndef COEFFUNCTIONOPT_HH_
#define COEFFUNCTIONOPT_HH_

#include <memory>

#include "CoefFunction.hh"
#include "Optimization/Design/DesignElement.hh"

namespace CoupledField {

class DesignSpace;
class BiLinearForm;
class LinearForm;

/** Replaces a material data CoefFunction with an ersatz material optimization version.
 * In the SIMP case the original const material is scaled, in the bi-material it is interpolated with an additional
 * material, in multi-material it is the weighted sum of various materials and in parametric optimization the material (tensor)
 * is completely constructed out of optimization design variables.
 * See the state! */
class CoefFunctionOpt : public CoefFunction, public enable_shared_from_this<CoefFunctionOpt>
{
public:

  /** One of the three states for the object:
   * OPT -> optimization and we ask DesignSpace::ApplyPhysicalDesign()
   * ORG -> the original material CoefFunction is used as if this would be a standard CoefFunction
   * SHADOW -> the temporary material CoefFunction is used similar to the org material
   * DIRECTION -> gradient of ParamMat to be set in DesignMaterial */
  typedef enum { OPT = 0, ORG = 1, SHADOW = 2, DIRECTION = 3 } State;

  /** @param orgMat the CoefFunctionConst that would be originally used to provide the material data. */
  CoefFunctionOpt(DesignSpace* design, PtrCoefFct orgMat, MaterialType mt, SinglePDE* pde);

  virtual ~CoefFunctionOpt() { }

  virtual string GetName() const { return "CoefFunctionOpt"; }

  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<double>& coefMat, const LocPointMapped& lpm) {
    GetTensor<double>(coefMat, lpm);
  }

  void GetTensor(Matrix<Complex>& coefMat, const LocPointMapped& lpm) {
    GetTensor<Complex>(coefMat, lpm);
  }

  //! \copydoc CoefFunction::GetScalar
  void GetScalar(double& scal, const LocPointMapped& lpm) {
    GetScalar<double>(scal, lpm);
  }

  //! \copydoc CoefFunction::GetScalar
  void GetScalar(Complex& scal, const LocPointMapped& lpm) {
    GetScalar<Complex>(scal, lpm);
  }

  //! \copydoc CoefFunction::GetVector
  void GetVector(Vector<Double>& vec, const LocPointMapped& lpm) {
    GetVector<double>(vec, lpm);
  }


  //! \copydoc CoefFunction::GetVecSize
  unsigned int GetVecSize() const {
    return orgMat->GetVecSize();
  }

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
    orgMat->GetTensorSize(numRows, numCols);
  }

  //! Set associated coordinate system
  virtual void SetCoordinateSystem(CoordSystem* cSys) {
    orgMat->SetCoordinateSystem(cSys);
  }

  //! \copydoc CoefFunction::IsZero
  bool IsZero() const  {
    // TODO: add optimization
    return orgMat->IsZero();
  }

  /** set the form such that the proper transfer function can be found */
  void SetForm(BiLinearForm* form) {
    this->form = form;
  }

  /** set the form such that the proper transfer function can be found */
  void SetForm(LinearForm* formL) {
    this->formL = formL;
  }

  /** only to query the name to find the proper transfer function */
  const BiLinearForm* GetForm() const  {
    return form;
  }

  /** only to query the name to find the proper transfer function */
  const LinearForm* GetFormL() const  {
    return formL;
  }

  MaterialType GetMaterialType() const { return materialType; }

  SinglePDE* GetPDE() const { return pde; }


  /** enable optimization, which means that the design space is asked for the the proper material. state -> OPT */
  void SetToOptimization();

  /** only the original material CoefFunction does the job. state -> ORG */
  void SetToOrgMaterial();

  /** No optimization but use the shadow material as long as the state is switched. state -> SHADOW.
  /** Necessary for OptimizationMaterial to find basic local element matrices. */
  void SetToShadow(PtrCoefFct shadow);

  /** the direction (silly name :( ) is gradient of a tensor parameterization with respect to a design type.
   * Evaluated by DesignMaterial */
  void SetToMaterialDerivative(DesignElement::Type direction);

  /** Is DesignElement::NO_DERIVATIVE is state is not DIRECTION. For tensors and mass */
  DesignElement::Type GetMaterialDerivative() const { return direction; }

  State GetState() const { return state; }

  /** the original material. Required if no optimization available or for SIMP and bi-material optimization */
  PtrCoefFct orgMat;

  SubTensorType subTensor;

protected:

  template <class T>
  void GetTensor(Matrix<T>& coefMat, const LocPointMapped& lpm);

  template <class T>
   void GetMsfemElementMatrix(Matrix<T>& elemMat, const LocPointMapped& lpm);

  template <class T>
  void GetScalar(T& scal, const LocPointMapped& lpm);

  template <class T>
  void GetVector(Vector<T>& scal, const LocPointMapped& lpm);

  /** the original material functions seem not to save what they actually are?! */
  MaterialType materialType = NO_MATERIAL;

  /** the pde we work on */
  SinglePDE* pde;

  /** This is the DesignSpace we use -> could be also requested from domain */
  DesignSpace* design;

  /** we store the form such that we can identify the proper transfer function */
  BiLinearForm* form;

  /** we store the form such that we can identify the proper transfer function */
  LinearForm* formL;

  /** be may switch the query of optimization off to get the real material */
  State state;

  /** see DisableOpt() */
  PtrCoefFct shadowMat;

  /** the current tensor derivative, only for state == DIRECTION */
  DesignElement::Type direction;
};

}

#endif /* COEFFUNCTIONOPT_HH_ */
