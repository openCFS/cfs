#ifndef COEFFUNCTIONOPT_HH_
#define COEFFUNCTIONOPT_HH_

#include <boost/tr1/type_traits.hpp>
#include <boost/shared_ptr.hpp>

#include "CoefFunction.hh"


namespace CoupledField {

class DesignSpace;
class BiLinearForm;

/** Replaces a material data CoefFunctionConst with an ersatz material optimization version.
 * In the SIMP case this the original const material is scaled, in the bi-material it is interpolated with an additional
 * material, in multi-material it is the weighted sum of various materials and in parametric optimization the material (tensor)
 * is completely constructed out of optimization design variables.
 * See the state! */
class CoefFunctionOpt : public CoefFunction, public boost::enable_shared_from_this<CoefFunctionOpt>
{
public:

  /** One of the three states for the object:
   * OPT -> optimization and we ask DesignSpace::ApplyPhysicalDesign()
   * ORG -> the original material CoefFunction is used as if this would be a standard CoefFunction
   * SHADOW -> the temporary material CoefFunction is used similar to the org material */
  typedef enum { OPT = 0, ORG = 1, SHADOW = 2 } State;

  /** @param orgMat the CoefFunctionConst that would be originally used to provide the material data. */
  CoefFunctionOpt(DesignSpace* design, PtrCoefFct orgMat);

  virtual ~CoefFunctionOpt() { }

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

  /** only to query the name to finde the proper transfer function */
  const BiLinearForm* GetForm() const  {
    return form;
  }

  /** enable optimization, which means that the design space is asked for the the proper material. state -> OPT */
  void SetToOptimization();

  /** only the original material CoefFunction does the job. state -> ORG */
  void SetToOrgMaterial();

  /** No optimization but use the shado material as long as the state is switched. state -> SHADOW.
  /** Necessary for OptimizationMaterial to find basic local element matrices. */
  void SetToShadow(PtrCoefFct shadow);

  /** the original material. Required if no optimization available or for SIMP and bi-material optimization */
  PtrCoefFct orgMat;

protected:

  template <class T>
  void GetTensor(Matrix<T>& coefMat, const LocPointMapped& lpm);

  template <class T>
  void GetScalar(T& scal, const LocPointMapped& lpm);

  /** This is the DesignSpace we use -> could be also requested from domain */
  DesignSpace* design;

  /** we store the form such that we can identify the proper transfer function */
  BiLinearForm* form;

  /** be may switch the query of optimization off to get the real material */
  State state;

  /** see DisableOpt() */
  PtrCoefFct shadowMat;
};

}

#endif /* COEFFUNCTIONOPT_HH_ */
