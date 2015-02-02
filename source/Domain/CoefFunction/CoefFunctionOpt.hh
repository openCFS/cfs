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
 * is completely constructed out of optimization design variables. */
class CoefFunctionOpt : public CoefFunction, public boost::enable_shared_from_this<CoefFunctionOpt>
{
public:

  /** @param orgMat the CoefFunctionConst that would be originally used to provide the material data. */
  CoefFunctionOpt(DesignSpace* design, PtrCoefFct orgMat);

  virtual ~CoefFunctionOpt() { }

  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<double>& coefMat, const LocPointMapped& lpm);

  //! \copydoc CoefFunction::GetScalar
  void GetScalar(double& scal, const LocPointMapped& lpm);


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

  /** flag the optimization off */
  void SetEnable(bool val) { enabled = val; }

  /** the original material. Required if no optimization available or for SIMP and bi-material optimizaton */
  PtrCoefFct orgMat;

protected:

  /** This is the DesignSpace we use -> could be also requested from domain */
  DesignSpace* design;

  /** we store the form such that we can identify the proper transfer function */
  BiLinearForm* form;

  /** be may switch the query of optimization off to get the real material */
  bool enabled;
};

}

#endif /* COEFFUNCTIONOPT_HH_ */
