#ifndef COEFFUNCTIONOPT_HH_
#define COEFFUNCTIONOPT_HH_

#include <boost/tr1/type_traits.hpp>
#include <boost/shared_ptr.hpp>

#include "CoefFunction.hh"


namespace CoupledField {

/** Replaces a material data CoefFunctionConst with an ersatz material optimization version.
 * In the SIMP case this the original const material is scaled, in the bi-material it is interpolated with an additional
 * material, in multi-material it is the weighted sum of various materials and in parametric optimization the material (tensor)
 * is completely constructed out of optimization design variables. */
template<typename T>
class CoefFunctionOpt : public CoefFunction,
                          public boost::enable_shared_from_this<CoefFunctionOpt<T> >
{
public:

  /** @param orgMat the CoefFunctionConst that would be originally used to provide the material data. */
  CoefFunctionOpt(PtrCoefFct orgMat) : CoefFunction()
  {
    // this type of coefficient is always constant
    dependType_ = GENERAL;
    isAnalytic_ = false;
    isComplex_ = std::tr1::is_same<T,Complex>::value;
    supportDerivative_ = false;
    dimType_ = orgMat->GetDimType();

    this->orgMat = orgMat;
  }

  virtual ~CoefFunctionOpt() { }

  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<T>& coefMat, const LocPointMapped& lpm )
  {
    assert(this->dimType_ == TENSOR);
    // if no coordinate system is set, just
    // use internal vector
    if( !coordSys_ ) {
      orgMat->GetTensor(coefMat, lpm);
    }
    else
      EXCEPTION("the rotation is not fully finished ':-(\n");

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

protected:

  /** the original material. Required if no optimization available or for SIMP and bi-material optimizaton */
  PtrCoefFct orgMat;
};

}

#endif /* COEFFUNCTIONOPT_HH_ */
