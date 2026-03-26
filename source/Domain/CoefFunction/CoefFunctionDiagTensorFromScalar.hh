#ifndef COEFFUNCTION_TENSOR_FROM_SCALAR_HH
#define COEFFUNCTION_TENSOR_FROM_SCALAR_HH

// #include <boost/tr1/type_traits.hpp>
#include <memory>

#include "CoefFunction.hh"

namespace CoupledField {

//! Provide a diagonal tensorial coefficient function defined by scalar ones

//! This class generates a tensorial coefficient function, where the diagonal
//! elements are defined in terms of scalar-valued coefficient functions.
class CoefFunctionDiagTensorFromScalar : public CoefFunction,
                                         public enable_shared_from_this<CoefFunctionDiagTensorFromScalar> {
public:
  //! Constructor
  CoefFunctionDiagTensorFromScalar(const StdVector<PtrCoefFct> &scalVals, std::string subType = "");

  //! Destructor
  virtual ~CoefFunctionDiagTensorFromScalar()
  {
    ;
  }

  virtual string GetName() const { return "CoefFunctionDiagTensorFromScalar"; }

  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<Double> &coefMat, const LocPointMapped &lpm)
  {
    GetTensor_<Double>(coefMat, lpm);
  }
  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<Complex> &coefMat, const LocPointMapped &lpm)
  {
    GetTensor_<Complex>(coefMat, lpm);
  }

  //! \copydoc CoefFunction::GetTensor
  template <typename TYPE>
  void GetTensor_(Matrix<TYPE> &coefMat, const LocPointMapped &lpm)
  {
    coefMat.Resize(size_, size_);
    coefMat.Init();
    for (UInt i = 0; i < size_; ++i) {
      scalVals_[i]->GetScalar(coefMat[i][i], lpm);
    }
  }

  //! \copydoc CoefFunction::GetVector
  void GetVector(Vector<Double> &coefVec, const LocPointMapped &lpm)
  {
    GetVector_<Double>(coefVec, lpm);
  }
  //! \copydoc CoefFunction::GetVector
  void GetVector(Vector<Complex> &coefVec, const LocPointMapped &lpm)
  {
    GetVector_<Complex>(coefVec, lpm);
  }

  //! \copydoc CoefFunction::GetVector
  template <typename TYPE>
  void GetVector_(Vector<TYPE> &coefVec, const LocPointMapped &lpm)
  {
    if (subType_ == "3d") {
      // Components: "xx", "yy", "zz", "yz", "xz", "xy"
      coefVec.Resize(size_, 1);
      coefVec.Init();
      scalVals_[0]->GetScalar(coefVec[0], lpm);
      scalVals_[1]->GetScalar(coefVec[1], lpm);
      scalVals_[2]->GetScalar(coefVec[2], lpm);
    }
    else if (subType_ == "plane") {
      // Components: "xx", "yy", "xy"
      coefVec.Resize(size_, 1);
      coefVec.Init();
      scalVals_[0]->GetScalar(coefVec[0], lpm);
      scalVals_[1]->GetScalar(coefVec[1], lpm);
    }
    else if (subType_ == "axi") {
      // Components = "rr", "zz", "rz", "phiphi"
      coefVec.Resize(size_, 1);
      coefVec.Init();
      scalVals_[0]->GetScalar(coefVec[0], lpm);
      scalVals_[1]->GetScalar(coefVec[1], lpm);
      scalVals_[3]->GetScalar(coefVec[3], lpm);
    }
    else {
      EXCEPTION("Unkown subtype, can't convert to voigt notation")
    }
  }

  //! \copydoc CoefFunction::GetVecSize
  virtual UInt GetVecSize() const
  {
    assert(this->dimType_ == VECTOR);
    return size_;
  }

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize(UInt &numRows, UInt &numCols) const
  {
    assert(this->dimType_ == TENSOR);
    numRows = size_;
    numCols = size_;
  }

  //! \copydoc CoefFunction::IsZero
  bool IsZero() const;

  //! \copydoc CoefFunction::ToString
  std::string ToString() const;

protected:
  //! Size of square tensor
  UInt size_;

  //! Vector with diagonal scalar values
  StdVector<PtrCoefFct> scalVals_;

  //! Subtype of the solution
  std::string subType_;
};

} // namespace CoupledField
#endif
