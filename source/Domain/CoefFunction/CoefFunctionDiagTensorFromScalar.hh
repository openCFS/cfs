#ifndef COEFFUNCTION_TENSOR_FROM_SCALAR_HH
#define COEFFUNCTION_TENSOR_FROM_SCALAR_HH

//#include <boost/tr1/type_traits.hpp>
#include <boost/shared_ptr.hpp>

#include "CoefFunction.hh"


namespace CoupledField {

//! Provide a diagonal tensorial coefficient function defined by scalar ones

//! This class generates a tensorial coefficient function, where the diagonal
//! elements are defined in terms of scalar-valued coefficient functions.
class CoefFunctionDiagTensorFromScalar :
   public CoefFunction,
   public boost::enable_shared_from_this<CoefFunctionDiagTensorFromScalar >
{
public:

  //! Constructor
  CoefFunctionDiagTensorFromScalar(const StdVector<PtrCoefFct>& scalVals);

     //! Destructor
  virtual ~CoefFunctionDiagTensorFromScalar(){
    ;
  }

  virtual string GetName() const { return "CoefFunctionDiagTensorFromScalar"; }


  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<Double>& coefMat,
                 const LocPointMapped& lpm ) {
    coefMat.Resize(size_, size_);
    coefMat.Init();
    for(UInt i = 0; i < size_; ++i) {
      scalVals_[i]->GetScalar(coefMat[i][i], lpm);
    }
  }

  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<Complex>& coefMat,
                 const LocPointMapped& lpm ) {
    coefMat.Resize(size_, size_);
    coefMat.Init();
    for(UInt i = 0; i < size_; ++i) {
      scalVals_[i]->GetScalar(coefMat[i][i], lpm);
    }
  }

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
    assert(this->dimType_ == TENSOR );
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
};

}
#endif
