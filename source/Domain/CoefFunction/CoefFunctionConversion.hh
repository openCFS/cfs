// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
//==============================================================================
/*!
 *       \file     CoefFunctionConversion.hh
 *       \brief    Collection for CoefFunctions which are used to convert one base CoefFunction to one
 *                 of a different type. This works by re-implementing the query methods GetXXX in whicht the
 *                 base data from the original CoefFunction is manipulated. Possible use cases are e.g.
 *                 - Vector to Tensor (every vector is a tensor)
 *
 *                 Since all of this implementations should be pretty small, we do them in the header directly.
 *
 *       \date     10/02/2020
 *       \author   Florian Toth
 */
//==============================================================================

#ifndef COEF_FUNCTION_CONVERSION_HH
#define COEF_FUNCTION_CONVERSION_HH

#define _USE_MATH_DEFINES
#include <cmath>

#include <string>
#include <map>

#include "CoefFunction.hh"
#include "FeBasis/BaseFE.hh"
#include "Materials/BaseMaterial.hh"
#include "FeBasis/FeFunctions.hh"
#include "Utils/mathParser/mathParser.hh"


namespace CoupledField  {

template <class T>
class CoefFunctionVectorToTensor : public CoefFunction {

public:
  //! Constructor
  CoefFunctionVectorToTensor(PtrCoefFct coefToChange, bool transpose=false){
    assert( coefToChange->GetDimType() == VECTOR);
    coefToChange_ = coefToChange;

    // Variables from CoefFunction base-class
    isAnalytic_ = false;
    isComplex_ = false;

    dimType_ = TENSOR;
    len_ = coefToChange_->GetVecSize();
    transpose_= transpose;
    if (transpose_){
      nrows_ = 1;
      ncols_ = len_;
    } else {
      nrows_ = len_;
      ncols_ = 1;
    }
  }

  //! Destructor
  virtual ~CoefFunctionVectorToTensor(){;};


  //! \copydoc CoefFunction::GetScalar
  virtual void GetScalar(T& coefScal, const LocPointMapped& lpm ){
    EXCEPTION("Should never be called ...");
  }


  //! \copydoc CoefFunction::GetVector
  virtual void GetVector(Vector<T>& coefVec, const LocPointMapped& lpm){
    EXCEPTION("Should never be called ...");
  }

  //! \copydoc CoefFunction::GetTensor
  virtual void GetTensor(Matrix<T>& coefMat, const LocPointMapped& lpm){
    coefMat.Resize(nrows_,ncols_);
    Vector<T> vec;
    coefToChange_->GetVector(vec,lpm);
    if (transpose_){
      for (UInt i=0; i<len_; i++) {
          coefMat[0][i] = vec[i];
      }
    } else {
      for (UInt i=0; i<len_; i++) {
          coefMat[i][0] = vec[i];
      }
    }
  }

  //! \copydoc CoefFunction::ToString
  virtual std::string ToString() const{
    return coefToChange_->ToString();
  };

  //! Return row and columns size of tensor if coefficient function is a tensor
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
      numRows = nrows_;
      numCols = ncols_;
  };

private:

  //! CoefFunction to change it's type
  PtrCoefFct coefToChange_ = NULL;

  //! Dimensions of the tensor
  UInt nrows_ = 0;
  UInt ncols_ = 0;
  UInt len_ = 0;

  //! Transpose or not (not transposed means we create a column vector with ncols_=len_ columns and one row);
  bool transpose_ = false;
};

} //end of namespace

#endif
