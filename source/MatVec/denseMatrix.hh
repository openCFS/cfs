// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_DENSEMATRIX_HH
#define FILE_DENSEMATRIX_HH
 
#include <string>

namespace CoupledField
{
  //! Abstract interface class for a general dense matrix
  /*! Abstract interface class for a general dense matrix
    \note Although this class provides a general interface
    to the matrix class , one should always prefer a cast into
    the current type of matrix, e.g. 
    dynamic_cast<Matrix<Double>*>(ptDenseMatrix)
  */
  class DenseMatrix
  {
  public:
    //! Constructor 
    DenseMatrix() {};

    //! Destructor
    virtual ~DenseMatrix() {};

    /** @see Matrix::ToString() */
    virtual std::string ToString(const int level, const bool newline) const = 0;
    
    /** @see Matrix::ToXML() */
    virtual std::string ToXML(int offset = 0) const = 0;
  };
} // namespace
#endif
