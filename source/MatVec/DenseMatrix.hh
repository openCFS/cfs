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
    virtual std::string ToXMLFormat(const std::string& name, const int offset) const = 0;
  };
} // namespace
#endif
