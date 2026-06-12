#ifndef FILE_DENSEMATRIX_HH
#define FILE_DENSEMATRIX_HH
 
#include <string>
#include "General/EnvironmentTypes.hh"

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
    virtual std::string ToString(ToStringFormat format = TS_PYTHON, const std::string& linesep="", int digits=-1) const = 0;

    
    /** @see Matrix::ToXML() */
    virtual std::string ToXMLFormat(const std::string& name, const int offset,
                                    bool typeTag = true, const std::string& extra_attr = "") const = 0;
  };
} // namespace
#endif
