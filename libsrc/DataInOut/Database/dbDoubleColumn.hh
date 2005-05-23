#ifndef FILE_DBDOUBLECOLUMN_2004
#define FILE_DBDOUBLECOLUMN_2004

#include "dbColumn.hh"
#include "dbMatrix.hh"
#include <vector>

namespace CoupledField
{

  class dbDoubleColumn :  public dbColumn
  {
  public:
    //! Append value to the end of the integer-vector. 
    /*!
      \param v Value to append at the end
    */
    void append (int v);

    //! Append value to the end of the double-vector. 
    /*!
      \param v Value to append at the end
    */
    void append (double v);

    //! Get value at idx and store it in result. 
    /*!
      \param result Pointer, where to store the result
      \param idx Index of vector
    */
    void get (double &result, int idx);

    //! Return size of vector
    int size();

  private:
    //! Vector for one column of double values
    std::vector<double> value;
  };

} // End of namespace

#endif
