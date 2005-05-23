#ifndef FILE_DBINTCOLUMN_2004
#define FILE_DBINTCOLUMN_2004

#include "dbMatrix.hh"
#include "dbColumn.hh"
#include <vector>

namespace CoupledField
{

  class dbIntColumn : public dbColumn
  {
  public:
    //! Append value to the end of the integer-vector. 
    /*!
      \param v Value to append at the end
    */
    void append (int v);

    //! Get value at idx and store it in result. 
    /*!
      \param result Pointer, where to store the result
      \param idx Index of vector
    */
    void get (int &result, int idx);

    //! Get size of vector
    int size();

  private:
    //! Vector for one column of integer values
    std::vector<int> value;
  };

} // End of namespace

#endif
