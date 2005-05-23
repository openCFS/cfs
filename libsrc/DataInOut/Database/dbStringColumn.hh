#ifndef FILE_DBSTRINGCOLUMN_2004
#define FILE_DBSTRINGCOLUMN_2004

#include "dbColumn.hh"
#include "dbMatrix.hh"
#include <vector>

namespace CoupledField
{

  class dbStringColumn : public dbColumn
  {
  public:
    //! Append string to the end of the string-vector. 
    /*!
      \param v Value to append at the end
    */
    void append (std::string v);

    //! Get value at idx and store it in result. Implemented in dbStringColumn only.
    /*!
      \param result Pointer, where to store the result
      \param idx Index of vector
    */
    void get (std::string &result, int idx);

    //! Return size of vector
    int size();

  private:
    //! Vector for one column of string values
    std::vector<std::string> value;
  };

} // End of namespace

#endif
