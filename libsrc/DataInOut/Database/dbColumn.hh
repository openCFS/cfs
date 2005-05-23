#ifndef FILE_DBCOLUMN_2004
#define FILE_DBCOLUMN_2004

#include <iostream>
#include "General/environment.hh"   // for debug and trace information only
#include "Utils/tools.hh"

namespace CoupledField
{

  class dbColumn
  {
  public:
    //! Append value to the end of the integer-vector. Implemented in dbIntColumn and dbDoubleColumn only.
    /*!
      \param v Value to append at the end
    */
    virtual void append (int v);   

    //! Append value to the end of the double-vector. Implemented in dbDoubleColumn only.
    /*!
      \param v Value to append at the end
    */
    virtual void append (double v);

    //! Append string to the end of the string-vector. Implemented in dbStringColumn only.
    /*!
      \param v Value to append at the end
    */
    virtual void append (std::string v);

    //! Get value at idx and store it in result. Implemented in dbIntColumn only.
    /*!
      \param result [Out] Pointer, where to store the result
      \param idx [In] Index of vector
    */
    virtual void get (int &result, int idx);

    //! Get value at idx and store it in result. Implemented in dbDoubleColumn only.
    /*!
      \param result [Out] Pointer, where to store the result
      \param idx [In] Index of vector
    */
    virtual void get (double &result, int idx);

    //! Get value at idx and store it in result. Implemented in dbStringColumn only.
    /*!
      \param result [Out] Pointer, where to store the result
      \param idx [In] Index of vector
    */
    virtual void get (std::string &result, int idx);

    //! Get size of vector
    virtual int size()=0;
  };

} // End of namespace

#endif
