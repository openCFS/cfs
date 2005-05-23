#ifndef FILE_DBMATRIX_2004
#define FILE_DBMATRIX_2004

#include "General/environment.hh"    // for traces and debug information only
#include "dbColumn.hh"
#include "dbIntColumn.hh"
#include "dbDoubleColumn.hh"
#include "dbStringColumn.hh"
#include <vector>
#include <iostream>
#include <iomanip>

#define DBFIELDTYPEINT 1
#define DBFIELDTYPEDOUBLE 2
#define DBFIELDTYPESTRING 3

namespace CoupledField
{

  class dbMatrix
  {
  public:
    //! Default destructor, frees table
    ~dbMatrix();

    //! Appends one column for integer values
    /*!
      \param colname Name of column
    */
    void appendIntColumn(std::string colname);

    //! Appends one column for double values
    /*!
      \param colname Name of column
    */
    void appendDoubleColumn(std::string colname);

    //! Appends one column for string values
    /*!
      \param colname Name of column
    */
    void appendStringColumn(std::string colname);

    //! Returns column of table with given column name
    /*!
      \param colname column name
    */
    dbColumn *operator[] (std::string colname) ;

    //! Returns i-th column of table
    /*!
      \param i i-th column
    */
    dbColumn *operator[] (const int i) ;

    //! Returns type of i-th column
    /*!
      \param i i-th column
    */
    int getFieldType (int i);

    //! Returns name of i-th column
    /*!
      \param i i-th column
    */
    std::string getFieldName (int i);

    //! Get number of columns
    int getNoOfColumn();

    //! Get number of rows
    int getNoOfRow();

    //! print matrix in debug file (only #ifdef DEBUG)
    void printMatrixDebug();
  

  private:
    //! Representation of table from database
    std::vector <dbColumn*> table;

    //! Names of columns
    std::vector <std::string> field_name;

    //! Data type of columns
    std::vector <int> field_type;
  };

} // End of namespace

#endif
