#ifndef FILE_DATABASE
#define FILE_DATABASE

#include <mysql/mysql.h>
#include "DataInOut/writeresults.hh"
#include "DataInOut/Database/dbMatrix.hh"
#include "DataInOut/Database/dbLineData.hh"
#include <sstream>
#include <string>
#include <stdio.h>
#include <errno.h>

#define SIZEOFPARAMBLOCK 4096    // 4kB

namespace CoupledField
{

class Database
{
public:
  //! Default constructor
  Database();

  //! Default destructor
  ~Database();

  //! Connect to specified database
  /*!
    \param hostname Hostname of database server
    \param port TCP/IP-Port of database server
    \param user Username of MySQL-user
    \param passwd Password of chosen user
    \param db Name of database
  */
  int Connect (std::string hostname, unsigned int port, std::string user, std::string passwd, std::string db);

  //! Connect to database defined before
  int Connect ();

  //! Close connection to database
  void Close ();

  //! Send query to database
  /*!
    \param statement SQL-statement
  */
  int Query (std::string statement);

  //! Is connected to database?
  Boolean IsConnected;

  //! Insert data in specified table
  /*!
    \param d Dataset with the data to insert into table
  */
  int Insert (dbLineData &d);

  //! Insert data and get Index
  /*!
    \param d Dataset with the data to insert into table
  */
  int InsertAndGetIndex (dbLineData &d);

  //! Convert string to escaped string (use for strings with special characters
  /*!
    \param source String to mask with escape characters if needed
  */
  std::string EscapeString(std::string source);

  //! Get table with given name from database
  /*!
    \param matrix [Out] Matrix to store the table
    \param tablename [In] Name of table
  */
  Boolean GetTable(dbMatrix &matrix, std::string tablename);

  //! Get number of fields of last query
  unsigned int GetNoOfResultFields();

  //! Get number of resulting lines of last query
  int GetNoOfResultLines();

  //! Select fields from table. Call FetchFields afterwards to get the table
  /*!
    \param fields Fields to receive from database (comma separated, * for all)
    \param table Table name
    \param where Special selection (e.g. 'idx=123')
  */
  Boolean SelectFrom(std::string fields, std::string table, std::string where);

  //! Fetch results of last SQL statement
  /*!
    \param matrix [Out] Matrix to store the table
  */
  Boolean FetchFields(dbMatrix &matrix);

  //! Update selected fields
  /*!
    \param table Table name
    \param d New data
    \param wherestr WHERE-string, selection within table
  */
//  Boolean Update (dbLineData &d, dbLineData &where);

  //! Update selected fields
  /*!
    \param table Table name
    \param set Variables to set (e.g: "value=3")
    \param where Selection (e.g: "idx=123")
  */
  Boolean Update (std::string table,
                  std::string set,
                  std::string where);

protected:
  //! Connection to the database
  MYSQL *Conn_;

  //! Hostname of database server
  std::string Hostname_;

  //! TCP/IP-Port of database server
  unsigned int Port_;

  //! Name of database user
  std::string User_;

  //! Password of database user
  std::string Passwd_;

  //! Database name
  std::string Db_;

  //! Are the results of the last query stored?
  Boolean ResultStored_;

  //! Results of last query
  MYSQL_RES *Result_;

  //! Store the results of the last query
  Boolean StoreResults();

  //! Field information of last query
  MYSQL_FIELD *Fields_;

  //! Are the field information stored?
  Boolean FieldsStored_;

  //! Store field information of last query
  Boolean StoreFields();

  //! Free memory used for results
  void FreeResult();

}; // End of class Database

} // End of namespace

#endif
