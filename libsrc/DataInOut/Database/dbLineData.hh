#ifndef FILE_DBINSERTLINEDATA_2004 
#define FILE_DBINSERTLINEDATA_2004


#include <vector>
#include <string>


namespace CoupledField
{

class dbLineData
{
public: 
  //! Default constructor
  /*! 
    \param tname Table name
  */
  dbLineData(const std::string &tname);

  //! Set field to value
  /*!
    \param fieldname Field name
    \param val Value
  */
  template<class T> int Set (const char *fieldname, const T val)
  {
    field.push_back(std::string(fieldname));
    std::stringstream v;
    v<<"'"<<val<<"'";
    value.push_back(v.str());

  };

  //! Set field to value and do not quote value. Use it only if you give functions in 'val'
  /*!
    \param fieldname Field name
    \param val Value
  */
  template<class T> int SetUnquoted (const char *fieldname, const T val)
  {
    field.push_back(std::string(fieldname));
    std::stringstream v;
    v<<val;
    value.push_back(v.str());

  };

  //! Get size of vector
  int Size();

  //! Clear vector
  void Clear();

  //! Set table name
  void SetTableName(const std::string &tname);

  //! Table name 
  std::string	tablename;

  //! Get field name
  /*!
    \param idx Index
  */
  std::string GetFieldName(int idx);

  //! Get value (as string)
  /*!
    \param idx Index
  */
  std::string GetValue (int idx);

private:
  //! Field names
  std::vector<std::string> field;

  //! Values
  std::vector<std::string> value;
};

} // end of namespace

#endif
