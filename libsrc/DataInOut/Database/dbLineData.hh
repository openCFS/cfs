#ifndef FILE_DBINSERTLINEDATA_2004 
#define FILE_DBINSERTLINEDATA_2004


#include <vector>
#include <string>
#include "General/environment.hh"

#define STRINGBLOCKSIZE	 1024
#define ELEMENTBLOCKSIZE   32


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

  //! Default destructor
  ~dbLineData();

  //! Set field to value
  /*!
    \param fieldname Field name
    \param val Value
    \param quoted If value is sent to database quoted or unquoted
  */
  template<class T> int Set (const char *fieldname, 
                             const T val, 
                             Boolean quoted=TRUE)
  {
    ENTER_IFCN("dbLineData::Set"); 
    Integer fnamelength = strlen(fieldname);
    Integer curFLength  = strlen(fields);
    if (curFLength+fnamelength>=fieldlength)	// Allocate new memory
    {
#ifdef DEBUG
    (*debug)<<"dbLineData::Set : Allocate new memory"<<std::endl;
#endif
      char *newFields = new char[2*fieldlength];
      strcpy(newFields,fields);
      delete[] fields;
      fields = newFields;
    }
#ifdef DEBUG2
    (*debug)<<"dbLineData::Set : Fields="<<fields<<" (before) "<<std::endl;
#endif
    if (curFLength>0)
      strcat(fields,",");
    strcat(fields,fieldname); 
#ifdef DEBUG2
    (*debug)<<"dbLineData::Set : Fields="<<fields<<std::endl;
#endif

    char *valstr = new char[ELEMENTBLOCKSIZE];
    toString(val,valstr);
    Integer vallength = strlen(valstr);
    Integer curVLength= strlen(values);
    if (curVLength+vallength>=valuelength)	// Allocate new memory
    {
      char *newValues = new char[2*valuelength];
      strcpy(newValues,values);
      delete[] values;
      values = newValues;
    }
    if (curVLength>0)
      strcat(values,",");
    if (quoted)
      strcat(values,"'");
    strcat(values,valstr);
    delete[] valstr;
    if (quoted)
      strcat(values,"'");
#ifdef DEBUG2
    (*debug)<<"dbLineData::Set : Values="<<values<<std::endl;
#endif

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
  char *GetFieldNames();

  //! Get value (as string)
  /*!
    \param idx Index
  */
  char *GetValues ();

  //! Field names
//  std::vector<std::string> field;
  char *fields;

  //! Values
//  std::vector<std::string> value;
  char *values;

  template<class T> void toString(const T val, char *dst)
  {
    ENTER_IFCN("dbLineData::toString(class T, char*)");
    std::stringstream sstream;
    sstream<<val;
    strcpy(dst,sstream.str().c_str());
  };


  int fieldlength;
  int valuelength;

private:
  
  dbLineData(){};

  dbLineData(const dbLineData &x)
  {
    ENTER_FCN("dbLineData::copy-constructor");
  };
};

} // end of namespace

#endif
