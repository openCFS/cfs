#ifndef FILE_DBINSERTLINEDATA_2004 
#define FILE_DBINSERTLINEDATA_2004


#include <vector>
#include <string>
#include "General/environment.hh"
#include "DataInOut/WriteInfo.hh"

#define STRINGBLOCKSIZE  1024
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
      if (curFLength+fnamelength>=fieldlength)    // Allocate new memory
        {
          char *newFields;
          if (strlen(fieldname)<fieldlength)
            {
              newFields = new char[2*fieldlength];
              fieldlength=fieldlength*2;
            }
          else
            {
              newFields = new char[fieldlength+2*strlen(fieldname)];
              fieldlength += 2*strlen(fieldname);
            }
          strcpy(newFields,fields);
          delete[] fields;
          fields = newFields;
        }
      if (curFLength>0)
        strcat(fields,",");
      strcat(fields,fieldname); 

      Integer valstrlen = ELEMENTBLOCKSIZE;
      char *valstr = new char[valstrlen];
      toString(val,valstr,valstrlen);
      Integer vallength = strlen(valstr);
      Integer curVLength= strlen(values);
      if (curVLength+vallength>=valuelength)      // Allocate new memory
        {
          char *newValues; 
          if (strlen(valstr)<valuelength)
            {
              newValues = new char[2*valuelength];
              valuelength=2*valuelength;
            }
          else
            {
              newValues = new char[valuelength+2*strlen(valstr)];
              valuelength += 2*strlen(valstr);
            }
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
    };

    //! Get size of vector
    int Size();

    //! Clear vector
    void Clear();

    //! Set table name
    void SetTableName(const std::string &tname);

    //! Table name 
    std::string   tablename;

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
    char *fields;

    //! Values
    char *values;

    //! Convert any type to a C-string
    /*!
      \param val [In] Value to convert
      \param dst [Out] Where to store the string
    */
    template<class T> void toString(const T val, char* &dst, int &strsize)
    {
      ENTER_IFCN("dbLineData::toString(class T, char*)");
      std::stringstream sstream;
      sstream<<val;
      Integer len = sstream.str().length();
      if (len>=strsize)
        {
          delete[] dst;
          strsize = len+1;
          dst = new char[strsize];
        }
      strcpy(dst,sstream.str().c_str());
    };


    //! Length of fields
    int fieldlength;
 
    //! Length of values
    int valuelength;

  private:
  
    //! Empty default constructor - use dbLineData(std::string) instead
    dbLineData(){};

    //! Copy constructor - not implemented - should not be used at all
    dbLineData(const dbLineData &x)
    {
      ENTER_FCN("dbLineData::copy-constructor");
      Error("dbLineData copy constructor not implemented - should not be used at all",__FILE__,__LINE__);
    };
  };

} // end of namespace

#endif
