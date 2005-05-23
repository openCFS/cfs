#include "dbLineData.hh"


namespace CoupledField
{

  dbLineData::dbLineData(const std::string &tname)
  {
    ENTER_FCN("dbLineData::dbLineData");
    SetTableName(tname);
    fields        = new char[STRINGBLOCKSIZE];
    strcpy(fields,"");
    fieldlength   = STRINGBLOCKSIZE;
    values        = new char[STRINGBLOCKSIZE];
    strcpy(values,"");
    valuelength   = STRINGBLOCKSIZE;
  }

  dbLineData::~dbLineData()
  {
    ENTER_FCN("dbLineData::~dbLineData");
    delete[] fields;
    delete[] values;
  }

  char *dbLineData::GetFieldNames()
  {
    ENTER_IFCN("dbLineData::GetFieldNames");
    return fields;
  }

  char *dbLineData::GetValues()
  {
    ENTER_IFCN("dbLineData::GetValues");
    return values;
  }

  void dbLineData::SetTableName(const std::string &tname)
  {
    ENTER_IFCN("dbLineData::SetTableName");
    tablename = tname;
  }

  void dbLineData::Clear()
  {
    ENTER_FCN("dbLineData::Clear");
    strcpy(values,"");
    strcpy(fields,"");
  }

  template<> void dbLineData::toString(const char *val, char* &dst, int &strsize)
  {
    ENTER_IFCN("dbLineData::toString(char*,char*)");
    if (strlen(val)>=strsize)
      {
        delete[] dst;
        strsize=strlen(val)+1;
        dst = new char[strsize];
      }
    strcpy(dst,val);
  }

  template<> void dbLineData::toString(std::string val, char* &dst, int &strsize)
  {
    ENTER_IFCN("dbLineData::toString(std::string,char*)");
    Integer vallen=val.length();
    if (vallen>=strsize)
      {
        delete[] dst;
        strsize=vallen+1;
        dst = new char[strsize];
      }
    strcpy(dst,val.c_str());
  }

  template<> void dbLineData::toString(Double val, char* &dst, int &strsize)
  {
    ENTER_IFCN("dbLineData::toString(Double,char*)");
    if (strsize<20)
      {
        delete[] dst;
        strsize = 25;
        dst = new char[strsize];
      }
    sprintf(dst,"%e",val);
  }

  template<> void dbLineData::toString(float val, char* &dst, int &strsize)
  {
    ENTER_IFCN("dbLineData::toString(float,char*)");
    if (strsize<20)
      {
        delete[] dst;
        strsize = 25;
        dst = new char[strsize];
      }
    sprintf(dst,"%e",val);
  }

  template<> void dbLineData::toString(Integer val, char* &dst, int &strsize)
  {
    ENTER_IFCN("dbLineData::toString(Integer,char*)");
    if (strsize<20)
      {
        delete[] dst;
        strsize = 25;
        dst = new char[strsize];
      }
    sprintf(dst,"%d",val);
  }


} // end of namespace
