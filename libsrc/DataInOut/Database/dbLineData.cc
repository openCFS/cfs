#include "dbLineData.hh"


namespace CoupledField
{

dbLineData::dbLineData(const std::string &tname)
{
  ENTER_FCN("dbLineData::dbLineData");
  SetTableName(tname);
  fields 	= new char[STRINGBLOCKSIZE];
  strcpy(fields,"");
  fieldlength	= STRINGBLOCKSIZE;
  values 	= new char[STRINGBLOCKSIZE];
  strcpy(values,"");
  valuelength	= STRINGBLOCKSIZE;
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

/*int dbLineData::Size()
{
  return field.size();
}*/

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
//  value.clear();
//  field.clear();
}

template<> void dbLineData::toString(const char *val, char *dst)
{
  ENTER_IFCN("dbLineData::toString(char*,char*)");
  if (strlen(val)>ELEMENTBLOCKSIZE)
  {
    delete dst;
    dst = new char[strlen(val)+1];
  }
  strcpy(dst,val);
}

template<> void dbLineData::toString(std::string val, char *dst)
{
  ENTER_IFCN("dbLineData::toString(std::string,char*)");
  Integer vallen=strlen(val.c_str());
  if (vallen>ELEMENTBLOCKSIZE)
  {
    delete dst;
    dst = new char[vallen+1];
  }
  strcpy(dst,val.c_str());
}

template<> void dbLineData::toString(Double val, char *dst)
{
  ENTER_IFCN("dbLineData::toString(Double,char*)");
  sprintf(dst,"%e",val);
}

template<> void dbLineData::toString(float val, char *dst)
{
  ENTER_IFCN("dbLineData::toString(float,char*)");
  sprintf(dst,"%e",val);
}

template<> void dbLineData::toString(Integer val, char *dst)
{
  ENTER_IFCN("dbLineData::toString(Integer,char*)");
  sprintf(dst,"%d",val);
}


} // end of namespace
