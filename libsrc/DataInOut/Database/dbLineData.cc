#include "dbLineData.hh"


namespace CoupledField
{

dbLineData::dbLineData(const std::string &tname)
{
  SetTableName(tname);
}

std::string dbLineData::GetFieldName(int i)
{
  return field[i];
}

std::string dbLineData::GetValue(int i)
{
  return value[i];
}

int dbLineData::Size()
{
  return field.size();
}

void dbLineData::SetTableName(const std::string &tname)
{
  tablename = tname;
}

void dbLineData::Clear()
{
  value.clear();
  field.clear();
}

} // end of namespace
