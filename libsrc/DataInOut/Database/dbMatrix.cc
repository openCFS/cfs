#include "dbMatrix.hh"

namespace CoupledField
{

dbColumn * dbMatrix::operator[] (const int i) 
{
  if (i>getNoOfColumn())
    Error("dbMatrix::operator[] (int): index too high",__FILE__,__LINE__);
  return table[i];
}

dbColumn * dbMatrix::operator[] (std::string colname)
{
  for (int i=0; i<getNoOfColumn(); i++)
  {
    if (colname==getFieldName(i))
    {
      return table[i];
    }
  }
  std::string err1("dbMatrix::operator[] (string) : Column '");
  err1+=colname; err1+="' not found in matrix";
  Error(err1.c_str(),__FILE__,__LINE__);
}

dbMatrix::~dbMatrix()
{
  int i;
  for (i=0; i<table.size(); i++)
    delete table[i];
  table.clear();
}

void dbMatrix::appendIntColumn(std::string colname)
{
  ENTER_FCN("dbMatrix::appendIntColumn");
  dbColumn* col;
  col = new dbIntColumn;
  table.push_back(col);
  field_type.push_back(DBFIELDTYPEINT);	
  field_name.push_back(colname);
}

void dbMatrix::appendDoubleColumn(std::string colname)
{
  ENTER_FCN("dbMatrix::appendDoubleColumn");
  dbColumn* col;
  col = new dbDoubleColumn();
  table.push_back(col);
  field_type.push_back(DBFIELDTYPEDOUBLE);	
  field_name.push_back(colname);
}

void dbMatrix::appendStringColumn(std::string colname)
{
  ENTER_FCN("dbMatrix::appendStringColumn");
  dbColumn* col;
  col = new dbStringColumn();
  table.push_back(col);
  field_type.push_back(DBFIELDTYPESTRING);
  field_name.push_back(colname);
}

void dbMatrix::printMatrixDebug()
{
#ifdef DEBUG
  for (int i=0; i<getNoOfColumn(); i++)
  {
    (*debug)<< std::setw(10) << getFieldName(i);
  }
  (*debug)<< std::endl<<std::endl;
  for (int j=0; j<getNoOfRow(); j++)
  {
    for (int i=0; i<getNoOfColumn(); i++)
    {
      if (getFieldType(i)==DBFIELDTYPEINT)
      {
        int val;
        table[i]->get(val,j);
        (*debug)<< std::setw(10) <<val;
      }
      else if (getFieldType(i)==DBFIELDTYPEDOUBLE)
      {
        double val;
        table[i]->get(val,j);
        (*debug)<< std::setw(10) <<val;
      }
      else if (getFieldType(i)==DBFIELDTYPESTRING)
      {
        std::string val;
        table[i]->get(val,j);
        (*debug)<< std::setw(10) <<val;
      }
      else
      {
        Warning ("dbMatrix::printMatrixDebug : Unknown DBFIELDTYPE",__FILE__,__LINE__);
        std::string val("Error");
        (*debug)<< std::setw(10) <<val;
      }
    }
  }
  (*debug)<<std::endl;

#endif
}

int dbMatrix::getFieldType(int i)
{
  return field_type[i];
}

std::string dbMatrix::getFieldName(int i)
{
  return field_name[i];
}

int dbMatrix::getNoOfColumn()
{
  return table.size();
}

int dbMatrix::getNoOfRow()
{
  if (getNoOfColumn()>0)
    return table[0]->size();
  else
    return 0;
}

} // End of namespace
