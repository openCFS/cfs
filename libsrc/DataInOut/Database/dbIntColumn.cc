#include "dbIntColumn.hh"

namespace CoupledField
{

void dbIntColumn::append (int v)
{
  value.push_back(v);
}

void dbIntColumn::get (int &result, int idx)
{
  if (idx<size())
  {
    result = value[idx];
  }
  else
  {
#ifdef DEBUG
    (*debug)<<"dbIntColumn::get(int,int): returns -1 as idx="<<idx<<" and size()="<<size()<<std::endl;
#endif
    result = -1;
  }
}

int dbIntColumn::size()
{
  return value.size();
}

} // End of namespace
