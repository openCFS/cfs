#include "dbDoubleColumn.hh"

namespace CoupledField
{

  void dbDoubleColumn::append (int v)
  {
    value.push_back(v);
  }

  void dbDoubleColumn::append (double v)
  {
    value.push_back(v);
  }

  void dbDoubleColumn::get (double &result, int idx)
  {
    result = value[idx];
  }

  int dbDoubleColumn::size()
  {
    return value.size();
  }

} // End of namespace
