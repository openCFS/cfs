#include "dbStringColumn.hh"

namespace CoupledField
{

  void dbStringColumn::append (std::string v)
  {
    value.push_back(v);
  }

  void dbStringColumn::get (std::string &result, int idx)
  {
    result = value[idx];
  }

  int dbStringColumn::size()
  {
    return value.size();
  }

} // end of namespace
