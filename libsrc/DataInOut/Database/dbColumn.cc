#include "dbColumn.hh"

namespace CoupledField
{

  void dbColumn::append (int v)
  {
    Warning("Wrong call of append (int)",__FILE__,__LINE__);
  }

  void dbColumn::append (double v)
  {
    Warning("Wrong call of append (double),__FILE__,__LINE__");
  }

  void dbColumn::append (std::string v)
  {
    Warning("Wrong call of append (string)",__FILE__,__LINE__);
  }

  void dbColumn::get (int &result, int idx)
  {
    Warning("Wrong call of get (int*,int)",__FILE__,__LINE__);
  }

  void dbColumn::get (double &result, int idx)
  {
    Warning("Wrong call of get (double*,int)",__FILE__,__LINE__);
  }

  void dbColumn::get (std::string &result, int idx)
  {
    Warning("Wrong call of get (string*,int)",__FILE__,__LINE__);
  }

} // End of namespace
