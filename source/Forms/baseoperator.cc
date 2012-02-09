// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "baseoperator.hh"

namespace CoupledField {
class EqnMap;
}  // namespace CoupledField

namespace CoupledField
{

  BaseOperator::BaseOperator(Grid * ptGrid, StdPDE * ptPDE, 
                             boost::shared_ptr<EqnMap> eqnMap,
                             bool isaxi, bool coordUpdate )
    : isaxi_(isaxi)
  {

    this->ptGrid_ = ptGrid;
    this->ptPDE_ = ptPDE;
    this->eqnMap_ = eqnMap;
    this->coordUpdate_ = coordUpdate;
  }

} // end of namespace
