#include "baseoperator.hh"
#include "Domain/grid.hh"
namespace CoupledField
{

  BaseOperator::BaseOperator(Grid * ptGrid, StdPDE * ptPDE, 
                             shared_ptr<EqnMap> eqnMap,
                             bool isaxi, bool coordUpdate )
    : isaxi_(isaxi)
  {
    ENTER_FCN( "BaseOperator::BaseOperator" );

    this->ptGrid_ = ptGrid;
    this->ptPDE_ = ptPDE;
    this->eqnMap_ = eqnMap;
    this->coordUpdate_ = coordUpdate;
  }


  BaseOperator::~BaseOperator()
  {
    ENTER_FCN( "BaseOperator::BaseOperator" );

    ;
  }

} // end of namespace
