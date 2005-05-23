#include "baseoperator.hh"
#include "Domain/grid.hh"
namespace CoupledField
{

BaseOperator::BaseOperator(Grid * ptGrid, StdPDE * ptPDE, 
			   NodeEQN * ptEQN, 
			   Boolean isaxi)
  : isaxi_(isaxi)
{
  ENTER_FCN( "BaseOperator::BaseOperator" );

  this->ptGrid_ = ptGrid;
  this->ptPDE_ = ptPDE;
  this->ptEQN_ = ptEQN;
}


BaseOperator::~BaseOperator()
{
  ENTER_FCN( "BaseOperator::BaseOperator" );

  ;
}

} // end of namespace
