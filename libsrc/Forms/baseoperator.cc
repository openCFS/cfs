#include "baseoperator.hh"
#include "Domain/grid.hh"
namespace CoupledField
{

BaseOperator::BaseOperator(Grid * ptGrid, BasePDE * ptPDE, 
			   std::vector<Integer> *ptMesh2PDENode, 
			   const Integer level, 
			   Boolean isaxi)
  : isaxi_(isaxi)
{
#ifdef TRACE
  (*trace) << "entering BaseOperator::BaseOperator" << std::endl;
#endif

  this->ptGrid_ = ptGrid;
  this->ptPDE_ = ptPDE;
  this->level_ = level;
  this->ptMesh2PDENode_ = ptMesh2PDENode;
}


BaseOperator::~BaseOperator()
{
#ifdef TRACE
  (*trace) << "entering BaseOperator::BaseOperator" << std::endl;
#endif
  ;
}

} // end of namespace
