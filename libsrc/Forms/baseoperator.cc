#include "baseoperator.hh"
#include "Domain/grid.hh"
namespace CoupledField
{

BaseOperator::BaseOperator(Grid * ptGrid, std::vector<Integer> *ptMesh2PDENode, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering BaseOperator::BaseOperator" << std::endl;
#endif
  this->ptGrid_ = ptGrid;
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
