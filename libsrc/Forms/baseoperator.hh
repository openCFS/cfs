#ifndef FILE_BASEOPERATOR_2003
#define FILE_BASEOPERATOR_2003

#include <iostream>
#include <General/environment.hh>
#include <vector>

namespace CoupledField
{

  class Grid;

 //! Base class for operators working on elements
  /*! Class BaseOperator is base class from which different kinds of operators are derived
   */
class BaseOperator
{
public:

   //! Constructor
  BaseOperator(Grid * ptGrid, std::vector<Integer> *ptMesh2PDENode, const Integer level);
  
  //! Destructor
  virtual ~BaseOperator() = 0;

protected:

  Grid * ptGrid_;                       //!< pointer to grid
  std::vector<Integer> *ptMesh2PDENode_; //!< pointer to node transformation array
  Integer level_;                       //!< current level 

};

} // end of namespace

#endif
