#ifndef FILE_BASEOPERATOR_2003
#define FILE_BASEOPERATOR_2003

#include <General/environment.hh>
#include <vector>

namespace CoupledField
{

  class Grid;
  class BasePDE;

 //! Base class for operators working on elements
  /*! Class BaseOperator is base class from which different kinds of operators are derived
   */
class BaseOperator
{
public:

   //! Constructor
  BaseOperator(Grid * ptGrid, 
	       BasePDE * ptPDE,  
	       std::vector<Integer> *ptMesh2PDENode, 
	       const Integer level,
	       Boolean isaxi=FALSE);
  
  //! Destructor
  virtual ~BaseOperator() = 0;

protected:

  Grid * ptGrid_;                       //!< pointer to grid
  BasePDE * ptPDE_;
  std::vector<Integer> *ptMesh2PDENode_; //!< pointer to node transformation array
  Integer level_;                       //!< current level 
  Boolean isaxi_;
};

} // end of namespace

#endif
