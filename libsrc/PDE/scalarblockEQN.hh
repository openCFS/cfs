#ifndef FILE_SCALARBLOCKEQN_2004
#define FILE_SCALARBLOCKEQN_2004

#include "nodeEQN.hh"

namespace CoupledField
{

//! a base class for quation data handling

class ScalarBlockEQN : public NodeEQN 
{
public:
  
  //! Constructor
  ScalarBlockEQN(Grid * aptGrid, 
		 BCs * aptBCs,
		 StdVector<std::string>& asubdoms, 
		 Integer actlevel, 
		 Integer dofsPerNode);
  
  //! Destructor
  virtual ~ScalarBlockEQN();
  
  //! Calculate the mapping after Dirichlet and
  //! constraint nodes were set
  void CalcMapping();
  
  //! Print the mapping nodes <-> EQNs
  void Print(std::ostream & out) const;
  
  // -----------------------------------------------------------------------
  // Functions for mapping node numbers <-> EQN numbers
  // -----------------------------------------------------------------------
  
  //! Map node number and dof to according equation number
  void Node2EQN(const Integer nodeNr, 
		const Integer dof,
		Integer & eqnNr,
		Integer & eqnDof) const;
  
  //! Map node number to according equation number(s)
  void Node2EQN(const Integer nodeNr, StdVector<Integer> &eqns) const;
  
  //! Map vector of node numbers to according
  //! vector of equiation numbers
  void Node2EQN(const StdVector<Integer> &nodeNr,
		StdVector<Integer> &eqnNr) const;

  //! Maps the equation numbers according to the reordering

  //! If the algebraic system has re-ordered the equations, e.g. to improve
  //! the performance of a sparse direct solver, this method must be used to
  //! incorporate the resulting change into the dof to eqn map.
  //! \param order permutation vector containing the re-ordering, i.e.
  //!              order[k-1] is the new equation number for the dof that
  //!              was previously associated to the k-th equation number.
  //! \note It is safe to pass a NULL pointer to the method. In this case
  //!       the method assumes that no re-ordering was performed by the
  //!       algebraic system.
  void ReorderMapping(Integer *order);

private:

  //@{ 
  //! Mapping for Node->EQN. 
  //! For multi-dof PDEs one node
  //! corresponds to variuos EQN numbers,
  //! therefore a compressed row storage format
  //! is chosen. (numPDENodes x dofsPerNode)

  //! Contains the equation numbers
  Matrix<Integer> pdeNode2EQN_;
  //@}
};

  
}  // end of namespace

#endif // FILE_SCALARBLOCKEQN
