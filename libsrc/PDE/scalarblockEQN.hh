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
  
  //! Map vector of equation numbers to 
  //! positions in global solution vector
  void EQN2SolVectorPos(const StdVector<Integer> &eqnNr, 
			StdVector<Integer> &pos) const;
  
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
