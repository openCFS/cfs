#ifndef FILE_SUPERBLOCKEQN_2004
#define FILE_SUPERBLOCKEQN_2004

#include "nodeEQN.hh"

namespace CoupledField
{

//! a base class for quation data handling

class SuperBlockEQN : public NodeEQN 
{
public:
  
  //! Constructor
  SuperBlockEQN(Grid * aptGrid, 
		BCs * aptBCs,
		StdVector<std::string>& asubdoms, 
		Integer actlevel, 
		Integer dofsPerNode);
  
  virtual ~SuperBlockEQN();
  
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
  Integer Node2EQN(const Integer nodeNr, 
		   const Integer dof) const;
  
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
  StdVector<Integer> node2EQNVal_;
  
  //! Contains the according column numbers
  StdVector<Integer> node2EQNcol_;

  //! Contains the starting positions of each row
  StdVector<Integer> node2EQNrow_;
  //@}
};

  
}  // end of namespace

#endif // FILE_SUPERBLOCKEQN
