#ifndef FILE_BLOCKNODEEQN_2004
#define FILE_BLOCKNODEEQN_2004

#include "nodeEQN.hh"

namespace CoupledField
{

//! Equation handling class for PDEs with number
//! of degrees of freedom > 1

class BlockNodeEQN : public NodeEQN
{
public:
  
  //! Constructor
  BlockNodeEQN(Grid * aptgrid, 
	       BCs * aptBCs,
	       StdVector<std::string>& asubdoms, 
	       Integer actlevel, 
	       Integer dofsPerNode);
  
  //! Destructor
  virtual ~BlockNodeEQN();
  
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

  //! Mapping global node number->EQN
  StdVector<Integer> pdeNode2EQN_;
  
};
  
  
}  // end of namespace

#endif // FILE_BLOCKNODEEQN
