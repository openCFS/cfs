#ifndef FILE_SCALARNODEEQN_2004
#define FILE_SCALARNODEEQN_2004

#include "nodeEQN.hh"

namespace CoupledField
{

//! Euation handling class for scalar valued
//! PDEs.

class ScalarNodeEQN : public NodeEQN
{
public:
  
  //! Constructor
  ScalarNodeEQN(Grid * aptGrid, 
		BCs * aptBCs,
		std::vector<std::string>& asubdoms, 
		Integer actlevel, 
		Integer dofsPerNode);
  
  //! Destructor
  virtual ~ScalarNodeEQN();
  
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
  void EQN2SolVectorPos(const std::vector<Integer> &eqnNr, 
			std::vector<Integer> &pos) const;
  
  //! Map node number to according equation number(s)
  void Node2EQN(const Integer nodeNr, std::vector<Integer> &eqns) const;
  
  //! Map vector of node numbers to according
  //! vector of equiation numbers
  void Node2EQN(const std::vector<Integer> &nodeNr,
		std::vector<Integer> &eqnNr) const;

private:

  //! Mapping global node number->EQN
  std::vector<Integer> pdeNode2EQN_;
  
  
};

  
}  // end of namespace

#endif // FILE_SCALARNODEEQN
