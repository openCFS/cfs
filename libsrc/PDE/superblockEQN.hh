#ifndef FILE_SUPERBLOCKEQN_2004
#define FILE_SUPERBLOCKEQN_2004

#include "nodeEQN.hh"

namespace CoupledField
{

//! Class for numbering equations in a superblock way.
//! This class is hardcoded for use with piezoPDE,
//! since in summer there (hopefully) will be a general
//! mechanism to handel direct-coupled PDEs, so that
//! this class will superfluous.

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
  
  //! 
  Integer GetNumMechEQNs()
  {return numMechEQNs_;}

  //!
  Integer GetNumElecEQNs()
  {return numElecEQNs_;}

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
    void ReorderMapping(Integer *order) {
      Error("ReorderMapping not working for SuperBlockEQN");
    }
  

private:

  //! number of mechanic eqns
  Integer numMechEQNs_;

  //! number of electric eqns
  Integer numElecEQNs_;

  //! Contains the equation numbers
  Matrix<Integer> mechNode2EQN_;
  
  //! Contains the according column numbers
  StdVector<Integer> elecNode2EQN_;

};

  
}  // end of namespace

#endif // FILE_SUPERBLOCKEQN
