#ifndef FILE_NODEEQN_2004
#define FILE_NODEEQN_2004

#include "baseEQN.hh"

#include "Utils/vector.hh"
#include "DataInOut/WriteInfo.hh"
#include <map>

namespace CoupledField
{

//! a base class for quation data handling

class NodeEQN : public BaseEQN 
{
public:
  
  //! Constructor
  NodeEQN(Grid * aptGrid, 
	  BCs * aptBCs,
	  std::vector<std::string>& asubdoms, 
	  Integer actlevel, 
	  Integer dofsPerNode);
  
  //! Destructor
  virtual ~NodeEQN();
  
  //! Map equation number to position in 
  //! global solution vector
  inline void EQN2SolVectorPos(const Integer eqnNr, Integer &pos) const;

  //! Map vector of equation numbers to 
  //! positions in global solution vector
  virtual void EQN2SolVectorPos(const std::vector<Integer> &eqnNr, 
				std::vector<Integer> &pos) const = 0;
  
  //! Map node number to according equation number(s)
  virtual void Node2EQN(const Integer nodeNr, std::vector<Integer> &eqns) const = 0;
  
  //! Map vector of node numbers to according
  //! vector of equiation numbers
  virtual void Node2EQN(const std::vector<Integer> &nodeNr,
			std::vector<Integer> &eqnNr) const = 0;
 
protected:
  //! Map all element nodes to local numbering
  void CalcMesh2PDENode(std::vector<Integer> & mesh2PDENode,
			std::vector<Integer> & pde2MeshNode);

  //! Mapping Local -> Global node numbering
  std::vector<Integer> pde2MeshNode_;

  //! Mappuing Global -> Local node numbering
  std::vector<Integer> mesh2PDENode_;
  
  //! Mapping for EQN->position in solution vector

  //! This mapping can already be done at 
  //! this level, since the mapping of
  //! an equation to a postition in the 
  //! solution vector is uniquely defined
  std::vector<Integer> eqn2Pos_;

  

};

  // -----------------------------------------------------------------------
  // Inline function definition
  // -----------------------------------------------------------------------
  
  void NodeEQN::EQN2SolVectorPos(const Integer eqnNr, Integer &pos) const
  {
    
    ENTER_IFCN( "NodeEQN::EQN2SolVectorPos(Integer)" );
#ifdef CHECK_INITIALIZED
    if(!isInitialized_)
      Info->Error("Use of uninitialized object. Call 'CalcMapping' before use");
#endif
    
#ifdef CHECK_INDEX
    if (eqnNr > numEqns_)
      Info->Error("Index out of bounds");
#endif
    pos = eqn2Pos_[eqnNr];
  }
  
}  // end of namespace

#endif // FILE_SCALARNODEEQN
