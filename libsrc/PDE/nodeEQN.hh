#ifndef FILE_NODEEQN_2004
#define FILE_NODEEQN_2004

#include "baseEQN.hh"

#include "Utils/vector.hh"
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
	  StdVector<std::string>& asubdoms, 
	  Integer actlevel, 
	  Integer dofsPerNode);
  
  //! Destructor
  virtual ~NodeEQN();
  
  //! Map equation number to position in 
  //! global solution vector
  inline void EQN2SolVectorPos(const Integer eqnNr, Integer &pos) const;

  //! Map vector of equation numbers to 
  //! positions in global solution vector
  virtual void EQN2SolVectorPos(const StdVector<Integer> &eqnNr, 
				StdVector<Integer> &pos) const = 0;
  
  //! Map node number and dof to according equation number
  virtual Integer Node2EQN(const Integer nodeNr, 
			   const Integer dof = 1) const = 0;
  
  //! Map node number to according equation number(s)
  virtual void Node2EQN(const Integer nodeNr, StdVector<Integer> &eqns) const = 0;
  
  //! Map vector of node numbers to according
  //! vector of equiation numbers
  virtual void Node2EQN(const StdVector<Integer> &nodeNr,
			StdVector<Integer> &eqnNr) const = 0;
 
  //! Map global to local node numbers
  //! (needed for nodal displacement of grid)
  void Mesh2PDENode(StdVector<Integer> & PDENodes,
		    const StdVector<Integer> & MeshNodes) const;

  //! /! Map global to local node number
  //! (needed for nodal displacement of grid)
  inline Integer Mesh2PDENode(const Integer meshNode) const;

protected:
  //! Map all element nodes to local numbering
  void CalcMesh2PDENode(StdVector<Integer> & mesh2PDENode,
			StdVector<Integer> & pde2MeshNode);

  //! Mapping Local -> Global node numbering
  StdVector<Integer> pde2MeshNode_;

  //! Mappuing Global -> Local node numbering
  StdVector<Integer> mesh2PDENode_;
  
  //! Mapping for EQN->position in solution vector

  //! This mapping can already be done at 
  //! this level, since the mapping of
  //! an equation to a postition in the 
  //! solution vector is uniquely defined
  StdVector<Integer> eqn2Pos_;

  

};

  // -----------------------------------------------------------------------
  // Inline function definition
  // -----------------------------------------------------------------------
  
  void NodeEQN::EQN2SolVectorPos(const Integer eqnNr, Integer &pos) const
  {
    
    ENTER_IFCN( "NodeEQN::EQN2SolVectorPos(Integer)" );
#ifdef CHECK_INITIALIZED
    if(!isInitialized_)
      Error("Use of uninitialized object. Call 'CalcMapping' before use");
#endif
    
#ifdef CHECK_INDEX
    if (eqnNr -1 > numEqns_)
      Error("Index out of bounds", __FILE__, __LINE__);
#endif
    pos = eqn2Pos_[eqnNr-1];
  }
  
  Integer NodeEQN::Mesh2PDENode(const Integer meshNode) const
  {
    ENTER_FCN( "NodeEQN::Mesh2PDENode" );
#ifdef CHECK_INDEX
    if (meshNode > mesh2PDENode_.GetSize())
      Error( "Index out of bound", __FILE__, __LINE__ );
#endif
    return mesh2PDENode_[meshNode-1];
  }
  
}  // end of namespace

#endif // FILE_SCALARNODEEQN
