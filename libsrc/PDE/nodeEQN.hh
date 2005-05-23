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
            StdVector<RegionIdType>& asubdoms, 
            Integer dofsPerNode);
  
    //! Destructor
    virtual ~NodeEQN();
  
    //! Get number of local nodes
    inline Integer GetNumLocalNodes()
    {return pde2MeshNode_.GetSize();}

    //! Get number of global nodes
    inline Integer GetNumGlobalNodes()
    {return mesh2PDENode_.GetSize();}
  
    //! Get number of local elements
    inline Integer GetNumLocalElems()
    {return pde2MeshElem_.GetSize();}

    //! Get number of global elements
    inline Integer GetNumGlobalElems() 
    {return mesh2PDEElem_.GetSize();}
  
    //! Map node number and dof to according equation number
    virtual void Node2EQN(const Integer nodeNr, 
                          const Integer dof,
                          Integer & eqnNr,
                          Integer & eqnDof) const = 0;
  
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

    //! Map global to local node number
    //! (needed for nodal displacement of grid)
    inline Integer Mesh2PDENode(const Integer meshNode) const;

    //! Map local to global node number
    void PDE2MeshNode(StdVector<Integer> & meshNodes,
                      const StdVector<Integer> & pdeNodes) const;

    //! Map local to global node number
    inline Integer PDE2MeshNode(const Integer pdeNode) const;


    //! Map global to local elem number
    inline Integer Mesh2PDEElem(const Integer elemNumGlob) const;

    //! Map local to global elem number
    inline Integer PDE2MeshElem(const Integer elemNumLoc) const;

  protected:
    //! Calculate mapping local<->global for nodes and elems
    void CalcLocalGlobalMapping(StdVector<Integer> & mesh2PDENode,
                                StdVector<Integer> & pde2MeshNode,
                                StdVector<Integer> & mesh2PDEElem,
                                StdVector<Integer> & pde2MeshElem);

    //! Mapping Local -> Global node numbering
    StdVector<Integer> pde2MeshNode_;

    //! Mapping Global -> Local node numbering
    StdVector<Integer> mesh2PDENode_;
 
    //! Element mapping from local->global
    StdVector<Integer> pde2MeshElem_;
  
    //! Element mapping from global->local
    StdVector<Integer> mesh2PDEElem_;
  
  protected:
  
    //! Default constructor is disallowed
    NodeEQN() {};

  };

  // -----------------------------------------------------------------------
  // Inline function definition
  // -----------------------------------------------------------------------
  
  Integer NodeEQN::Mesh2PDENode(const Integer meshNode) const
  {
    ENTER_FCN( "NodeEQN::Mesh2PDENode" );
#ifdef CHECK_INDEX
    if (meshNode > mesh2PDENode_.GetSize())
      Error( "Index out of bounds", __FILE__, __LINE__ );
#endif
    return mesh2PDENode_[meshNode-1];
  }
  
  //! Map local to global node number
  Integer NodeEQN::PDE2MeshNode(const Integer pdeNode) const {
    ENTER_FCN( "NodeEQN::PDEN2Meshode" );
#ifdef CHECK_INDEX
    if (pdeNode > pde2MeshNode_.GetSize())
      Error( "Index out of bounds", __FILE__, __LINE__ );
#endif
    return pde2MeshNode_[pdeNode-1];
  }



  Integer NodeEQN::Mesh2PDEElem(const Integer elemNumGlob) const
  {
    ENTER_IFCN( "NodeEQN::Mesh2PDEElem" );
#ifdef CHECK_INDEX
    if (elemNumGlob > mesh2PDEElem_.GetSize())
      Error( "Index out of bounds", __FILE__, __LINE__ );
#endif
    return mesh2PDEElem_[elemNumGlob-1];
  }

  Integer NodeEQN::PDE2MeshElem(const Integer elemNumLoc) const
  {
    ENTER_IFCN( "NodeEQN::PDE2MeshElem" );
#ifdef CHECK_INDEX
    if (elemNumLoc > pde2MeshElem_.GetSize())
      Error( "Index out of bounds", __FILE__, __LINE__ );
#endif
    return pde2MeshElem_[elemNumLoc-1];
  }
}  // end of namespace

#endif // FILE_SCALARNODEEQN
