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
            UInt dofsPerNode);
  
    //! Destructor
    virtual ~NodeEQN();
  
    //! Get number of local nodes
    inline UInt GetNumLocalNodes()
    {return pde2MeshNode_.GetSize();}

    //! Get number of global nodes
    inline UInt GetNumGlobalNodes()
    {return mesh2PDENode_.GetSize();}
  
    //! Get number of local elements
    inline UInt GetNumLocalElems()
    {return pde2MeshElem_.GetSize();}

    //! Get number of global elements
    inline UInt GetNumGlobalElems() 
    {return mesh2PDEElem_.GetSize();}
  
    //! Map node number and dof to according equation number
    virtual void Node2EQN(const UInt nodeNr, 
                          const UInt dof,
                          Integer & eqnNr,
                          UInt & eqnDof) const = 0;
  
    //! Map node number to according equation number(s)
    virtual void Node2EQN(const UInt nodeNr, StdVector<Integer> &eqns) const = 0;
  
    //! Map vector of node numbers to according
    //! vector of equiation numbers
    virtual void Node2EQN(const StdVector<UInt> &nodeNr,
                          StdVector<Integer> &eqnNr) const = 0;
 
    //! Map global to local node numbers
    //! (needed for nodal displacement of grid)
    void Mesh2PDENode(StdVector<UInt> & PDENodes,
                      const StdVector<UInt> & MeshNodes) const;

    //! Map global to local node number
    //! (needed for nodal displacement of grid)
    inline UInt Mesh2PDENode(const UInt meshNode) const;

    //! Map local to global node number
    void PDE2MeshNode(StdVector<UInt> & meshNodes,
                      const StdVector<UInt> & pdeNodes) const;

    //! Map local to global node number
    inline UInt PDE2MeshNode(const UInt pdeNode) const;


    //! Map global to local elem number
    inline UInt Mesh2PDEElem(const UInt elemNumGlob) const;

    //! Map local to global elem number
    inline UInt PDE2MeshElem(const UInt elemNumLoc) const;

    //! Print the mapping nodes <-> EQNs
    virtual void Print( std::ostream &out ) const = 0;

  protected:
    //! Calculate mapping local<->global for nodes and elems
    void CalcLocalGlobalMapping(StdVector<Integer> & mesh2PDENode,
                                StdVector<UInt> & pde2MeshNode,
                                StdVector<Integer> & mesh2PDEElem,
                                StdVector<UInt> & pde2MeshElem);

    //! Mapping Local -> Global node numbering
    StdVector<UInt> pde2MeshNode_;

    //! Mapping Global -> Local node numbering
    StdVector<Integer> mesh2PDENode_;
 
    //! Element mapping from local->global
    StdVector<UInt> pde2MeshElem_;
  
    //! Element mapping from global->local
    StdVector<Integer> mesh2PDEElem_;
  
  protected:
  
    //! Default constructor is disallowed
    NodeEQN() {};

  };

  // -----------------------------------------------------------------------
  // Inline function definition
  // -----------------------------------------------------------------------
  
  UInt NodeEQN::Mesh2PDENode(const UInt meshNode) const
  {
    ENTER_FCN( "NodeEQN::Mesh2PDENode" );
#ifdef CHECK_INDEX
    if (meshNode > mesh2PDENode_.GetSize())
      Error( "Index out of bounds", __FILE__, __LINE__ );
#endif
    if ( mesh2PDENode_[meshNode-1] < 0 ) {
      (*error) << "MeshNode Nr. " << meshNode 
               << " has no local node number!";
      Error( __FILE__, __LINE__);
    }

    return abs(mesh2PDENode_[meshNode-1]);
  }
  
  //! Map local to global node number
  UInt NodeEQN::PDE2MeshNode(const UInt pdeNode) const {
    ENTER_FCN( "NodeEQN::PDEN2Meshode" );
#ifdef CHECK_INDEX
    if (pdeNode > pde2MeshNode_.GetSize())
      Error( "Index out of bounds", __FILE__, __LINE__ );
#endif
    
    return pde2MeshNode_[pdeNode-1];
  }



  UInt NodeEQN::Mesh2PDEElem(const UInt elemNumGlob) const
  {
    ENTER_IFCN( "NodeEQN::Mesh2PDEElem" );
#ifdef CHECK_INDEX
    if (elemNumGlob > mesh2PDEElem_.GetSize())
      Error( "Index out of bounds", __FILE__, __LINE__ );
#endif
    if ( mesh2PDEElem_[elemNumGlob-1] < 0 ) {
      (*error) << "MeshElem Nr. " << elemNumGlob 
               << " has no local elem number!";
      Error( __FILE__, __LINE__);
    }
    
    return abs(mesh2PDEElem_[elemNumGlob-1]);
  }

  UInt NodeEQN::PDE2MeshElem(const UInt elemNumLoc) const
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
