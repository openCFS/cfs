// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_DISCONT_EQNMAP_HH
#define FILE_CFS_DISCONT_EQNMAP_HH

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "eqnMap.hh"
//#include "Utils/StdVector.hh"
//#include "Domain/resultInfo.hh"
//#include "Utils/vector.hh"

namespace CoupledField {
  
  //class EntityList;
  //class NodeList;
  //class ElemList;
  //class EntityIterator;
  //struct ResultInfo;
  //class ParamNode; 

  //! Class for mapping entities and continuous ansatz functions to equation numbers
class EntityIterator;
class EntityList;
class Grid;
struct ResultInfo;

  class DiscontinuousEqnMap : public EqnMap {
    
  public:
    
    //! Constructor
    DiscontinuousEqnMap( Grid* ptGrid, PdeIdType, bool usePenalty);

    //! Constructor
    DiscontinuousEqnMap(Grid* grid, PdeIdType pdeId, bool usePenalty,EqnMap* startMap);
    
    //! Destructor
    virtual ~DiscontinuousEqnMap();
    
    //! Maps the equation numbers according to the reordering

    //! If the algebraic system has re-ordered the equations, e.g. to improve
    //! the performance of a sparse direct solver, this method must be used to
    //! incorporate the resulting change into the dof to eqn map.
    //! \param order permutation vector containing the re-ordering, i.e.
    //!              order[k-1] is the new equation number for the dof that
    //!              was previously associated to the k-th equation number.
    //! \note
    //! - It is safe to pass a NULL pointer to the method. In this case
    //!   the method assumes that no re-ordering was performed by the
    //!   algebraic system.
    //! - Since the re-ordering is private property of the equation data
    //!   object, we re-set order to NULL in this method.
    //! - For the sake of a low memory foot-print we delete the permutation
    //!   buffer once the reordering was incorporated
    virtual void ReorderMapping( const StdVector<UInt>& order );
    // ======================================================================
    // EQUATION MAPPING
    // ======================================================================
    
    //@{
    //! \name Equation Mapping
    //! Map result and entities to equation numbers
    virtual void GetEqns( StdVector<Integer>& eqns,
                  const ResultInfo& result, const EntityIterator& it ) const;

    //! Map result, entities and dof numbers to equation numbers
    virtual void GetEqns( StdVector<Integer>& eqns,
                  const ResultInfo& result, const EntityIterator& it,
                  UInt dof ) const;
    
    //! Mpa combination of result, entity and dof to single equation number
    virtual Integer GetEqn( const ResultInfo& result, 
                    const EntityIterator& it, UInt dof ) const;


    //! Name mapping function for obtaining a nodal equation
    virtual Integer GetNodeEqn( const ResultInfo& result, UInt nodeNr, UInt dof );

    //! Name mapping function for obtaining a nodal equation
    virtual Integer GetNodeEqn( UInt nodeNr, UInt dof );

    virtual void GetNodeEqn( const StdVector<UInt>& nodeNrs, 
                     StdVector<Integer>& eqnNrs );

    //@}
    
    // ======================================================================
    // LOCAL/GLOBAL MAPPING OF MESH ENTITIES
    // ======================================================================
    
    //@{ 
    //! \name Local/Global Mapping of Mesh Entities
    //! Map global to local node numbers
    //! (needed for nodal displacement of grid)
    virtual void Mesh2PdeNode(StdVector<UInt> & PdeNodes,
                      const StdVector<UInt> & MeshNodes) const;

    //! Map global to local node number
    //! (needed for nodal displacement of grid)
    virtual UInt Mesh2PdeNode(const UInt meshNode) const;
    
    //@}

    // ======================================================================
    // MISCELLANEOUS
    // ======================================================================
    //@{
    //! \name Miscellaneous

    //! Print the mapping nodes <-> Eqns
    virtual void ToInfo(PtrParamNode in) const;

    //@}
    
  private:

    // ======================================================================
    // PRIVATE HELPER METHODS
    // ======================================================================

    //! Trigger local/global mapping of element/nodal numbers
    virtual void CalcNodeElemMapping();

    //! Trigger local/global mapping of edge numbers
    virtual void CalcEdgeMapping();

    //! Trigger local/global mapping of face numbers
    virtual void CalcFaceMapping();

    //! Calculate equation numbers for nodes

    //! Triggers the mapping for nodal equations
    //! \param phase Parameter indicating if numbering of free equations (1)
    //!              has to be performed or if the fixed equations have to
    //!              be mapped (2)
    virtual void CalcNodalEquations( UInt phase );

    //! Calculate equation numbers for elements interior (pfem)

    //! Triggers the mapping for element equations
    //! \param phase Parameter indicating if numbering of free equations (1)
    //!              has to be performed or if the fixed equations have to
    //!              be mapped (2)
    virtual void CalcElemInteriorEquations( UInt phase );

    //! Calculate equation numbers for elements (constants)

    //! Triggers the mapping for element equations
    //! \param phase Parameter indicating if numbering of free equations (1)
    //!              has to be performed or if the fixed equations have to
    //!              be mapped (2)
    virtual void CalcElemConstEquations( UInt phase );

    //! Calculate equation numbers for edges

    //! Triggers the mapping for edge equations
    //! \param phase Parameter indicating if numbering of free equations (1)
    //!              has to be performed or if the fixed equations have to
    //!              be mapped (2)
    virtual void CalcEdgeEquations( UInt phase );

    //! Calculate equation numbers for faces

    //! Triggers the mapping for face equations
    //! \param phase Parameter indicating if numbering of free equations (1)
    //!              has to be performed or if the fixed equations have to
    //!              be mapped (2)
    virtual void CalcFaceEquations( UInt phase );

    //! Extract node numbers from given entitylist
    virtual void GetNodesOfEntities( StdVector<UInt>& nodeNr, 
                             shared_ptr<EntityList> ent );
    
    // ======================================================================
    // GLOBAL/LOCAL MAPPING
    // ======================================================================

    //! Mapping Global -> Local node numbering for discontinuous elements
    StdVector< StdVector< Integer* > > mesh2DisPdeNode_;
 
    //! Edge mapping from global->local discontiuous elements
    StdVector< StdVector< Integer* > > mesh2DisPdeEdge_;

    //! Edge mapping from global->local discontiuous elements
    StdVector< StdVector< Integer* > > mesh2DisPdeFace_;
  };
  
  
#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class EqnMap
  //! 
  //! \purpose
  //!
  //! This class is the base class which defines the interface of all
  //! objects that deal with the assignment of equation numbers to degrees
  //! of freedom. The following gives an overview on the meaning of the
  //! different equation numbers
  //! <center><img src="../AddDoc/splitting.png"></center>
  //!
  //! \collab 
  //!
  //! \implement 
  //!
  //! \status In use
  //!
  //! \unused 
  //!
  //! \improve
  //! 

#endif
}

#endif
