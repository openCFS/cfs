// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     FeSpaceL2.hh
 *       \brief    <Description>
 *
 *       \date     Feb 2, 2012
 *       \author   ahueppe
 */
//================================================================================================

#ifndef FESPACEL2_HH_
#define FESPACEL2_HH_

#include "FeBasis/FeSpaceNodal.hh"
#include <unordered_map>

namespace CoupledField {

/*!  \class FeSpaceL2
 *! \brief Class for the disconitnuous space L2
 *     @author A. Hueppe
 *     @date 02/2012
 *
 *    Perhaps in the near future we should recombine the functionaliy with
 *    the FeSpaceH1 as its basically the same
 *
 */
class FeSpaceL2 : public FeSpaceNodal {
public:

  //! Struct containing all virtual nodes for on entity type (Vertex, Face
  struct EntityTypeNodes {
    //! Nodes for all numbers of entity (edgeNodes, faceNodes, innerNodes)
    StdVector<UInt> vNodes;
    //! Offset to vNodes array
    StdVector<UInt> offset;
  };
  
  //! Enum which stores the (Virtual) Nodes of an element according to
    //! their definition on vertices,edges,faces and interior
  typedef std::map< BaseFE::EntityType , EntityTypeNodes> ElemVirtualNodes;
  
  //! Constructor
  FeSpaceL2(PtrParamNode aNode, PtrParamNode infoNode, Grid* ptGrid );

  //! Destructor
  virtual ~FeSpaceL2();
  
  //! Returns the number of (vectorial) unknowns on the element
  virtual UInt GetNumFunctions( const EntityIterator ent );

  //! Return equation numbers
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent );

  //! Return equation numbers for a specific dof
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
                        , UInt dof );
  //! Return equation numbers for a specific dof and entitytype
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
                        , UInt dof, BaseFE::EntityType );

  //! get a nodal equation number
  virtual UInt GetNodeEqn(UInt nodeNr, UInt dof);


  //! Get Equation numbers for a specific element
  virtual void GetElemEqns(StdVector<Integer>& eqns,const Elem* elem);

  //! Get Equation numbers for a specific element
  virtual void GetElemEqns(StdVector<Integer>& eqns,const Elem* elem, UInt dof);


  //! \copydoc FeSpace::GetOlasMappings
  virtual void GetOlasMappings(  StdVector<AlgebraicSys::SBMBlockDef>& sbmBlocks,
                                std::map<UInt,StdVector<std::set<Integer> > >&
                                minorBlocks );

  //! \copydoc FeSpace::GetInteriorSurfaceElems
  virtual void GetInteriorSurfaceElems(RegionIdType region,
                                    shared_ptr<EntityList> & surfElems,
                                    shared_ptr<EntityList> & opposingElems);

  //! \copydoc FeSpace::GetInteriorSurfaceElems
  virtual void GetInteriorSurfaceElems(RegionIdType region,
                                       shared_ptr<EntityList> & surfElems);

  //! \copydoc FeSpace::GetExteriorSurfaceElemsOfFeSpace
  virtual void GetExteriorSurfaceElems(RegionIdType region, shared_ptr<EntityList>   & surfElems );

  //! Add result
  virtual void AddFeFunction( shared_ptr<BaseFeFunction> fct );


  /** @see FeSpace::PrintEqnMap() */
  virtual void PrintEqnMap(std::ostream* file = NULL);

  
  // =====================================================
  //  Nodal section
  // =====================================================

  //! Get all Node Numbers of a specific element according to the mapping GRID based or
  //! POLYNOMIAL Based according to the Type requested
  virtual void GetNodesOfElement( StdVector<UInt>& nodes,
                                  const Elem* ptElem,
                                  BaseFE::EntityType entType = BaseFE::ALL);

  virtual void CreatePolynomialNodes();
  
 protected:

  //! Map Nodal BC Equation NUmbers
  virtual void MapNodalBCs();

  //! Map Nodal Equation Numbers
  virtual void MapNodalEqns(UInt phase);

  //! fill the requested EntityLists and push them to fefunction
  virtual void CreateSurfaceElems();

  // ====================================================================
  // PROCESS USER INPUT
  // ====================================================================

  //! Nodal Equation Map
  //SingleEqnMap nodeMap_;

  //! This is the virtual node Map for standard element it just contains
  //! the connectivity of the element, for higher order elements it contains also 
  //! the virtual node numbers in the correct ordering
  //! This Variable could be extended to store also the coordinates of all nodes
  //! created
  std::unordered_map< UInt, ElemVirtualNodes > virtualNodes_;
  
  // ====================================================================
  // Store surface elements
  // ====================================================================

  //!map that stores interor surface elements
  std::map< RegionIdType, shared_ptr<NcSurfElemList> > interiorElemMap_;

  //!map that stores interor surface elements opposing elements
  std::map< RegionIdType, shared_ptr<NcSurfElemList> > interiorElemMapOpposite_;

  //! map storing exterior elements of regions excluding those who have a
  //! neighbor in another region!
  std::map< RegionIdType, shared_ptr<NcSurfElemList> > exteriorElements_;
};
}


#endif /* FESPACEL2_HH_ */
