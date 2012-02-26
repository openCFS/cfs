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

#include "FeBasis/H1/FeSpaceH1.hh"

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
class FeSpaceL2 : public FeSpaceH1 {
public:

  //! Constructor
  FeSpaceL2(PtrParamNode aNode, PtrParamNode infoNode);

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


    //! Dump information to screen
  virtual void PrintEqnMap();


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

  // ====================================================================
  // Store surface elements
  // ====================================================================

  //!map that stores interor surface elements
  std::map< RegionIdType, shared_ptr<NcElemList> > interiorElemMap_;

  //!map that stores interor surface elements opposing elements
  std::map< RegionIdType, shared_ptr<NcElemList> > interiorElemMapOpposite_;

  //! map storing exterior elements of regions excluding those who have a
  //! neighbor in another region!
  std::map< RegionIdType, shared_ptr<NcElemList> > exteriorElements_;
};
}


#endif /* FESPACEL2_HH_ */
