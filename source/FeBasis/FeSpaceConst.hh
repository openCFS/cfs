// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
//=================================
/*
 * \file   FeSpaceConst.cc
 * \brief  This FeSpace provides an interface for adding single equations to the system,
 *         e.g. equations of an equivalent circuit of lumped elements.
 *         Basically, any type of entitiy list can be passed to the corresponding
 *         FeFunction. For each unique id of the entries of the list, one equation is added.
 *         This id is obtained from EntityIterator::GetIdString().
 *         This can be dangerous if multiple entity types are added to the space with a
 *         non-unique id. Hence, the entity types are restricted. However, basically any
 *         entities could be connected by equations by this space.
 *
 * \date   May 6, 2014
 * \author dperchto
 */
//=================================

#ifndef FESPACECONST_HH
#define FESPACECONST_HH

#include "FeBasis/FeSpace.hh"
#include <unordered_map>

namespace CoupledField {

class FeSpaceConst : public FeSpace {
public:
  
  //! Constructor
  FeSpaceConst(PtrParamNode paramNode, PtrParamNode infoNode, Grid* ptGrid, bool isAVExc );

  //! Destructor
  virtual ~FeSpaceConst();
  
  //! \copydoc FeSpace::Init
  virtual void Init( shared_ptr<SolStrategy> solStrat );

  //! This space does not have elements.
  virtual BaseFE* GetFe( const EntityIterator ent );
  
  //! This space does not have elements.
  virtual BaseFE* GetFe( const EntityIterator ent ,
                         IntScheme::IntegMethod& method,
                         IntegOrder & order );
  
  //! This space does not have elements.
  virtual BaseFE* GetFe( UInt elemNum );

  //! \copydoc FeSpace::GetNumFunctions
  virtual UInt GetNumFunctions( const EntityIterator ent );

  //! Return equation numbers
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent );

  //! Return equation numbers, does the same as called with 2 arguments
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
                        , UInt dof );

  //! Return equation numbers, does the same as called with 2 arguments
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent,
                        UInt dof, BaseFE::EntityType ); 

  //! Return equation numbers, does the same as called with 2 arguments
  virtual void GetEqns( StdVector<Integer>& eqns, const EntityIterator ent,
                        BaseFE::EntityType ); 

  //! This space does not have elements.
  virtual void GetElemEqns(StdVector<Integer>& eqns,const Elem* elem);

  //! This space does not have elements.
  virtual void GetElemEqns(StdVector<Integer>& eqns,const Elem* elem, UInt dof);

  //! Maps the entity ids to equation numbers
  virtual void Finalize();

  //! This FeSpace does not approximate space.
  virtual void MapCoefFctToSpace(StdVector<shared_ptr<EntityList> > support, 
                                 shared_ptr<CoefFunction> coefFct,
                                 shared_ptr<BaseFeFunction> feFct,
                                 std::map<Integer, Double>& vals,
                                 bool cache,
                                 const std::set<UInt>& comp = std::set<UInt>(),
                                 bool updatedGeo = true);

  //! This FeSpace does not approximate space.
  virtual void MapCoefFctToSpace(StdVector<shared_ptr<EntityList> > support, 
                                 shared_ptr<CoefFunction> coefFct,
                                 shared_ptr<BaseFeFunction> feFct,
                                 std::map<Integer, Complex>& vals,
                                 bool cache,
                                 const std::set<UInt>& comp = std::set<UInt>(),
                                 bool updatedGeo = true);

  //! \copydoc FeSpace::IsSameEntityApproximation
  virtual bool IsSameEntityApproximation( shared_ptr<EntityList> list,
                                          shared_ptr<FeSpace> space );

  virtual void InsertElemsToCoilList(shared_ptr<ElemList> eL, shared_ptr<CoilList> cL);

protected:

  //! This FeSpace does not approximate space.
  virtual void SetRegionElements( RegionIdType region, MappingType mType,
                                  const ApproxOrder& order,
                                  PtrParamNode infoNode );
  
  //! \copydoc FeSpace::CheckConsistency
  virtual void CheckConsistency();

  //! \copydoc FeSpace::SetDefaultIntegration
  virtual void SetDefaultIntegration(PtrParamNode infoNode );

  //! \copydoc FeSpace::SetDefaultElements
  virtual void SetDefaultElements(PtrParamNode infoNode);

  //! does not do anything
  virtual void MapNodalBCs();



private:

  //! can be removed if this space is generalized
  std::set<EntityList::ListType> allowedEntities_;

  bool isAVExc_;

  //! maps every inserted ElemList to a CoilList (only needed for MagEdgeSpecialAVPDE)
  std::unordered_map<UInt, EntityIterator> elemToCoilMap_;

  //! checks the passed entity iterator if it is allowed
  void CheckEntityType(const EntityIterator ent) const;

  //! maps entity id provided by EntityIterator::GetIdString() to equation number
  std::unordered_map<std::string,Integer> equationMap_;

};

}

#endif
