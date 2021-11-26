// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
//================================================================================================
/*!
 *       \file     FeSpaceNodal.hh 
 *       \brief    This space implements basic functionality shared by all spaces which rely on 
 *                 nodal elements of arbitrary order (i.e. elements with delta property). 
 *                 Currently just Lagrangian types implemented.
 *
 *       \date     01/01/2013
 *       \author   Andreas Hueppe
 */
//================================================================================================

#ifndef FESPACENODAL_HH 
#define FESPACENODAL_HH

#include "FeSpace.hh"


namespace CoupledField {

//! Space for nodal finite elements
class FeSpaceNodal : public FeSpace {
public:

  //! Constructor
  FeSpaceNodal (PtrParamNode paramNode, PtrParamNode infoNode, Grid* ptGrid );

  //! Destructor
  virtual ~FeSpaceNodal();

  //! \copydoc FeSpace::MapCoefFctToSpace
  virtual void MapCoefFctToSpace(StdVector<shared_ptr<EntityList> > support,
                                shared_ptr<CoefFunction> coefFct,
                                shared_ptr<BaseFeFunction> feFct,
                                std::map<Integer, Double>& vals,
                                bool cache,
                                const std::set<UInt>& comp = std::set<UInt>(),
                                bool updatedGeo = true){
    MapCoefFctToSpacePriv<Double>(support,coefFct,feFct,vals,cache,comp,updatedGeo);
  }

  //! \copydoc FeSpace::MapCoefFctToSpace
  virtual void MapCoefFctToSpace(StdVector<shared_ptr<EntityList> > support,
                                 shared_ptr<CoefFunction> coefFct,
                                 shared_ptr<BaseFeFunction> feFct,
                                 std::map<Integer, Complex>& vals,
                                 bool cache,
                                 const std::set<UInt>& comp = std::set<UInt>(),
                                 bool updatedGeo = true){
    MapCoefFctToSpacePriv<Complex>(support,coefFct,feFct,vals,cache,comp,updatedGeo);
  }

protected:

  ///In a nodal space we can gather the coordinates of the DOFs explicitly
  void GetCoordinateRepresentation( shared_ptr<EntityList>   eList,
                                        StdVector< StdVector<UInt> > & idxMap,
                                        StdVector< Vector<Double> > & coords,
                                        bool updatedGeo);

  ///Get actual coordinates
  void GetNodalCoords(StdVector<Vector<Double> > & coords,
                        StdVector<UInt> & vNodeNums,
                        const shared_ptr<EntityList> & entities,
                        bool updatedGeo);


  //! get a nodal equation number
  StdVector<Integer> GetVNodeEqns(UInt nodeNr){
    return this->nodeMap_[nodeNr];
  }
private:

  template<typename T>
  void MapCoefFctToSpacePriv(StdVector<shared_ptr<EntityList> > support,
                             shared_ptr<CoefFunction> coefFct,
                             shared_ptr<BaseFeFunction> feFct,
                             std::map<Integer, T>& vals,
                             bool cache,
                             const std::set<UInt>& comp = std::set<UInt>(),
                             bool updatedGeo = true);

};
}
#endif
