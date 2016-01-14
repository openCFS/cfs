// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     TriaIntersect.hh
 *       \brief    Classes for intersetion of volume elements based on
 *                 triangulation
 *
 *       \date     Dec 11, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_TRIAINTERSECT_HH_
#define SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_TRIAINTERSECT_HH_


#include "ElemIntersect.hh"
#include "Utils/StdVector.hh"
#include "FeBasis/BaseFE.hh"
#include "CoordTetra.hh"

namespace CoupledField{

//! class prociding a intersection based on triagnulation of input elements

//! each elements is converted to a set of tetrahedrons
//! for each set of tetrahedrons, an intersection is performed
class TriaIntersect : public ElemIntersect{

public:

  //! \copydoc ElemIntersect::TriaIntersect
  TriaIntersect(Grid* trgGrid, Grid* srcGrid)
   : ElemIntersect(trgGrid,srcGrid){
    tElem_ = NULL;
  }

  //! copy constructor
  TriaIntersect(const TriaIntersect& inter)
   : ElemIntersect(inter){

    tElem_ = inter.tElem_;
    tTets_ = inter.tTets_;
    intersectingTets_ = inter.intersectingTets_;
    lastTetIdx_ =  inter.lastTetIdx_;
  }

  virtual ~TriaIntersect(){

  }

  //! deep pointer copy
  virtual TriaIntersect* Clone(){
     return new TriaIntersect(*this);
  }

  //! \copydoc ElemIntersect::SetTElem
  virtual void SetTElem( UInt tNum );

  //! \copydoc ElemIntersect::Intersect
  virtual bool Intersect(UInt sNum);

  //! \copydoc ElemIntersect::GetIntersectionElems
  virtual void GetIntersectionElems(StdVector<IntersectionElem*>& interElems);

  //! \copydoc ElemIntersect::GetVolumeAndCenters
  virtual void GetVolumeAndCenters(StdVector<VolCenterInfo>& infos);

private:

  //! Sprits a given tetrahedron by the active hyperplane indicated by tetIdx
  inline void SplitAndDecompose(UInt tetIdx, UInt planeIdx, CoordTetra& tetra,
                                StdVector<CoordTetra>& genTets);

  //! Triangulates a given element according to the baseFE
  //!\param(in) newTElem element to be triangulated
  //!\param(in) aGrid grid pointer for given element
  //!\return vector of Tetrahedrons in coordinate format
  inline StdVector<CoordTetra> GetTetsFromElem(const Elem* newTElem, Grid* aGrid);

  //! pointer to currenty active base element
  const Elem* tElem_;

  //! coordinateTetra corresponding to tElem_
  StdVector<CoordTetra> tTets_;

  //! coordinateTetra from last intersection
  StdVector<CoordTetra> intersectingTets_;

  //! cached triangulated baseFE
  StdVector< StdVector<UInt> > lastTetIdx_;
};

}

#endif /* SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_TRIAINTERSECT_HH_ */
