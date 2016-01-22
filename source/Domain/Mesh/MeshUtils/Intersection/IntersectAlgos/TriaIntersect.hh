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
#include "Utils/ThreadLocalStorage.hh"

namespace CoupledField{

//! class prociding a intersection based on triagnulation of input elements

//! each elements is converted to a set of tetrahedrons
//! for each set of tetrahedrons, an intersection is performed
class TriaIntersect : public ElemIntersect{

public:

  //! \copydoc ElemIntersect::TriaIntersect
  TriaIntersect(Grid* trgGrid, Grid* srcGrid)
   : ElemIntersect(trgGrid,srcGrid){
    tElemNum_ = 0;
    sElemNum_ = 0;
    InitElemMap();
  }

  //! copy constructor
  TriaIntersect(const TriaIntersect& inter)
   : ElemIntersect(inter){

    tElemNum_ = inter.tElemNum_;
    tTets_ = inter.tTets_;
    sElemNum_ = inter.sElemNum_;
    intersectingTets_ = inter.intersectingTets_;
    lastTetIdx_ =  inter.lastTetIdx_;
    InitElemMap();
  }

  //! assignment
  TriaIntersect & operator=(const TriaIntersect& inter){
    this->sGrid_ = inter.sGrid_;
    this->tGrid_ = inter.tGrid_;
    tElemNum_ = inter.tElemNum_;
    tTets_ = inter.tTets_;
    sElemNum_ = inter.sElemNum_;
    intersectingTets_ = inter.intersectingTets_;
    lastTetIdx_ =  inter.lastTetIdx_;
    InitElemMap();
    return *this;
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

  //!exports a list of Tetrahedrons, each as an off file
  //! given file will be overwritten
  //!\param(in) tetList List of tetrahedrons
  //!\param(in) baseFName base name of generated files
  void ExportTetras(StdVector<CoordTetra> tetList,std::string baseFName);


  //! initialize reference element map for triangulation
  void InitElemMap();

  //! map of reference elements for access to triangulation
  std::map<Elem::FEType, FeH1LagrangeExpl* > refFeMap;

  //! coordinateTetra corresponding to tElem_
  StdVector<CoordTetra> tTets_;

  //! storage of target element number
  UInt tElemNum_;

  //! storage of source element number
  UInt sElemNum_;

  //! coordinateTetra from last intersection
  StdVector<CoordTetra> intersectingTets_;

  //! cached, triangulated baseFE
  StdVector< StdVector<UInt> > lastTetIdx_;
};

}

#endif /* SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_TRIAINTERSECT_HH_ */
