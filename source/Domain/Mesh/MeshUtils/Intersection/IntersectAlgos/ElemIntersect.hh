// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     ElementIntersection.hh
 *       \brief    <Description>
 *
 *       \date     Dec 18, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_ELEMINTERSECT_HH_
#define SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_ELEMINTERSECT_HH_

#include "Domain/ElemMapping/Elem.hh"
#include "Utils/ThreadLocalStorage.hh"
#include "Domain/Mesh/Grid.hh"

namespace CoupledField{

///Base Interface for volume element intersection algorithms

//! As the creation of full intersection elements is usually costly
//! this class provides also a little cheaper access to the intersection
//! volumes and centroids
class ElemIntersect : public CfsCopyable{

public:

  //! struct fo basic intersection information
  struct VolCenterInfo{
    Double volume;
    Vector<Double> center;
    UInt oElem1;
    UInt oElem2;
  };

  //!constructor
  //! \param(in) g1 the base grid. SetTElem function refers to this
  //! \param(in) g2 operating grid element passed to Intersect method refers to this
  ElemIntersect(Grid* g1, Grid* g2) : g1_(g1),g2_(g2){

  }

  //! copy constructor
  ElemIntersect(const ElemIntersect& inter){
    g1_ = inter.g1_;
    g2_ = inter.g2_;
  }

  //! deep pointer copy
  virtual ElemIntersect* Clone()=0;

  //! destructor
  virtual ~ElemIntersect(){

  }

  //! sets currently active target element
  //! \param(in) tNum new target element number
  virtual void SetTElem(UInt tNum)=0;

  //! intersects currently active target element with given element index
  //! \param(in) sNum element of g2 to intersect
  virtual bool Intersect(UInt sNum)=0;

  //! returns all intersection elements of last call to intersect
  //! \param(out) interElems list of intersection elements
  virtual void GetIntersectionElems(StdVector<IntersectionElem*>& interElems)=0;

  //! returns only volumes and centroids of call to intersect
  //! \param(out) infos list of intersection infos
  virtual void GetVolumeAndCenters(StdVector<VolCenterInfo>& infos)=0;


protected:
  Grid* g1_;
  Grid* g2_;
};

}

#endif /* SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_ELEMINTERSECT_HH_ */
