// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     EntityAssociation.cc
 *       \brief    <Description>
 *
 *       \date     Jan 6, 2016
 *       \author   ahueppe
 */
//================================================================================================

#include "EntityAssociation.hh"
#include "boost/scoped_array.hpp"

namespace CoupledField{

#ifdef USE_CGAL

  // ========================================================================
  //  C G A L  -  S P E C I F I C   I M P L E M E N T A T I O N
  // ========================================================================


  //! Define box handler, which additionally stores an index
  typedef CGAL::Box_intersection_d
      ::Box_with_handle_d<double,3,const UInt*> HandleBox;


  typedef CGAL::Bbox_3 BBox3D;

  // Iterator reporter class, returning the two ids of the CGAL-Boxes,
  // the first being the element number, the second being the node index
  template <class OutputIterator>
  struct CGAL_ElemElemIdReporter {
    OutputIterator it;
    CGAL_ElemElemIdReporter(OutputIterator i  )
    : it(i) {} // store iterator in object

    // We write the id-number of box a to the output iterator assuming
    // that box b (the query box) is not interesting in the result.
    void operator()( const HandleBox& a, const HandleBox& b) {
      UInt aNum = *a.handle();
      UInt bNum = *b.handle();
      std::pair<UInt, UInt > pair;
      pair.first = aNum;
      pair.second = bNum;
      *it++ = pair;
    }
  };
  // helper function to create the function object
  template <class Iter>
  CGAL_ElemElemIdReporter<Iter> elemElemIdReporter(Iter it)
  { return CGAL_ElemElemIdReporter<Iter>(it); }


std::vector< std::pair<UInt, UInt > > BoundingBoxAssociate::AssociateEntities(const StdVector<Elem*>& entList1,
                                                                       const StdVector<Elem*>& entList2,
                                                                       Grid* g1, Grid* g2,
                                                                       Double tolerance){
  std::vector< std::pair<UInt, UInt > > assocs;

  std::array<Double,6> bbox;
  std::vector<HandleBox> boxes1(entList1.GetSize());
  std::vector<HandleBox> boxes2(entList2.GetSize());

  //create both bounding box lists
  for(UInt aEl=0;aEl<entList1.GetSize();++aEl){
    Elem* aElem = entList1[aEl];
    g1->CreateBBoxFromElement(aElem, tolerance, &bbox[0]);
    HandleBox hbox(BBox3D(bbox[0], bbox[1], bbox[2],
                          bbox[3], bbox[4], bbox[5]),
                          &aElem->elemNum);
    boxes1[aEl] = hbox;
  }

  //create both bounding box lists
  for(UInt aEl=0;aEl<entList2.GetSize();++aEl){
    Elem* aElem = entList2[aEl];
    g2->CreateBBoxFromElement(aElem, tolerance, &bbox[0]);
    HandleBox hbox(BBox3D(bbox[0], bbox[1], bbox[2],
                          bbox[3], bbox[4], bbox[5]),
                          &aElem->elemNum);
    boxes2[aEl] = hbox;
  }

  CGAL::box_intersection_d( boxes1.begin(), boxes1.end(),
                            boxes2.begin(), boxes2.end(),
                            elemElemIdReporter( std::back_inserter( assocs )));

  return assocs;
}

std::vector< std::pair<UInt, UInt > > BoundingBoxAssociate::AssociateEntities(const StdVector< Vector<Double> >& points,
                                                                       const StdVector<Elem*>& elemList,
                                                                       Grid* elGrid, Double tolerance){
  std::vector< std::pair<UInt, UInt > > assocs;

  std::array<Double,6> bbox;
  std::vector<HandleBox> elBoxes(elemList.GetSize());
  std::vector<HandleBox> pointBoxes(points.GetSize());

  //create both bounding box lists
  for(UInt aEl=0;aEl<elemList.GetSize();++aEl){
    Elem* aElem = elemList[aEl];
    elGrid->CreateBBoxFromElement(aElem, tolerance, &bbox[0]);
    HandleBox hbox(BBox3D(bbox[0], bbox[1], bbox[2],
                          bbox[3], bbox[4], bbox[5]),
                          &aElem->elemNum);
    elBoxes[aEl] = hbox;
  }

  // create also temporary index array (will be automatically deleted)
  boost::scoped_array<UInt> nodeIndices(new UInt[points.GetSize()]);

  for( UInt i = 0; i < points.GetSize(); ++i ) {
    nodeIndices[i] = i;
    pointBoxes[i] = BoundingBoxAssociate::CreateBoxFromCoord( points[i],
                                                      &nodeIndices[i], tolerance );
  }

  CGAL::box_intersection_d( elBoxes.begin(), elBoxes.end(),
                            pointBoxes.begin(), pointBoxes.end(),
                            elemElemIdReporter( std::back_inserter( assocs )));

  return assocs;
}

HandleBox BoundingBoxAssociate::CreateBoxFromCoord( const Vector<double>& coords, UInt* id,
                                            Double tol )
{
  if(coords.GetSize()==2){
    return HandleBox(BBox3D(coords[0]-tol/2.0, coords[1]-tol/2.0, 0.0,
                            coords[0]+tol/2.0, coords[1]+tol/2.0, 0.0), id);
  }else{
    return HandleBox(BBox3D(coords[0]-tol/2.0, coords[1]-tol/2.0,
                            coords[2]-tol/2.0, coords[0]+tol/2.0,
                            coords[1]+tol/2.0, coords[2]+tol/2.0), id);
  }
}

#elif USE_LIBFBI

std::vector< std::pair<UInt, UInt > > BoundingBoxAssociate::AssociateEntities(const StdVector<Elem*>& entList1,
                                                                       const StdVector<Elem*>& entList2,
                                                                       Grid* g1, Grid* g2,
                                                                       Double tolerance){
  std::vector< std::pair<UInt, UInt > > assocs;

  EXCEPTION("Activate CGAL in cmake to use this feature");
  return assocs;
}

std::vector< std::pair<UInt, UInt > > BoundingBoxAssociate::AssociateEntities(const StdVector< Vector<Double> >& points,
                                                                       const StdVector<Elem*>& elemList,
                                                                       Grid* elGrid, Double tolerance){
  std::vector< std::pair<UInt, UInt > > assocs;

  EXCEPTION("Activate CGAL in cmake to use this feature");
  return assocs;
}
#else
std::vector< std::pair<UInt, UInt > > BoundingBoxAssociate::AssociateEntities(const StdVector<Elem*>& entList1,
                                                                       const StdVector<Elem*>& entList2,
                                                                       Grid* g1, Grid* g2,
                                                                       Double tolerance){
  std::vector< std::pair<UInt, UInt > > assocs;

  EXCEPTION("Activate CGAL in cmake to use this feature");
  return assocs;
}

std::vector< std::pair<UInt, UInt > > BoundingBoxAssociate::AssociateEntities(const StdVector< Vector<Double> >& points,
                                                                       const StdVector<Elem*>& elemList,
                                                                       Grid* elGrid, Double tolerance){
  std::vector< std::pair<UInt, UInt > > assocs;

  EXCEPTION("Activate CGAL in cmake to use this feature");
  return assocs;
}
#endif

}
