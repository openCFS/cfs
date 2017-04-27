// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     EntityAssociation.hh
 *       \brief    <Description>
 *
 *       \date     Jan 6, 2016
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SOURCE_DOMAIN_MESH_MESHUTILS_ENTITYASSOCIATION_HH_
#define SOURCE_DOMAIN_MESH_MESHUTILS_ENTITYASSOCIATION_HH_


#include <def_use_cgal.hh>
#include <def_use_libfbi.hh>

#include "Domain/Mesh/Grid.hh"

#if defined(USE_CGAL) && defined(USE_LIBFBI)
#error "Either USE_CGAL or USE_LIBFBI can be active, but not both!"
#endif

#ifdef USE_CGAL
#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2_algorithms.h>
#endif

namespace CoupledField{

class BoundingBoxAssociate{
public:
  static std::vector< std::pair<UInt, UInt > > AssociateEntities(const StdVector<Elem*>& entList1, const StdVector<Elem*>& entList2,
                                                                 Grid* g1, Grid* g2,Double tolerance=1e-6);

  static std::vector< std::pair<UInt, UInt > > AssociateEntities(const StdVector< Vector<Double> >& points, const StdVector<Elem*>& elemList,
                                                                 Grid* elGrid, Double tolerance=1e-6 );


#ifdef USE_CGAL

  protected:
    //! Define 3-dimensional bounding box
    typedef CGAL::Bbox_3 BBox3D;

    //! Define box handler, which additionally stores an index
    typedef CGAL::Box_intersection_d
        ::Box_with_handle_d<double,3,const UInt*> HandleBox;

    //! Define box handler just with an ID
    typedef CGAL::Box_intersection_d
        ::Box_d<double,3> IdBox;

    //! \param coords(in) The vector of global coordinates
    //! \param id(in) An identifier for this specific coordinate (e.g. index in a vector)
    //! \param tol(in) Tolerance in meters which determines the size of the bounding box
    static HandleBox CreateBoxFromCoord( const Vector<double>& coords,
                                         UInt *id, Double tol = 0.0 );

#elif USE_LIBFBI
    //ToBeCompleted
#endif

};

}



#endif /* SOURCE_DOMAIN_MESH_MESHUTILS_ENTITYASSOCIATION_HH_ */
