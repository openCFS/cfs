// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     ResultManager.hh
 *       \brief    Nearest Neighbour Search with CGAL. Usage:
 * 1. Add some points with an unique index to be found by using the AddPoint() functions
 * 2. BuiltTree(); build the search tree
 * 3. Use QueryPoint() functions to get the indices of the nearest points and their distances to the queried point
 *
 *       \date     May 22, 2018
 *       \author   mtautz
 */
//================================================================================================

#ifndef KNNINDEXSEARCH_HH_
#define KNNINDEXSEARCH_HH_



#include <utility>
#include <vector>

//CFS
#include "MatVec/Vector.hh"
#include "Utils/StdVector.hh"


#include <def_use_cgal.hh>
#ifdef USE_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/basic.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <boost/iterator/zip_iterator.hpp>
#endif


#define DEBUG



namespace CFSDat{

#ifdef USE_CGAL  
// cgal copy and paste stuff for association of points to CF::UInt from
// http://doc.cgal.org/Manual/4.2/doc_html/cgal_manual/Spatial_searching/Chapter_main.html#Subsection_61.3.1
//
// later on this should be refactored into KNNSearch.hh and KNNSearch.cc
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef boost::tuple<Point_3,UInt> Point_and_int;

//definition of the property map
struct My_point_property_map{
  typedef Point_3 value_type;
  typedef const value_type& reference;
  typedef const Point_and_int& key_type;
  typedef boost::lvalue_property_map_tag category;
};

//get function for the property map
inline My_point_property_map::reference
get(My_point_property_map,My_point_property_map::key_type p)
{return boost::get<0>(p);}

typedef CGAL::Random_points_in_cube_3<Point_3>                                          Random_points_iterator;
typedef CGAL::Search_traits_3<Kernel>                                                   Traits_base;
typedef CGAL::Search_traits_adapter<Point_and_int,My_point_property_map,Traits_base>    Traits;


typedef CGAL::Orthogonal_k_neighbor_search<Traits>                      K_neighbor_search;
typedef K_neighbor_search::Tree                                         Tree;
typedef K_neighbor_search::Distance                                     Distance;
#else
   typedef UInt Point_3;
   typedef UInt Tree;
#endif 
  
class KNNIndexSearch {

public:

  KNNIndexSearch();

  ~KNNIndexSearch();
  
  void AddPoint(Vector<Double>& coord, UInt idx);

  void AddPoint(Double x, Double y, Double z, UInt idx);
  
  void BuildTree();
  
  void QueryPoint(Vector<Double>& coord, StdVector<UInt>& nearestIndices, StdVector<Double>& distances, UInt numNeighbours);
  
  void QueryPoint(Double x, Double y, Double z, StdVector<UInt>& nearestIndices, StdVector<Double>& distances, UInt numNeighbours);
  
  UInt QueryPoint(Vector<Double>& coord);
  
  UInt QueryPoint(Double x, Double y, Double z);
  
private:


  StdVector<Point_3> points_;
  StdVector<UInt> indices_;

  bool builtTree_;
  Tree* tree_;  

};

}



#endif /* KNNINDEXSEARCH_HH_ */
