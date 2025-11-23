// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     KNNIndexSearch.cc
 *       \brief    Implementation für nearest neighbour based search with index association
 *
 *       \date     May 22, 2018
 *       \author   mtautz
 */
//================================================================================================

#include "KNNIndexSearch.hh"
#include "MatVec/Vector.hh"
#include <set>
// #include <algorithm>
#include <cmath>


namespace CFSDat{

KNNIndexSearch::KNNIndexSearch(){
  builtTree_ = false;
}

KNNIndexSearch::~KNNIndexSearch(){
  if (builtTree_) {
    delete tree_; 
  }
}

void KNNIndexSearch::AddPoint(Vector<Double>& coord, UInt idx) {
  AddPoint(coord[0],coord[1],coord[2],idx);
}

void KNNIndexSearch::AddPoint(Double x, Double y, Double z, UInt idx) {
#ifdef USE_CGAL 
  Point_3 point(x,y,z);
  points_.Push_back(point);
  indices_.Push_back(idx);
#else
  EXCEPTION("KNNIndexSearch only works with CGAL=on");
#endif  
}
  
void KNNIndexSearch::BuildTree() {
#ifdef USE_CGAL 
  tree_ = new Tree(boost::make_zip_iterator(boost::make_tuple( points_.Begin(),indices_.Begin() )),
                   boost::make_zip_iterator(boost::make_tuple( points_.End(),indices_.End() ) )
  );
  builtTree_ = true;
#else
  EXCEPTION("KNNIndexSearch only works with CGAL=on");
#endif  
}

void KNNIndexSearch::QueryPoint(Vector<Double>& coord, StdVector<UInt>& nearestIndices, StdVector<Double>& distances, UInt numNeighbours) {
  QueryPoint(coord[0],coord[1],coord[2],nearestIndices,distances,numNeighbours);
}

void KNNIndexSearch::QueryPoint(Double x, Double y, Double z, StdVector<UInt>& nearestIndices, StdVector<Double>& distances, UInt numNeighbours) {
  nearestIndices.Clear();
  distances.Clear();
#ifdef USE_CGAL 
  Point_3 query(x,y,z);
  Distance tr_dist;
  K_neighbor_search search(*tree_, query, numNeighbours);
  for(K_neighbor_search::iterator it = search.begin(); it != search.end(); it++) {
    nearestIndices.Push_back(boost::get<1>(it->first));
    distances.Push_back(tr_dist.inverse_of_transformed_distance(it->second));
  }
#else
  EXCEPTION("KNNIndexSearch only works with CGAL=on");
#endif  
}

UInt KNNIndexSearch::QueryPoint(Vector<Double>& coord) {
  return QueryPoint(coord[0],coord[1],coord[2]);
}
  
UInt KNNIndexSearch::QueryPoint(Double x, Double y, Double z) {
  StdVector<UInt> nearestIndices;
  StdVector<Double> distances;
  QueryPoint(x, y, z, nearestIndices, distances, 1);
  return nearestIndices[0];
}

  
}

