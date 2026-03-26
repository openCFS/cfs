#ifndef GRID_KD_TREE_HH
#define GRID_KD_TREE_HH

#include <array>
#include <nanoflann.hpp>

#include "MatVec/Vector.hh"
#include "Utils/StdVector.hh"
#include "Domain/ElemMapping/Elem.hh"

// we are already in the namespace CoupledField (and need to be)

/** When we have nanoflann (USE_NANOFLANN is ON by default).
 * We use it for nearest neighbor searches for nodes coordinates and Elem*.
 * nanoflann is header only and by the adaptors the data does not need to be copied,
 * we just build the index. */

/** Adaptor for nodes (GridCFS::coords_) */
struct NodeCoordsAdaptor 
{
  const StdVector<Vector<Double>>& pts;
  const int dim;

  inline size_t kdtree_get_point_count() const { return pts.GetSize(); }
  inline double kdtree_get_pt(size_t idx, int d) const { 
    return pts[idx][(unsigned int) d]; 
  }

  template<class BBOX>
  bool kdtree_get_bbox(BBOX&) const { return false; }
};

/** Adaptor for elements. Takes GridCFS::orderedElems_. 
 * Uses extended->barycenter, initialized by SetElementBarycenters() */
struct ElemCoordsAdaptor {
  const StdVector<Elem*>& elems;
  const int dim;

  inline size_t kdtree_get_point_count() const { return elems.GetSize(); }
  inline double kdtree_get_pt(size_t idx, int d) const {
    if(Elem::shapes[elems[idx]->type].dim == (unsigned int) dim)
      return elems[idx]->extended->barycenter.data[(unsigned int) d];
    else
      return 1e13; // far away point for surface elements 
  }

  template<class BBOX>
  bool kdtree_get_bbox(BBOX&) const { return false; }
};


template<typename Adaptor>
struct GridKDTreeT {
  using Index = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, Adaptor>, Adaptor>;

  GridKDTreeT(Adaptor adaptor, unsigned int dim)
    : adaptor_(adaptor), // 50 as max leaf size searches still quite fast and saves indexing time 
      index_((int) dim, adaptor_, nanoflann::KDTreeSingleIndexAdaptorParams(50))
  {
    index_.buildIndex();
  }

  size_t FindNearest(const Vector<double>& c, unsigned int dim) const
  {
    return FindNearest(c.GetPointer());
  }

private:

  /** boilerplate code */ 
  size_t FindNearest(const double* q) const
  {
    size_t idx;
    double dist2;
    nanoflann::KNNResultSet<double> result(1);
    result.init(&idx, &dist2);
    index_.findNeighbors(result, q, nanoflann::SearchParameters());
    return idx;
  }

  Adaptor adaptor_;
  Index   index_;
};

/** there are dummy implementations of Node/ElemGridKDTree when compiled without nanoflann. */
struct NodeGridKDTree : GridKDTreeT<NodeCoordsAdaptor> {
  NodeGridKDTree(const StdVector<Vector<double>>& coords, unsigned int dim)
    : GridKDTreeT<NodeCoordsAdaptor>(NodeCoordsAdaptor{coords, (int) dim}, dim) {}
};

struct ElemGridKDTree : GridKDTreeT<ElemCoordsAdaptor> {
  ElemGridKDTree(const StdVector<Elem*>& elems, unsigned int dim)
    : GridKDTreeT<ElemCoordsAdaptor>(ElemCoordsAdaptor{elems, (int) dim}, dim) {}
};

#endif // GRID_KD_TREE_HH
