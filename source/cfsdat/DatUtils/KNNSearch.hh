// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     KNNSearch.hh
 *       \brief    Class for nearest-neighbor search based on knn-trees, using CGAL
 *
 *       \date     Nov 23, 2016
 *       \author   kroppert
 */
//================================================================================================

#pragma once

#include <def_use_cgal.hh>
#include <def_use_flann.hh>
#include <cfsdat/DatUtils/Point.hh>
#include <Filters/BaseMeshFilterType.hh>
#include "MatVec/Vector.hh"

namespace CFSDat{

class KNNSearch{
public:
  KNNSearch();
  ~KNNSearch();


#ifdef USE_CGAL
  boost::shared_ptr<Tree> searchTree_;

  // This knn-search method is used by RBF, NearestNeighborInterpolator
  void KNN_CGAL_Interpolation(const Vector<Double> globPoint,
                                StdVector< Vector<Double> >& neighbors,
                                StdVector< Double >& l2Distances,
                                StdVector< Vector<Double> >& vectors,
                                UInt numNeighbors);

  // This knn-search method is used by Div-, Curl- and Grad-Differentiator
  void KNN_CGAL_Differentiation(const Vector<Double> globPoint,
                                StdVector< Vector<Double> >& neighbors,
                                StdVector< Double >& l2Distances,
                                StdVector< Vector<Double> >& vectors,
                                UInt numNeighbors);



  // Read scattered data and prepare for KNN search
  // Different methods for Interpolators and Differentiator, due to the
  // different handling of scalar and vector values
    void ReadScatteredData_Interpolation(const StdVector< Vector<double> > sourceCoords,
                                         const UInt& inDim,
                                         const Grid* trgGrid,
                                         StdVector< Vector<Double> > scatteredData);

    void ReadScatteredData_Grad(const StdVector< Vector<double> > sourceCoords,
                                const UInt& inDim,
                                const Grid* trgGrid,
                                StdVector< Vector<Double> > scatteredData);

    void ReadScatteredData_Div(const StdVector< Vector<double> > sourceCoords,
                                const UInt& inDim,
                                const Grid* trgGrid,
                                StdVector< Vector<Double> > scatteredData);


    void ReadScatteredData_Curl(const StdVector< Vector<double> > sourceCoords,
                                const UInt& inDim,
                                const Grid* trgGrid,
                                StdVector< Vector<Double> > scatteredData);


    void ReadScatteredData_Lighthill(const StdVector< Vector<double> > sourceCoords,
                                     const UInt& inDim,
                                     const Grid* trgGrid,
                                     StdVector< Vector<Double> > scatteredData);


#else
    EXCEPTION("CoefFunctionScatteredData needs to be compiled with USE_CGAL=ON!");
#endif

private:
    Grid* trgGrid;
};






}
