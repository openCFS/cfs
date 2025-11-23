// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     KNNSearch.cc
 *       \brief    Class for nearest-neighbor search based on knn-trees, using CGAL or FLANN
 *
 *       \date     Nov 23, 2016
 *       \author   kroppert
 */
//================================================================================================

#include "KNNSearch.hh"
#include <def_use_cgal.hh>
#include <cfsdat/DatUtils/Point.hh>
#include "MatVec/Vector.hh"

namespace CFSDat{
KNNSearch::KNNSearch(): trgGrid(){}
KNNSearch::~KNNSearch(){}

#ifdef USE_CGAL

void KNNSearch::KNN_CGAL_Interpolation(const Vector<Double> globPoint,
                                        StdVector< Vector<Double> >& neighbors,
                                        StdVector< Double >& l2Distances,
                                        StdVector< Vector<Double> >& vectors,
                                        UInt numNN){
  CGAL::Point query(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  if(globPoint.GetSize() == 2)
  {
    CGAL::Point query2(globPoint[0], globPoint[1], 0.0, 0.0, 0.0, 0.0);
    query = query2;
  }
  else
  {
    CGAL::Point query3(globPoint[0], globPoint[1], globPoint[2], 0.0, 0.0, 0.0);
    query = query3;
  }

  K_neighbor_search search(*searchTree_.get(), query, numNN);

  K_neighbor_search::iterator it = search.begin();

  if(it == search.end())
  {
    EXCEPTION("Could not find a nearest neighbor for " << globPoint << "!");
  }

  UInt nn = std::distance(search.begin(), search.end());
  neighbors.Resize(nn);
  l2Distances.Resize(nn);
  vectors.Resize(nn);

  for(UInt i=0 ; it != search.end(); ++it, i++) {
    l2Distances[i] = std::sqrt(it->second);
    neighbors[i].Resize(globPoint.GetSize());
    vectors[i].Resize(globPoint.GetSize());

    if(globPoint.GetSize() == 2)
    {
      it->first.vx(vectors[i][0]);
      it->first.vy(vectors[i][1]);
      neighbors[i][0] = it->first.x();
      neighbors[i][1] = it->first.y();
    }
    else
    {
      it->first.vx(vectors[i][0]);
      it->first.vy(vectors[i][1]);
      it->first.vz(vectors[i][2]);

      neighbors[i][0] = it->first.x();
      neighbors[i][1] = it->first.y();
      neighbors[i][2] = it->first.z();
    }

  }
}


void KNNSearch::KNN_CGAL_Differentiation(const Vector<Double> globPoint,
                                         StdVector< Vector<Double> >& neighbors,
                                         StdVector< Double >& l2Distances,
                                         StdVector< Vector<Double> >& vectors,
                                         UInt numNN){
  CGAL::Point query(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  if(globPoint.GetSize() == 2)
  {
    CGAL::Point query2(globPoint[0], globPoint[1], 0.0, 0.0, 0.0, 0.0);
    query = query2;
  }
  else
  {
    CGAL::Point query3(globPoint[0], globPoint[1], globPoint[2], 0.0, 0.0, 0.0);
    query = query3;
  }

  K_neighbor_search search(*searchTree_.get(), query, numNN);

  K_neighbor_search::iterator it = search.begin();

  if(it == search.end())
  {
    EXCEPTION("Could not find a nearest neighbor for " << globPoint << "!");
  }

  UInt nn = std::distance(search.begin(), search.end());
  neighbors.Resize(nn);
  l2Distances.Resize(nn);
  vectors.Resize(nn);

  for(UInt i=0 ; it != search.end(); ++it, i++) {
    l2Distances[i] = std::sqrt(it->second);
    neighbors[i].Resize(globPoint.GetSize());
    vectors[i].Resize(globPoint.GetSize());

    if(globPoint.GetSize() == 2)
    {
      it->first.vx(vectors[i][0]);
      it->first.vy(vectors[i][1]);
      neighbors[i][0] = it->first.x();
      neighbors[i][1] = it->first.y();
    }
    else
    {
      it->first.vx(vectors[i][0]);
      it->first.vy(vectors[i][1]);
      it->first.vz(vectors[i][2]);

      neighbors[i][0] = it->first.x();
      neighbors[i][1] = it->first.y();
      neighbors[i][2] = it->first.z();
    }

  }
}



void KNNSearch::ReadScatteredData_Interpolation(const StdVector< Vector<double> > sourceCoords,
                                                const UInt& inDim,
                                                const Grid* trgGrid,
                                                StdVector< Vector<Double> > scatteredData){

  UInt n = sourceCoords.GetSize();


  UInt i=0;

  std::vector<CGAL::Point> points;
  points.resize(n);

  for( i=0; i<n; i++)
  {

    if(trgGrid->GetDim() == 2){
      if(inDim == 0){
        points[i] = (CGAL::Point(sourceCoords[i][0],
                                  sourceCoords[i][1],
                                  0.0,
                                  scatteredData[i][0],
                                  0.0,
                                  0.0));
      }else{
        points[i] = (CGAL::Point(sourceCoords[i][0],
                                  sourceCoords[i][1],
                                  0.0,
                                  scatteredData[i][0],
                                  scatteredData[i][1],
                                  0.0));
      }
    }
    if(trgGrid->GetDim() == 3){
      if(inDim == 0){
        points[i] = (CGAL::Point(sourceCoords[i][0],
                                  sourceCoords[i][1],
                                  sourceCoords[i][2],
                                  scatteredData[i][0],
                                  0.0,
                                  0.0));
      }else{
        points[i] = (CGAL::Point(sourceCoords[i][0],
                                sourceCoords[i][1],
                                sourceCoords[i][2],
                                scatteredData[i][0],
                                scatteredData[i][1],
                                scatteredData[i][2]));
      }
    }
  }
  searchTree_.reset(new Tree(points.begin(), points.end()));

}


void KNNSearch::ReadScatteredData_Grad(const StdVector< Vector<double> > sourceCoords,
                                        const UInt& inDim,
                                        const Grid* trgGrid,
                                        StdVector< Vector<Double> > scatteredData)
{
  UInt n = sourceCoords.GetSize();

  UInt i=0;

  std::vector<CGAL::Point> points;
  points.resize(n);

  for( i=0; i<n; i++)
  {

    if(trgGrid->GetDim() == 2){
      if(inDim == 0){
        points[i] = (CGAL::Point(sourceCoords[i][0],
            sourceCoords[i][1],
            0.0,
            scatteredData[i][0],
            0.0,
            0.0));
      }else{
        points[i] = (CGAL::Point(sourceCoords[i][0],
            sourceCoords[i][1],
            0.0,
            scatteredData[i][0],
            0.0,
            0.0));
      }
    }
    if(trgGrid->GetDim() == 3){
      if(inDim == 0){
        points[i] = (CGAL::Point(sourceCoords[i][0],
            sourceCoords[i][1],
            sourceCoords[i][2],
            scatteredData[i][0],
            0.0,
            0.0));
      }else{
        points[i] = (CGAL::Point(sourceCoords[i][0],
            sourceCoords[i][1],
            sourceCoords[i][2],
            scatteredData[i][0],
            0.0,
            0.0));
      }
    }

  }
  searchTree_.reset(new Tree(points.begin(), points.end()));
}



void KNNSearch::ReadScatteredData_Div(const StdVector< Vector<double> > sourceCoords,
                                        const UInt& inDim,
                                        const Grid* trgGrid,
                                        StdVector< Vector<Double> > scatteredData)
{

  UInt n = sourceCoords.GetSize();

  UInt i=0;

  std::vector<CGAL::Point> points;
  points.resize(n);

  for( i=0; i<n; i++)
  {

    if(trgGrid->GetDim() == 2){
      if(inDim == 0){
        points[i] = (CGAL::Point(sourceCoords[i][0],
            sourceCoords[i][1],
            0.0,
            scatteredData[i][0],
            0.0,
            0.0));
      }else{
        points[i] = (CGAL::Point(sourceCoords[i][0],
            sourceCoords[i][1],
            0.0,
            scatteredData[i][0],
            scatteredData[i][1],
            0.0));
      }
    }
    if(trgGrid->GetDim() == 3){
      if(inDim == 0){
        points[i] = (CGAL::Point(sourceCoords[i][0],
            sourceCoords[i][1],
            sourceCoords[i][2],
            scatteredData[i][0],
            0.0,
            0.0));
      }else{
        points[i] = (CGAL::Point(sourceCoords[i][0],
            sourceCoords[i][1],
            sourceCoords[i][2],
            scatteredData[i][0],
            scatteredData[i][1],
            scatteredData[i][2]));
      }
    }

  }
    searchTree_.reset(new Tree(points.begin(), points.end()));
  }



void KNNSearch::ReadScatteredData_Curl(const StdVector< Vector<double> > sourceCoords,
                                        const UInt& inDim,
                                        const Grid* trgGrid,
                                        StdVector< Vector<Double> > scatteredData)
{

  UInt n = sourceCoords.GetSize();

  UInt i=0;

  std::vector<CGAL::Point> points;
  points.resize(n);

  for( i=0; i<n; i++)
  {

    if(trgGrid->GetDim() == 2){
      if(inDim == 0){
        points[i] = (CGAL::Point(sourceCoords[i][0],
            sourceCoords[i][1],
            0.0,
            scatteredData[i][0],
            0.0,
            0.0));
      }else{
        points[i] = (CGAL::Point(sourceCoords[i][0],
            sourceCoords[i][1],
            0.0,
            scatteredData[i][0],
            scatteredData[i][1],
            0.0));
      }
    }
    if(trgGrid->GetDim() == 3){
      if(inDim == 0){
        points[i] = (CGAL::Point(sourceCoords[i][0],
            sourceCoords[i][1],
            sourceCoords[i][2],
            scatteredData[i][0],
            0.0,
            0.0));
      }else{
        points[i] = (CGAL::Point(sourceCoords[i][0],
            sourceCoords[i][1],
            sourceCoords[i][2],
            scatteredData[i][0],
            scatteredData[i][1],
            scatteredData[i][2]));
      }
    }

  }
  searchTree_.reset(new Tree(points.begin(), points.end()));

  }




void KNNSearch::ReadScatteredData_Lighthill(const StdVector< Vector<double> > sourceCoords,
                                            const UInt& inDim,
                                            const Grid* trgGrid,
                                            StdVector< Vector<Double> > scatteredData){
  UInt n = sourceCoords.GetSize();
  UInt i=0;

    std::vector<CGAL::Point> points;
    points.resize(n);

    for( i=0; i<n; i++)
    {

      if(trgGrid->GetDim() == 2){
        if(inDim == 0){
          points[i] = (CGAL::Point(sourceCoords[i][0],
                                   sourceCoords[i][1],
                                   0.0,
                                   scatteredData[i][0],
                                   0.0,
                                   0.0));
          }else{
            points[i] = (CGAL::Point(sourceCoords[i][0],
                                     sourceCoords[i][1],
                                     0.0,
                                     scatteredData[i][0],
                                     scatteredData[i][1],
                                     0.0));
            }
      }
      if(trgGrid->GetDim() == 3){
        if(inDim == 0){
          points[i] = (CGAL::Point(sourceCoords[i][0],
                                   sourceCoords[i][1],
                                   sourceCoords[i][2],
                                   scatteredData[i][0],
                                   0.0,
                                   0.0));
        }else{
          points[i] = (CGAL::Point(sourceCoords[i][0],
                                   sourceCoords[i][1],
                                   sourceCoords[i][2],
                                   scatteredData[i][0],
                                   scatteredData[i][1],
                                   scatteredData[i][2]));
      }
    }

    }
    searchTree_.reset(new Tree(points.begin(), points.end()));
  }






#else
EXCEPTION("CoefFunctionScatteredData needs to be compiled with USE_CGAL=ON!");
#endif


}
