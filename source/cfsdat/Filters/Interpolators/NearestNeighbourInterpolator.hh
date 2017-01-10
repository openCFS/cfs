// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     NearesNeighbourInterpolator.hh
 *       \brief    <Description>
 *
 *       \date     Aug 8, 2016
 *       \author   kroppert
 */
//================================================================================================

#pragma once //instead of the #ifndef #def ...


#include "MeshBasedInterpolator.hh"
#include "DataInOut/SimInput.hh"
#include <cfsdat/Utils/Point.hh>



namespace CFSDat{

//! Class for calculating a nearest neighbour interpolation using CGAL
//! for neighbour search


class NearestNeighbourInterpolator : public MeshBasedInterpolator{

  struct InpolationStruct{
    CF::Vector<Double> localCoords;
    UInt tENum;
    UInt srcEqn;
    Double volume;

    InpolationStruct() : tENum(0),srcEqn(0),volume(.0){
      localCoords.Resize(3);
    }

    bool operator < (const InpolationStruct& str) const
    {
      return (srcEqn < str.srcEqn);
    }
  };

public:

  NearestNeighbourInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~NearestNeighbourInterpolator();

  virtual bool Run();

protected:

  virtual void PrepareInterpolation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

  // Read scattered data
  void ReadScatteredData(CF::StdVector< CF::Vector<Double> > elemCentroids, CF::StdVector< CF::Vector<Double> > scatteredData);


private:
  void Interpolation(Vector<Double>& returnVec, CF::StdVector< CF::Vector<Double> >  scatteredData,
		  CF::Vector<Double> vec, str1::shared_ptr<EqnMapSimple> downMap);

  std::string CheckFilterResults();

  // Coordinates of input data
  CF::StdVector< CF::Vector<double> > sourceCoords_;

  // Coordinates of target data
  CF::StdVector< CF::Vector<double> > targetCoords_;

  //! Dimension of input values (0=scalar, 1=two-dim vector, 2=three-dim vector).
  UInt inDim_;

  //! Number of neighbor points to include in interpolation.
  UInt numNeighbors_;

  //! Exponent for calculation of interpolation weight function.
  Double p_;

  std::vector<InpolationStruct> interpolData_;
};

}
