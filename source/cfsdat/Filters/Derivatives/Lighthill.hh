// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     Lighthill.hh
 *       \brief    <Description>
 *
 *       \date     Oct 4, 2016
 *       \author   kroppert
 */
//================================================================================================

#pragma once


#include "MeshBasedDerivative.hh"
#include "DataInOut/SimInput.hh"
#include <boost/tr1/type_traits.hpp>
#include <def_use_cgal.hh>
#include <def_use_flann.hh>
#include <cfsdat/Utils/Point.hh>

namespace CFSDat{



class Lighthill : public MeshBasedDerivative{

  struct DifferentiationStruct{
    CF::Vector<Double> localCoords;
    UInt tENum;
    UInt srcEqn;
    Double volume;

    DifferentiationStruct() : tENum(0),srcEqn(0),volume(.0){
      localCoords.Resize(3);
    }

    bool operator < (const DifferentiationStruct& str) const
    {
        return (srcEqn < str.srcEqn);
    }
  };

public:

  Lighthill(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~Lighthill();

  virtual bool Run();



protected:

  virtual void PrepareDifferentiation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();


  void CalcGradUScalarU(Vector<Double>& retVec,
                const Vector<Double> inVec);


  void CalcLocGradDerivativeCoefs(CF::Matrix<Double>& vec,
                                 CF::Vector<Double>& globPoint,
                                 CF::StdVector< Vector<Double> >& neighbors,
                                 CF::StdVector< Double >& l2Distances,
                                 CF::StdVector< Vector<Double> >& vectors,
                                 UInt numNN,
                                 Double alpha);

  void CalcRotU(Vector<Double>& retVec,
                const Vector<Double> inVec);


  void CalcLocRotDerivativeCoefs(CF::Matrix<Double>& vec,
                                 CF::Vector<Double>& globPoint,
                                 CF::StdVector< Vector<Double> >& neighbors,
                                 CF::StdVector< Double >& l2Distances,
                                 CF::StdVector< Vector<Double> >& vectors,
                                 UInt numNN,
                                 Double alpha);

  void InterpolNodeToCenter(Vector<Double>& retVec,
                            const Vector<Double> inVec);


  void OmegaCrossU(const Vector<Double> Omega,
                              const Vector<Double> U,
                              Vector<Double>& retVec);


private:

  std::vector<DifferentiationStruct> derivData_;

  // Coordinates of input data
  CF::StdVector< CF::Vector<double> > sourceCoords_;

  // Coordinates of target data
  CF::StdVector< CF::Vector<double> > targetCoords_;

  // Dimension of input values (0=scalar, 1=two-dim vector, 2=three-dim vector).
  UInt inDim_;

  // Density, if not specified in xml-scheme it is automatically set to one
  Double density_;
};

}
