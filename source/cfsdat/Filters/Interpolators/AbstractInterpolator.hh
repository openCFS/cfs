// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     AbstractInterpolator.hh
 *       \brief    Abstract interpolation filter, the interpolation matrix has to be filledby implementing FillMatrix in a dpeneding class 
 *
 *       \date     Apr 4, 2018
 *       \author   matutz
 */
//================================================================================================

#pragma once //instead of the #ifndef #def ...


#include "DataInOut/SimInput.hh"
#include <Filters/MeshFilter.hh>
#include "DatUtils/InterpolationMatrix.hh"




namespace CFSDat{

class AbstractInterpolator : public MeshFilter {
  
public:

  AbstractInterpolator(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan);

  virtual ~AbstractInterpolator();

protected:

  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);
  
  virtual void PrepareCalculation();
  
  //! Fill the matrix using the given lists of global source and target entity numbers
  virtual void FillMatrix(StdVector<CF::UInt>& globSrcEntity, StdVector<CF::UInt>& globTrgEntity) = 0;

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();



protected:

  //! name of the interpolator
  std::string interpolatorName_;

  //! number of euqations per entity
  UInt numEquPerEnt_;
  
  //! Entity map used for source values
  shared_ptr<EqnMapSimple> scrMap_;

  //! if true, then element values are the interpolation target
  bool useElemAsSource_;

  //! Entity map used for target values
  shared_ptr<EqnMapSimple> trgMap_;
  
  //! if true, then element values are the interpolation target
  bool useElemAsTarget_;

  //! for the mesh-check this mesh also needs to be stored, trgGrid_ is
  //! stored in MeshFilter
  Grid* inGrid_;
  
  //! used interpolatioon matrix
  InterpolationMatrix matrix_;
  
  bool verboseSum_;
  
  virtual void PutIntoMatrix(UInt out, UInt in, Double factor);
  
  //! Global Factor for scaling the result
  Double globalFactor_;

};

}
