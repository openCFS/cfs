// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     GaussDerivative.hh
 *       \brief    Filter for the Gauss Derivative Computation using a finite volume representation
 *
 *       \date     Apr 23, 2018
 *       \author   mtautz
 */
//================================================================================================

#ifndef GAUSSDERIVATIVE_HH_
#define GAUSSDERIVATIVE_HH_

#include "cfsdat/Filters/BaseFilter.hh"
#include "cfsdat/DatUtils/DerivativeMatrix.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/ElemMapping/Elem.hh"

#include <boost/bimap.hpp>
#include <set>

namespace CFSDat{
class GaussDerivative : public BaseFilter{
public:


  //!  Constructor.
  GaussDerivative(UInt numWorkers, PtrParamNode config, PtrResultManager resMan);

  virtual ~GaussDerivative();


protected:

  virtual void PrepareCalculation();
  
  //! Calculate the finite element volume using the elemet replacements of polyhedrons
  void CalcFeVolume(UInt usedElems, StdVector<UInt>& masterElemIdx);

  // compute the derivative matrix and calculates the finite volumes on demand
  void CreateGaussComputation(Grid::FiniteVolumeRepresentation& fvr, 
      DerivativeMatrix&  derivMatrix, UInt usedElems, 
      std::vector<bool>& isFaceInterior, std::vector<bool>& isFaceBoundary, std::vector<bool>& isFaceSlaveBoundary,
      StdVector<UInt>& masterElemIdx, StdVector<UInt>& faceMasterIdx, StdVector<UInt>& faceSlaveIdx,
      bool calcFvVolume, std::vector<bool>& hasFvVolume, StdVector<Double>& fvVolume);
  
  // divide the results by the cell volume
  template<UInt dim>
  void VolumeDivide(Double* out);
  
  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);

  //! Setup of Results of the filter
  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

  uuids::uuid resId;

  std::string resName;
  std::string outName;
  uuids::uuid outId;

private:
  
  typedef enum { DIVERGENCE, GRADIENT, CURL } DerivType;
  
  Grid* grid_;
  boost::shared_ptr<std::set<std::string> > regNames_;
  
  ResultInfo::EntryType inputType_;
  DerivType derivType_;
  bool zeroBC_;
  
  bool useFvVolume_;
  std::vector<bool> hasFvVolume_;
  StdVector<Double> fvVolume_;
  std::vector<bool> hasFeVolume_;
  StdVector<Double> feVolume_;
  
  //! matrix for the interior derivative
  DerivativeMatrix interiorMatrix_;
  
  //! matrix for the boundary condition
  DerivativeMatrix boundaryMatrix_;
  
  // derivtaive containing the matrices to be used
  GaussDerivative* usedDerivative_;
  
  // all other instances of this class
  static StdVector<GaussDerivative*> allDerivatives_;

};

}
#endif /* GAUSSDERIVATIVE_HH_ */
