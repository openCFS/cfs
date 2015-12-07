// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CentroidInterpolator.hh
 *       \brief    <Description>
 *
 *       \date     Nov 14, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SOURCE_CFSDAT_FILTERS_INTERPOLATORS_CENTROIDINTERPOLATOR_HH_
#define SOURCE_CFSDAT_FILTERS_INTERPOLATORS_CENTROIDINTERPOLATOR_HH_


#include "BaseInterpolationFilter.hh"
#include "DataInOut/SimInput.hh"

namespace CFSDat{

//! Class for choosing a very simple but fast interpolation approach
//! Upon initialization we determine cell centroids and volume from input
//! additionally we set the local coordinates accoring to source
//! during traversal, we just apply those loads
class CentroidInterpolator : public BaseInterpolationFilter{

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

  CentroidInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~CentroidInterpolator();

  virtual bool Run();

  virtual void FinishInit();


protected:

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

private:

  void CreateDummyCfsParamNode();

  PtrParamNode dummyXMLNode_;

  //Siminput pointer to target mesh
  str1::shared_ptr<CoupledField::SimInput> trgMeshInp_;

  ///Pointer to new grid object
  Grid* trgGrid_;

  ///Set of input/src regions
  std::set<std::string> srcRegions_;

  ///Set of target regions
  std::set<std::string> trgRegions_;

  std::vector<InpolationStruct> interpolData_;

};

}
#endif /* SOURCE_CFSDAT_FILTERS_INTERPOLATORS_CENTROIDINTERPOLATOR_HH_ */
