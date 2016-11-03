// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     MeshBasedDerivative.hh
 *       \brief    <Description>
 *
 *       \date     Oct 4, 2016
 *       \author   kroppert
 */
//================================================================================================

#pragma once

#include "BaseDerivativeFilter.hh"

namespace CFSDat{


//! Base class for Grid based differentiation schemes

//! Many derivative procedures require mesh results for
//! input and output. This class serves as a base for algorithms
//! requiring this kind of result representation
class MeshBasedDerivative : public BaseDerivativeFilter{

public:
  MeshBasedDerivative(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~MeshBasedDerivative(){
    delete trgGrid_;
  }

  virtual bool Run() = 0;

  virtual void FinishInit();

protected:

  virtual void PrepareDifferentiation() = 0;

  virtual ResultIdList SetUpstreamResults()=0;

  virtual void AdaptFilterResults()=0;

  //TODO for now copy paste code from inputfilter... not very nice
  void CreateDummyCfsParamNode();

  /// XML node for creating a CFS input object
  PtrParamNode dummyXMLNode_;

  ///Set of input/src regions
  std::set<std::string> srcRegions_;

  ///Set of target regions
  std::set<std::string> trgRegions_;

  //Siminput pointer to target mesh
  str1::shared_ptr<CoupledField::SimInput> trgMeshInp_;

  ///Pointer to new grid object
  Grid* trgGrid_;

};

}
