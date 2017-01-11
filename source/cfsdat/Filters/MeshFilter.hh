// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     MeshFilter.hh
 *       \brief    <Description>
 *
 *       \date     Jan 11, 2017
 *       \author   kroppert
 */
//================================================================================================

#pragma once

#include <Filters/BaseMeshFilterType.hh>

namespace CFSDat{


//! Base class for Grid based interpolation- or differentiation-schemes

//! Some procedures require mesh results for
//! input and output. This class serves as a base for algorithms
//! requiring this kind of result representation
class MeshFilter : public BaseMeshFilterType{

public:

  //! Constructor, which reads the input- and output-result names from the xml-file and stores it
  //! in a string for further checking, if the "connectivity" of filters is correct
  //! Also the source- and target-regions are read from the xml-file
  //! Last step of this constructor is to read the target-mesh and create
  //! a new GridCFS-object
  MeshFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~MeshFilter(){
    delete trgGrid_;
  }

  virtual bool Run() = 0;

  virtual void FinishInit();

protected:

  virtual void PrepareCalculation() = 0;

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

