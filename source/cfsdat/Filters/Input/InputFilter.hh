// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     GenericInputFilter.hh
 *       \brief    Input filter for most openCFS input readers
 *
 *       \date     Aug 13, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef INPUTFILTER_HH_
#define INPUTFILTER_HH_



#include "cfsdat/Filters/BaseFilter.hh"
#include "DataInOut/SimInput.hh"
#include "cfsdat/DatUtils/ResultManager.hh"
#include "cfsdat/DatUtils/DataStructs.hh"

namespace CFSDat{

class InputFilter : public BaseFilter {

public:
  InputFilter(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan);

  virtual ~InputFilter();

  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);


protected:
  virtual void AddInput(shared_ptr<BaseFilter> filt){
    EXCEPTION("An input filter may not have a source! Filter id: " << this->filterId_);
  }

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

  void CreateDummyCfsParamNode();

  void CreateAvailableResultInfos();

  PtrParamNode dummyXMLNode;

  shared_ptr<CoupledField::SimInput> inFile_;

  CF::Grid* ptGrid;

  StdVector<ExtendedResultInfo> availInputResults_;

  std::map<std::string, std::map<Double,UInt> > stepNumbers_;

  std::set<std::string> volRegions;

  std::set<std::string> surfRegions;
private:
  bool ranAlready_;
  
  bool staticTimeType_;
};

}


#endif /* GENERICINPUTFILTER_HH_ */
