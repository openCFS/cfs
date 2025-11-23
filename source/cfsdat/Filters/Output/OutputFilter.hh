// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     OutputFilterGeneric.hh
 *       \brief    <Description>
 *
 *       \date     Oct 26, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef OUTPUTFILTER_HH_
#define OUTPUTFILTER_HH_

#include "cfsdat/Filters/BaseFilter.hh"
#include "DataInOut/SimOutput.hh"
#include "cfsdat/DatUtils/ResultManager.hh"

namespace CFSDat{

class OutputFilter : public BaseFilter {

public:
  OutputFilter(UInt numWorkers, CoupledField::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~OutputFilter();

  virtual bool Run();

  virtual void PrepareCalculation();

protected:

  virtual void AddOutput(str1::shared_ptr<BaseFilter> filt){
    EXCEPTION("Filter Error : " << this->filterId_ << "\n An output filter may not have another filter output. PassThrough is not supported, use branching.")
  }

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();


  str1::shared_ptr<CoupledField::SimOutput> outFile_;

  std::map<UInt,Double>::iterator aStepIter_;
private:
  //! dummy info node untill we get a clean log
  CoupledField::PtrParamNode infoNode_;
};

}


#endif /* OUTPUTFILTERGENERIC_HH_ */
