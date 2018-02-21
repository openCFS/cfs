// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     TemporalBlendingFilter.hh
 *       \brief    <Description>
 *
 *       \date     Apr 18, 2017
 *       \author   kroppert
 */
//================================================================================================

#ifndef TEMPBLENDFILTER_HH_
#define TEMPBLENDFILTER_HH_

#include "cfsdat/Filters/BaseFilter.hh"
#include <boost/bimap.hpp>

namespace CFSDat{
class TemporalBlendFilter : public BaseFilter{
public:


  //!  Constructor.
  TemporalBlendFilter(UInt numWorkers, PtrParamNode config, PtrResultManager resMan);

  virtual ~TemporalBlendFilter();


protected:

  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);

  //! Setup of Results of the filter
  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

  uuids::uuid resId;

  std::string resName;
  std::string outName;
  uuids::uuid outId;

private:

 //stores for each filter result the dependend upstream ids
 std::map<std::string, uuids::uuid> filtResToUpResIds_;

 std::string expr_;

 MathParser* mp_;

};

}
#endif /* TEMPBLENDFILTER_HH_ */
