// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     VolumeMultiplication.hh
 *       \brief    <Description>
 *
 *       \date     Apr 23, 2018
 *       \author   mtautz
 */
//================================================================================================

#ifndef VOLUMEMULTIPLICATION_HH_
#define VOLUMEMULTIPLICATION_HH_

#include "cfsdat/Filters/BaseFilter.hh"

#include <boost/bimap.hpp>
#include <set>

namespace CFSDat{
class VolumeMultiplication : public BaseFilter{
public:


  //!  Constructor.
  VolumeMultiplication(UInt numWorkers, PtrParamNode config, PtrResultManager resMan);

  virtual ~VolumeMultiplication();


protected:

  virtual void PrepareCalculation();
  
  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);

  //! Setup of Results of the filter
  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

  uuids::uuid resId;

  std::string resName;
  std::string outName;
  uuids::uuid outId;

private:

  UInt usedElems_;
  UInt numDofs_;
  StdVector<Double> volume_;
  
};

}
#endif /* VOLUMEMULTIPLICATION_HH_ */
