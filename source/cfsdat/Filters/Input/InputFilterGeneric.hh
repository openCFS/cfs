// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     GenericInputFilter.hh
 *       \brief    Input filter for most CFS++ input readers
 *
 *       \date     Aug 13, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef GENERICINPUTFILTER_HH_
#define GENERICINPUTFILTER_HH_



#include "cfsdat/Filters/BaseFilter.hh"
#include "DataInOut/SimInput.hh"

namespace CFSDat{

class InputFilterGeneric : public BaseFilter {

public:
  InputFilterGeneric();

  virtual ~InputFilterGeneric();

  virtual void Init();

  virtual void Run();

protected:

  shared_ptr<CoupledField::SimInput> inFile_;

};

}


#endif /* GENERICINPUTFILTER_HH_ */
