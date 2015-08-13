// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     BaseFilter.hh
 *       \brief    Basic interface for CFSDat filters
 *
 *       \date     Aug 13, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef BASEFILTER_HH_
#define BASEFILTER_HH_

namespace CFSDat{

class BaseFilter{

public:
  BaseFilter(){

  }

  virtual ~BaseFilter(){

  }

  virtual void Init() = 0;

  virtual void Run() = 0;

protected:
private:

};

}



#endif /* BASEFILTER_HH_ */
