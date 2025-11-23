// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     Defines.hh
 *       \brief    <Description>
 *
 *       \date     Nov 19, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef CFSDAT_DEFINES_HH_
#define CFSDAT_DEFINES_HH_

#include "General/Exception.hh"
#include "Utils/StdVector.hh"
#include "General/Environment.hh"
#define BOOST_UUID_RANDOM_PROVIDER_FORCE_POSIX
#include <boost/uuid/uuid.hpp>
#include <boost/shared_ptr.hpp>
#include <sstream>
#include <string>
#if !defined(WIN32)
  #include <unistd.h>
#endif
#include <ios>
#include <cmath>
#include <fstream>

namespace CFSDat{

class BaseFilter;
class ResultManager;

//we'll have to see if this is ok...
namespace CF = CoupledField;
namespace str1 = boost;
namespace uuids = boost::uuids;
using namespace str1;

//make some CoupledField variables available in this namespace
typedef CoupledField::UInt    UInt;
typedef CoupledField::Double  Double;
typedef CoupledField::Integer  Integer;
typedef str1::shared_ptr<BaseFilter> FilterPtr;
typedef str1::shared_ptr<ResultManager> PtrResultManager;

}



#endif /* DEFINES_HH_ */
