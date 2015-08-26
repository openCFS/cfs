// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     DataStructs.hh
 *       \brief    <Description>
 *
 *       \date     Aug 24, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef DATASTRUCTS_HH_
#define DATASTRUCTS_HH_

namespace CFSDat{

template<typename DType>
struct DataStruct{

  boost::shared_ptr< CF::Result<DType> > currentResult_;

  CF::Double timeFreqStep_;


};



}




#endif /* DATASTRUCTS_HH_ */
