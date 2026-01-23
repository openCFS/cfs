// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     mathParserOMP.cc
 *       \brief    <Description>
 *
 *       \date     Mar 3, 2015
 *       \author   ahueppe
 */
//================================================================================================


#include "mathParserOMP.hh"

namespace CoupledField{

MathParserOMP::MathParserOMP(){

  if(omp_get_num_threads()>1){
    EXCEPTION("MathParserOMP created from parallel region. This may not happen. Use standard Mathparser!")
  }

  //assume, constructor is called from
  //serial part of code and set isShaerd to true
  isShared_=true;
  threadHandles_[GLOB_HANDLER].Mine(0) = GLOB_HANDLER;
  for(UInt aT=1;aT<threadHandles_[GLOB_HANDLER].GetNumSlots();aT++){
    threadHandles_[GLOB_HANDLER].Mine(aT) = MathParser::GetNewHandle(false);
  }
}

MathParserOMP::~MathParserOMP(){
}

unsigned int MathParserOMP::GetNewHandle( bool setDefaults ){
  if(omp_get_num_threads()>1){
    EXCEPTION("MathParserOMP::GetNewHandle called from parallel region. Not safe!");
  }
  unsigned int masterHandle = MathParser::GetNewHandle(setDefaults);
  //assign thread 0 the first handle
  threadHandles_[masterHandle].Mine(0) = masterHandle;

  //now we create the corresponding slave handles
  for(UInt aT=1;aT<threadHandles_[masterHandle].GetNumSlots();aT++){
    threadHandles_[masterHandle].Mine(aT) = MathParser::GetNewHandle(setDefaults);
  }

  return masterHandle;
}

void MathParserOMP::ReleaseHandle( unsigned int handle ){
  if(omp_get_num_threads()>1){
    EXCEPTION("MathParserOMP::GetNewHandle called from parallel region. Not safe!");
  }

  //here we need to clear each mathParser handle
  for(UInt aT=0;aT<threadHandles_[handle].GetNumSlots();aT++){
    MathParser::ReleaseHandle(threadHandles_[handle].Mine(aT));
  }
}

void MathParserOMP::SetExpr( unsigned int handle, const std::string &expr ){
  if(omp_get_num_threads()>1){
    EXCEPTION("MathParserOMP::GetNewHandle called from parallel region. Not safe!");
  }

  for(UInt aT=0;aT<threadHandles_[handle].GetNumSlots();aT++){
    MathParser::SetExpr(threadHandles_[handle].Mine(aT),expr);
  }
}

Double MathParserOMP::Eval( unsigned int handle ){
  return MathParser::Eval(threadHandles_[handle].Mine());
}

void MathParserOMP::EvalVector( unsigned int handle, Vector<Double>& vec ){
  MathParser::EvalVector(threadHandles_[handle].Mine(),vec);
}

void MathParserOMP::EvalDivVector( unsigned int handle, Double& divergence ){
  MathParser::EvalDivVector(threadHandles_[handle].Mine(),divergence);
}

void MathParserOMP::EvalMatrix( unsigned int handle, Matrix<Double>& matrix,
             UInt numRows, UInt numCols ){
  MathParser::EvalMatrix(threadHandles_[handle].Mine(),matrix,numRows,numCols);
}

void MathParserOMP::Dump( std::ostream& os ){
  os << "=================================================================" << std::endl;
  os << "  SMP version of mathParser expressions copied to each thread    " << std::endl;
  os << "=================================================================" << std::endl;
  MathParser::Dump(os);
}

void MathParserOMP::SetValue( unsigned int handle,
           const std::string &varName,
           Double val ){
  if(omp_get_num_threads() == 1){
    //assume call from serial region and we set all!
    for(UInt aT=0;aT<threadHandles_[handle].GetNumSlots();aT++){
      MathParser::SetValue(threadHandles_[handle].Mine(aT),varName,val);
    }
  }else{
    //assume from parallel region
    MathParser::SetValue(threadHandles_[handle].Mine(),varName,val);
  }
}

void MathParserOMP::RegisterExternalVar( unsigned int handle,
                      const std::string& varName,
                      Double * ptVar ){
  if(omp_get_num_threads() == 1){
    //assume call from serial region
    for(UInt aT=0;aT<threadHandles_[handle].GetNumSlots();aT++){
      MathParser::RegisterExternalVar(threadHandles_[handle].Mine(aT),varName,ptVar);
    }
  }else{
    //assume from parallel region
    MathParser::RegisterExternalVar(threadHandles_[handle].Mine(),varName,ptVar);
  }
}

void MathParserOMP::SetCoordinates( unsigned int handle,
                         const CoordSystem &coosy,
                         const Vector<Double> &globCoord ){
  // Performance fix: Only set coordinates for the current thread's slot.
  // Coordinates are per-element and set right before evaluation, so each thread
  // will set its own coordinates when needed. The old code did N× work when
  // CFS_NUM_THREADS=N, causing massive slowdown in line search evaluations.
  MathParser::SetCoordinates(threadHandles_[handle].Mine(),coosy,globCoord);
}

boost::signals2::connection
MathParserOMP::AddExpChangeCallBack( const MathParserSignal::slot_function_type
                          &subscriber,
                          unsigned int handle ){
  //only master
  if(omp_get_num_threads() != 1)
    EXCEPTION("Trying to add mathParser callback from parallel region. Not Supported!");

  return MathParser::AddExpChangeCallBack(subscriber,threadHandles_[handle].Mine());

}

}

