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

  //check if constructor is called from serial code
  if(omp_get_num_threads() == 1){
    //assume, constructor is called from
    //serial part of code and set isShaerd to true
    isShared_=true;
    //obtain maximum number of threads
    //(maybe useful to define cfs variables for that)
    numThreads_ = NUM_CFS_THREADS;

    threadHandles_[GLOB_HANDLER].Resize(numThreads_);
    threadHandles_[GLOB_HANDLER][0] = GLOB_HANDLER;
    for(UInt aT=1;aT<numThreads_;aT++){
       threadHandles_[GLOB_HANDLER][aT] = MathParser::GetNewHandle(false);
    }

  }else if(omp_get_num_threads()>1){
    EXCEPTION("MathParserOMP created from parallel region. This may not happen.")
  }else{
    isShared_ = false;
    numThreads_ = 1;
    threadHandles_[GLOB_HANDLER].Resize(numThreads_);
    threadHandles_[GLOB_HANDLER][0] = GLOB_HANDLER;
  }
}

MathParserOMP::~MathParserOMP(){

}

MathParser::HandleType MathParserOMP::GetNewHandle( bool setDefaults ){
  MathParser::HandleType masterHandle = MathParser::GetNewHandle(setDefaults);

  if(omp_get_num_threads() != 1){
    EXCEPTION("MathParserOMP::GetNewHandle called from parallel region. Not safe!");
  }

  UInt aThread = omp_get_thread_num();

  //now we create the corresponding slave handles
  threadHandles_[masterHandle].Resize(numThreads_);
  threadHandles_[masterHandle][aThread] = masterHandle;
  for(UInt aT=0;aT<numThreads_;aT++){
    if(aT!=aThread)
      threadHandles_[masterHandle][aT] = MathParser::GetNewHandle(setDefaults);
  }

  return masterHandle;
}

void MathParserOMP::ReleaseHandle( HandleType handle ){
  for(UInt aT=0;aT<numThreads_;aT++){
    MathParser::ReleaseHandle(threadHandles_[handle][aT]);
  }
  threadHandles_[handle].Clear(true);
}

void MathParserOMP::SetExpr( HandleType handle, const std::string &expr ){
  for(UInt aT=0;aT<numThreads_;aT++){
    MathParser::SetExpr(threadHandles_[handle][aT],expr);
  }
}

Double MathParserOMP::Eval( HandleType handle ){
  UInt aThread = omp_get_thread_num();
  return MathParser::Eval(threadHandles_[handle][aThread]);
}

void MathParserOMP::EvalVector( HandleType handle, Vector<Double>& vec ){
  UInt aThread = omp_get_thread_num();
  MathParser::EvalVector(threadHandles_[handle][aThread],vec);
}

void MathParserOMP::EvalDivVector( HandleType handle, Double& divergence ){
  UInt aThread = omp_get_thread_num();
  MathParser::EvalDivVector(threadHandles_[handle][aThread],divergence);
}

void MathParserOMP::EvalMatrix( HandleType handle, Matrix<Double>& matrix,
             UInt numRows, UInt numCols ){
  UInt aThread = omp_get_thread_num();
  MathParser::EvalMatrix(threadHandles_[handle][aThread],matrix,numRows,numCols);
}

void MathParserOMP::Dump( std::ostream& os ){
  os << "=================================================================" << std::endl;
  os << "  SMP version of mathParser expressions copied to each thread    " << std::endl;
  os << "=================================================================" << std::endl;
  MathParser::Dump(os);
}

void MathParserOMP::SetValue( HandleType handle,
           const std::string &varName,
           Double val ){
  if(omp_get_num_threads() == 1){
    //assume call from serial region
    for(UInt aT=0;aT<numThreads_;aT++){
      MathParser::SetValue(threadHandles_[handle][aT],varName,val);
    }
  }else{
    //assume from parallel region
    UInt aThread = omp_get_thread_num();
    MathParser::SetValue(threadHandles_[handle][aThread],varName,val);
  }
}

void MathParserOMP::RegisterExternalVar( HandleType handle,
                      const std::string& varName,
                      Double * ptVar ){
  if(omp_get_num_threads() == 1){
    //assume call from serial region
    for(UInt aT=0;aT<numThreads_;aT++){
      MathParser::RegisterExternalVar(threadHandles_[handle][aT],varName,ptVar);
    }
  }else{
    //assume from parallel region
    UInt aThread = omp_get_thread_num();
    MathParser::RegisterExternalVar(threadHandles_[handle][aThread],varName,ptVar);
  }
}

void MathParserOMP::SetCoordinates( HandleType handle,
                         const CoordSystem &coosy,
                         const Vector<Double> &globCoord ){
  if(omp_get_num_threads() == 1){
    //assume call from serial region
    for(UInt aT=0;aT<numThreads_;aT++){
      MathParser::SetCoordinates(threadHandles_[handle][aT],coosy,globCoord);
    }
  }else{
    //assume from parallel region
    UInt aThread = omp_get_thread_num();
    MathParser::SetCoordinates(threadHandles_[handle][aThread],coosy,globCoord);
  }
}

boost::signals2::connection
MathParserOMP::AddExpChangeCallBack( const MathParserSignal::slot_function_type
                          &subscriber,
                          HandleType handle ){
  //only master
  if(omp_get_num_threads() != 1)
    EXCEPTION("Trying to add mathParser callback from parallel region. Not Supported!");

  return MathParser::AddExpChangeCallBack(subscriber,threadHandles_[handle][0]);

}

}

