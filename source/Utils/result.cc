// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "result.hh"

namespace CoupledField {

  
  BaseResult::BaseResult() {


  }

  BaseResult::~BaseResult() {

  }

  template<class TYPE>
  Result<TYPE>::Result() {
    ENTER_FCN("Result<TYPE>::Result()" );
  }
  
  template<class TYPE>
  Result<TYPE>::~Result() {
  }


 // explicit template instantiation for GCC compiler
#if defined(__GNUC__) 
  template class Result<Double>;
  template class Result<Complex>;
#endif 

}
