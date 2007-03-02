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
