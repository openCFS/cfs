// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "result.hh"
#include "Domain/resultInfo.hh"
#include <sstream>

namespace CoupledField {

  
  BaseResult::BaseResult() {


  }

  BaseResult::~BaseResult() {

  }

  std::string BaseResult::ToString() const
  {
    std::ostringstream os;
    os << "entity list: " << entities_->GetName()
       //<< " values: " << (const_cast<const CFSVector*>(GetCFSVector())->GetSize()
       << " result info: " << resultDof_->ToString();
    return os.str();
  }

  void BaseResult::Dump(StdVector<shared_ptr<BaseResult> >& resultList)
  {
    for(unsigned int i = 0; i < resultList.GetSize(); i++)
      std::cout << resultList[i]->ToString() << std::endl;
  }

  template<class TYPE>
  Result<TYPE>::Result() {
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
