#include "BaseResults.hh"
#include "ResultInfo.hh"

#include <sstream>

#include "DataInOut/ParamHandling/ParamNode.hh"


namespace CoupledField {


  BaseResult::BaseResult() {

  }

  BaseResult::~BaseResult() {
  }

  std::string BaseResult::ToString() const
  {
    std::ostringstream os;
    os << "entity list: " << entities_->GetName()
       //<< " values: " << (const_cast<const SingleVector*>(GetSingleVector())->GetSize()
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

  template<class TYPE>
  void Result<TYPE>::Init()
  {
    values_.Resize(entities_->GetSize());
    values_.Init(0);
  }

 // explicit template instantiation for GCC compiler
#if defined(__GNUC__)
  template class Result<Double>;
  template class Result<Complex>;
#endif

}
