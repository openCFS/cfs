// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef COEF_FUNCTION_CACHE_HH
#define COEF_FUNCTION_CACHE_HH


#include "CoefFunction.hh"
#include "FeBasis/BaseFE.hh"
#include "FeBasis/FeFunctions.hh"
#include <unordered_map>

namespace CoupledField  {

// ==========================================================================
//  COEFFICIENT FUNCTION CACHE
// ==========================================================================

//! This coefficient function caches values from another CoefFunction.
//
//! If GetScalar() or GetVector() or GetTensor() is called, a check is performed if
//! this has already been done during runtime for the element and the cached coefFct.
//! If not: Values are returned and saved in a std::unordered_map
//! Else: Cached Values are returned, no (repeated) evaluation done.
//! Tested with CoefFunctionEigen (principal strain, principal stress)

template <class TYPE>
class CoefFunctionCache : public CoefFunction {

public:
  //! Constructor
  CoefFunctionCache(shared_ptr<BaseFeFunction> feFct,
		            shared_ptr<ResultInfo> info,
		            shared_ptr<CoefFunction> cached_coefFct);

  //! Destructor
  virtual ~CoefFunctionCache();

  void EmptyCache();

  std::string GetName() const { return "CoefFunctionCache"; }

  //! \copydoc CoefFunction::GetVector
  virtual void GetVector(Vector<TYPE>& coefVec, const LocPointMapped& lpm);

  //! \copydoc CoefFunction::GetTensor
  virtual void GetTensor(Matrix<TYPE>& coefMat, const LocPointMapped& lpm);

  //! \copydoc CoefFunction::ToString
  virtual std::string ToString() const;

  //! \copydoc CoefFunction::GetVecSize
  virtual UInt GetVecSize() const;

  //! \copydoc CoefFunction::GetScalar
  virtual void GetScalar(TYPE& coefScal, const LocPointMapped& lpm );

protected:
  //! Coefficient function which should be cached
  shared_ptr<CoefFunction> cached_coefFct_;

  //! FeFunction containing the coefficients
  shared_ptr<FeFunction<TYPE> > feFct_;

  //! Result info object of cached result
  shared_ptr<ResultInfo> res_;

  //! Map which takes the cached data
  std::unordered_map<UInt, std::pair<Vector<Double> , TYPE > > cached_data_scal_;
  std::unordered_map<UInt, std::pair<Vector<Double> , Vector<TYPE> > > cached_data_vec_;
  std::unordered_map<UInt, std::pair<Vector<Double> , Matrix<TYPE> > > cached_data_mat_;

  //Iterator
  typename std::unordered_map<UInt, std::pair<Vector<Double>, TYPE > >::iterator it_scal_;
  typename std::unordered_map<UInt, std::pair<Vector<Double>, Vector<TYPE> > >::iterator it_vec_; //typename: the compiler likes it ;)
  typename std::unordered_map<UInt, std::pair<Vector<Double>, Matrix<TYPE> > >::iterator it_mat_;

  //std::unordered_map<UInt, std::pair<Vector<Double> , Vector<TYPE> > > cached_data_boost_;
  //typename std::unordered_map<UInt, std::pair<Vector<Double>, Vector<TYPE> > >::iterator it_boost_;

  MathParser* mp_ = nullptr;
  unsigned int mHandleTime_;

};
} //end of namespace

#endif
