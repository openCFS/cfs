
#include "Utils/mathParser/mathParser.hh"
#include "CoefFunctionCache.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "General/Environment.hh"
#include "PDE/StdPDE.hh"
#include "PDE/LatticeBoltzmannPDE.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include <boost/bind/bind.hpp>

namespace CoupledField
{

DEFINE_LOG(cfc, "coefFunctionCache")

// ==========================================================================
//  COEFFICIENT FUNCTION CACHE
// ==========================================================================
// Example for using Coefficient Function Cache:
// See MechPDE::DefinePostProcResults() MechStress, MechPrincipalStress
template<class TYPE> CoefFunctionCache<TYPE>::CoefFunctionCache(
		shared_ptr<BaseFeFunction> feFct,
		shared_ptr<ResultInfo> info,
		shared_ptr<CoefFunction> cached_coefFct
		)
        :CoefFunction() {

  cached_coefFct_ = cached_coefFct;
  feFct_ = dynamic_pointer_cast<FeFunction<TYPE> >(feFct);
  res_ = info;
  dimType_ = cached_coefFct_->GetDimType();
  isComplex_ =  std::is_same<TYPE,Complex>::value;

  mp_ = domain->GetMathParser();
  mHandleTime_ = mp_->GetNewHandle(true);
  mp_->SetExpr(mHandleTime_,"t"); // current step
  //This is called automatically at the beginning of each step
  mp_->AddExpChangeCallBack( boost::bind(&CoefFunctionCache::EmptyCache, this ), mHandleTime_ );
}

template<class TYPE> UInt CoefFunctionCache<TYPE>::GetVecSize() const{
  return cached_coefFct_->GetVecSize();
}

template<class TYPE> CoefFunctionCache<TYPE>::~CoefFunctionCache() {
  mp_->ReleaseHandle(mHandleTime_);
}

template<class TYPE> void CoefFunctionCache<TYPE>::EmptyCache() {
  cached_data_scal_.clear();
  cached_data_vec_.clear();
  cached_data_mat_.clear();
}

template<class TYPE> void CoefFunctionCache<TYPE>::GetVector(Vector<TYPE>& coefVec, const LocPointMapped& lpm) {

  // If coefs in element have already been calculated (determined by elemID) -> take values
  // else coefFct->GetVector

  it_vec_ = cached_data_vec_.find(lpm.ptEl->elemNum);

  if (it_vec_ != cached_data_vec_.end()) {
	//std::cout << "Element ID found:" << res_->resultName << std::endl;
	//Element ID found
	//Cached coefVec is second value in pair. For data structure see comment in else
    coefVec = cached_data_vec_[lpm.ptEl->elemNum].second;

    Vector<Double> actlp = lpm.lp.coord;

    //Distance saved LP to actual LP
    Double dist = actlp.NormL2(cached_data_vec_[lpm.ptEl->elemNum].first);

    //Give Warning if distance between actual LP and saved LP is > 1e-3
    if (dist > 0.001) {
      WARN("Cached Results and requested results of GetVector do not have the same local point. Local distance: " << dist);
    }
  }
  else {
	//Element ID not found
	//std::cout << "Element ID not found: " << res_->resultName << std::endl;
	cached_coefFct_->GetVector(coefVec, lpm);
	//Map structure
	//Key: element number
	//Value: pair consisting of
	// - local coordinates
	// - coefficient vector
	cached_data_vec_.insert(std::make_pair(lpm.ptEl->elemNum, std::make_pair(lpm.lp.coord, coefVec)));
  }

}

template<class TYPE> void CoefFunctionCache<TYPE>::GetTensor(Matrix<TYPE>& coefMat, const LocPointMapped& lpm) {
  //For a inline documentation of the code, see GetVector()
  it_mat_ = cached_data_mat_.find(lpm.ptEl->elemNum);

  if (it_mat_ != cached_data_mat_.end()) {

    coefMat = cached_data_mat_[lpm.ptEl->elemNum].second;

    Vector<Double> actlp = lpm.lp.coord;
    Double dist = actlp.NormL2(cached_data_mat_[lpm.ptEl->elemNum].first);

    if (dist > 0.001) {
      WARN("Cached Results and requested results of GetTensor do not have the same local point. Local distance: " << dist);
    }
  }
  else {
    cached_coefFct_->GetTensor(coefMat, lpm);
    cached_data_mat_.insert(std::make_pair(lpm.ptEl->elemNum, std::make_pair(lpm.lp.coord, coefMat)));
  }
}

template<class TYPE> void CoefFunctionCache<TYPE>::GetScalar(TYPE& coefScal, const LocPointMapped& lpm ) {
  //For a inline documentation of the code, see GetVector()
  it_scal_ = cached_data_scal_.find(lpm.ptEl->elemNum);

  if (it_scal_ != cached_data_scal_.end()) {

    coefScal = cached_data_scal_[lpm.ptEl->elemNum].second;

    Vector<Double> actlp = lpm.lp.coord;
    Double dist = actlp.NormL2(cached_data_scal_[lpm.ptEl->elemNum].first);

    if (dist > 0.001) {
      WARN("Cached Results and requested results of GetScalar do not have the same local point. Local distance: " << dist);
    }
  }
  else {
    cached_coefFct_->GetScalar(coefScal, lpm);
    cached_data_scal_.insert(std::make_pair(lpm.ptEl->elemNum, std::make_pair(lpm.lp.coord, coefScal)));
  }
}

template<class TYPE> std::string CoefFunctionCache<TYPE>::ToString() const {
  std::stringstream out;
  out << "CoefFunctionCache\n";
  out << "Result: " <<
      SolutionTypeEnum.ToString(feFct_->GetResultInfo()->resultType );
  return out.str();
}

template class CoefFunctionCache<Double>;
template class CoefFunctionCache<Complex>;


} // end of namespace
