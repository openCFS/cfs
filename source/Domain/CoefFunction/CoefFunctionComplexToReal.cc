// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include <limits>

#include "Utils/mathParser/mathParser.hh"
#include "CoefFunctionComplexToReal.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"

#include "Domain/CoefFunction/CoefXpr.hh"

#include "FeBasis/FeSpace.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField
{
DEFINE_LOG(coeffctcomplextoreal, "coeffctcomplextoreal")


  // ==========================================================================
  //  COEFFICIENT FUNCTION CHANGE TYPE
  // ==========================================================================
  template<class T>
CoefFunctionComplexToReal<T>::CoefFunctionComplexToReal(PtrCoefFct coefToChange,
                                                  Grid* ptGrid):CoefFunction() {

    ptGrid_ = ptGrid;

    coefToChange_ = coefToChange;

    // Variables from CoefFunction base-class
    isAnalytic_ = false;
    isComplex_ = false;

    // take care, it's changed throughout the class!!
    dimType_ = coefToChange->GetDimType();
  }

  template<class T>
  CoefFunctionComplexToReal<T>::
  ~CoefFunctionComplexToReal() {}

  template<class T>
  void CoefFunctionComplexToReal<T>::
  GetVector(Vector<T>& coefVec, const LocPointMapped& lpm){

      if(dimType_ == 0){
        EXCEPTION("CoefFunctionComplexToReal::GetVector is called with a scalar coeffunction!");
      }

      Vector<Complex> tmp;
      coefToChange_->GetVector(tmp, lpm);
      coefVec = tmp.GetPart(Global::REAL);
  }

  template<class T>
  void CoefFunctionComplexToReal<T>::
  GetTensor(Matrix<T>& coefMat, const LocPointMapped& lpm){
    EXCEPTION("CoefFunctionComplexToReal::GetTensor NOT IMPLEMENTED YET")
  }

  template<class T>
  void CoefFunctionComplexToReal<T>::
  GetScalar(T& coefScal, const LocPointMapped& lpm){
    if(dimType_ != 0){
      EXCEPTION("CoefFunctionComplexToReal::GetScalar is called with a vector coeffunction!");
    }

    Complex tmp;
    coefToChange_->GetScalar(tmp, lpm);
    coefScal = tmp.real();
  }

  template<class T>
  UInt CoefFunctionComplexToReal<T>::
  GetVecSize() const{
    // spatial dimension (dimension of B vector)
    return ptGrid_->GetDim();
  }

  template<class T>
  std::string CoefFunctionComplexToReal<T>::
  ToString() const {
    return "CoefFunctionComplexToReal";
  }


  // Explicit template instantiation
  template class CoefFunctionComplexToReal<Double>;
  template class CoefFunctionComplexToReal<Complex>;

} // end of namespace
