// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionMapping.cc
 *       \brief    <Description>
 *
 *       \date     Mar 12, 2016
 *       \author   sschoder
 */
//================================================================================================

#include <def_expl_templ_inst.hh>

#include "CoefFunctionMapping.hh"

#include "boost/bind.hpp"
#include "boost/lexical_cast.hpp"
#include "Domain/Mesh/Grid.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"

namespace CoupledField{

template<typename T>
CoefFunctionMapping<T>::CoefFunctionMapping(PtrParamNode pmlDef, PtrCoefFct speedOfSound,
                                    shared_ptr<EntityList> EntList,
                                    StdVector<RegionIdType> pdeDomains,
                                    bool isVector):
  CoefFunctionPML<T>(pmlDef, speedOfSound,EntList,
                  pdeDomains,isVector){

}

template<typename T>
CoefFunctionMapping<T>::~CoefFunctionPML(){

}


template<typename T>
void CoefFunctionMapping<T>::GetTensor(Matrix<Complex>& tensor,
                              const LocPointMapped& lpm ){
  EXCEPTION("Not implemented: Complex mapping!");
}

template<typename T>
void CoefFunctionMapping<T>::GetTensor(Matrix<Double>& tensor,
                              const LocPointMapped& lpm ){
  //this is diagonal tensor with the coefficients
  tensor.Resize(this->dim_,this->dim_);
  tensor.Init();
  Double locThick=0.0;
  Double position=0.0;
  Double sos=1.0;
  for(UInt i=0;i<this->dim_;++i){
    this->GetThicknessAtPoint(locThick,position,lpm,i);
    if(abs(locThick)>0.0){
      tensor[i][i] = this->dampFunction_->ComputeFactor(position,locThick) * sos;
    }else{
      tensor[i][i] = 1.0;
    }
  }
}

template<typename T>
void CoefFunctionMapping<T>::GetVector(Vector<Complex>& vec,
                              const LocPointMapped& lpm ){
  EXCEPTION("Not implemented: Complex mapping!");
}

template<typename T>
void CoefFunctionMapping<T>::GetVector(Vector<Double>& vec,
                              const LocPointMapped& lpm ){

  //first loop over every entry and determine the factor
  vec.Resize(this->dim_,0.0);
  Double locThick=0.0;
  Double position=0.0;
  Double sos = 1.0;
  for(UInt i=0;i<this->dim_;++i){
    this->GetThicknessAtPoint(locThick,position,lpm,i);
    if(abs(locThick)>0.0){
      vec[i] = this->dampFunction_->ComputeFactor(position,locThick) * sos;
    }else{
      vec[i] = 1.0;
    }
  }
}

template<typename T>
void CoefFunctionMapping<T>::GetScalar(Complex& val,
                              const LocPointMapped& lpm ){
  EXCEPTION("Not implemented: Complex mapping!");
}

template<typename T>
void CoefFunctionMapping<T>::GetScalar(Double& val,
                              const LocPointMapped& lpm ){

  //computes 1/(pml_x*pml_y*pml_z)
  //right now we ust return 1...
  EXCEPTION("GETSCALAR IS INVALID");
  val = 1.0;
  return;

  //UInt gridDim = entities_[0]->GetGrid()->GetDim();
  //Double locThick=0.0;
  //Double position=0.0;
  //val = 1.0;
  //for(UInt i=0;i<gridDim;++i){
  //  GetThicknessAtPoint(locThick,position,lpm,i);
  //  if(abs(locThick)>0.0){
  //    val *= dampFunction_->ComputeFactor(position,locThick);
  //  }else{
  //    val *= 1.0;
  //  }
  //}
}


// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class CoefFunctionMapping<Double>;
  template class CoefFunctionMapping<Complex>;
#endif
}
