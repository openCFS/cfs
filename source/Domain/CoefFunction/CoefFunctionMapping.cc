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



#include "CoefFunctionMapping.hh"

#include "boost/lexical_cast.hpp"
#include "Domain/Mesh/Grid.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"

namespace CoupledField{

template<typename T>
CoefFunctionMapping<T>::CoefFunctionMapping(PtrParamNode mapDef, PtrCoefFct matCoef,
                                    shared_ptr<EntityList> EntList,
                                    StdVector<RegionIdType> pdeDomains,
                                    CoefFunction::CoefDimType type):
  CoefFunctionPML<T>(mapDef, matCoef,EntList,pdeDomains,type){

}

template<typename T>
CoefFunctionMapping<T>::~CoefFunctionMapping(){

}


template<typename T>
void CoefFunctionMapping<T>::GetTensor(Matrix<Complex>& tensor,
                              const LocPointMapped& lpm ){
	  //this is diagonal tensor with the coefficients
	  tensor.Resize(this->dim_,this->dim_);
	  tensor.Init();
	  Double locThick=0.0;
	  Double position=0.0;
	  Double sos;
	  this->matCoef_->GetScalar(sos,lpm);
	  for(UInt i=0;i<this->dim_;++i){
	    this->GetThicknessAtPoint(locThick,position,lpm,i);
        Double z = position/locThick;
	    if(abs(locThick)>0.0){
          Complex fac(locThick*this->dampFunction_->ComputeFactor(z,sos),0.0);
	      tensor[i][i] = fac;
	    }else{
	      Complex one(1.0,0.0);
	      tensor[i][i] = one;
	    }
	  }
}

template<typename T>
void CoefFunctionMapping<T>::GetTensor(Matrix<Double>& tensor,
                              const LocPointMapped& lpm ){
  //this is diagonal tensor with the coefficients
  tensor.Resize(this->dim_,this->dim_);
  tensor.Init();
  Double locThick=0.0;
  Double position=0.0;
  Double sos;
  this->matCoef_->GetScalar(sos,lpm);
  for(UInt i=0;i<this->dim_;++i){
    this->GetThicknessAtPoint(locThick,position,lpm,i);
    Double z = position/locThick;
    if(abs(locThick)>0.0){
      tensor[i][i] = locThick*this->dampFunction_->ComputeFactor(z,sos) ;
    }else{
      tensor[i][i] = 1.0;
    }
  }
}

template<typename T>
void CoefFunctionMapping<T>::GetVector(Vector<Complex>& vec,
                              const LocPointMapped& lpm ){
	  //this is diagonal tensor with the coefficients
	vec.Resize(this->dim_);
	vec.Init();
	  Double locThick=0.0;
	  Double position=0.0;
	  Double sos;
	  this->matCoef_->GetScalar(sos,lpm);
	  for(UInt i=0;i<this->dim_;++i){
	    this->GetThicknessAtPoint(locThick,position,lpm,i);
        Double z = position/locThick;
	    if(abs(locThick)>0.0){
	    	Complex fac(locThick*this->dampFunction_->ComputeFactor(z,sos),0.0);
	    	vec[i] = fac;
	    }else{
	      Complex one(1.0,0.0);
	      vec[i] = one;
	    }
	  }

}

template<typename T>
void CoefFunctionMapping<T>::GetVector(Vector<Double>& vec,
                              const LocPointMapped& lpm ){
  //first loop over every entry and determine the factor
  vec.Resize(this->dim_,0.0);
  Double locThick=0.0;
  Double position=0.0;
  Double sos;
  this->matCoef_->GetScalar(sos,lpm);
  for(UInt i=0;i<this->dim_;++i){
    this->GetThicknessAtPoint(locThick,position,lpm,i);
    Double z = position/locThick;
    if(abs(locThick)>0.0){
      vec[i] = locThick*this->dampFunction_->ComputeFactor(z,sos);
    }else{
      vec[i] = 1.0;
    }
  }
}

template<typename T>
void CoefFunctionMapping<T>::GetScalar(Complex& val,
                              const LocPointMapped& lpm ){
  Double locThick=0.0;
  Double position=0.0;
  Complex one(1.0,0.0);
  val = one;
  Double sos;
  this->matCoef_->GetScalar(sos,lpm);
  for(UInt i=0;i<this->dim_;++i){
    this->GetThicknessAtPoint(locThick,position,lpm,i);
    Double z = position/locThick;
    if(abs(locThick)>0.0){
      Complex fac(locThick*this->dampFunction_->ComputeFactor(z,sos),0.0);
      val /= fac;
    }else{
      val *= one;
    }
  }
}

template<typename T>
void CoefFunctionMapping<T>::GetScalar(Double& val,
                              const LocPointMapped& lpm ){
  //computes 1/(map_x*map_y*map_z)
  Double locThick=0.0;
  Double position=0.0;
  val = 1.0;
  Double sos;
  this->matCoef_->GetScalar(sos,lpm);
  for(UInt i=0;i<this->dim_;++i){
    this->GetThicknessAtPoint(locThick,position,lpm,i);
    Double z = position/locThick;
    if(abs(locThick)>0.0){
      val /= (locThick*this->dampFunction_->ComputeFactor(z,sos));
    }else{
      val *= 1.0;
    }
  }
}


// Explicit template instantiation
template class CoefFunctionMapping<Double>;
template class CoefFunctionMapping<Complex>;
}
