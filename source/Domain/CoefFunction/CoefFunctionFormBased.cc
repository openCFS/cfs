#include "CoefFunctionFormBased.hh"

#include "Forms/BiLinForms/BDBInt.hh"

namespace CoupledField  {

// ==========================================================================
//  FORM BASED COEFFICIENT FUNCTION
// ==========================================================================

CoefFunctionFormBased::CoefFunctionFormBased( ) {

  
  // set inherited attributes
  dependType_ = CoefFunction::GENERAL;
  isAnalytic_ = false;
}

CoefFunctionFormBased::~CoefFunctionFormBased() {
  
}

void CoefFunctionFormBased::AddIntegrator( BaseBDBInt* form,  
                                           RegionIdType region ) {
  
  // check if region has already an integrator assigned
  if( forms_.find(region) != forms_.end() ) {
    EXCEPTION( "Multiply defined region");
  }
  
  forms_[region] = form;
}    
  
// ==========================================================================
//  COEFFICIENT FUNCTION BASED ON DIFFERENTIAL OPERATOR
// ==========================================================================
template<class TYPE> CoefFunctionBOp<TYPE>::
CoefFunctionBOp( shared_ptr<BaseFeFunction> feFct,
                 shared_ptr<ResultInfo> info,
                 TYPE factor ) 
                 : CoefFunctionFormBased() {
  
  feFct_ = dynamic_pointer_cast<FeFunction<TYPE> >(feFct);
  res_ = info;
  factor_ = factor;
  feSpace_ = feFct->GetFeSpace();
  
  isComplex_ =  std::tr1::is_same<TYPE,Complex>::value;
  
  // set inherited attributes
  dimType_ = CoefFunction::VECTOR;
}


template<class TYPE>  void CoefFunctionBOp<TYPE>::
AddIntegrator( BaseBDBInt* form,  
               RegionIdType region ) {

  // call method in base class
  CoefFunctionFormBased::AddIntegrator( form, region );
  
  // store B-operator as well
  AddBOperator( form->GetBOp(), region );
}


template<class TYPE>  void CoefFunctionBOp<TYPE>::
AddBOperator( BaseBOperator* bOp,
              RegionIdType region ) {
  
  // check for compatibility of dimension
  if( bOp->GetDimDMat() != res_->dofNames.GetSize() ) {
    EXCEPTION( "All B-operators must have the same vector size");
  }
  
  // check if region has already an integrator assigned
  if( bOps_.find(region) != bOps_.end() ) {
    EXCEPTION( "Multiply defined region");
  }
  
  bOps_[region] = bOp;

}

template<class TYPE>
CoefFunctionBOp<TYPE>::~CoefFunctionBOp() {

}

template<class TYPE>
void CoefFunctionBOp<TYPE>::GetVector(Vector<TYPE>& coefVec,
                                      const LocPointMapped& lpm ) {
  this->feFct_->GetElemSolution( elemSol_, lpm.ptEl);
  BaseFE* ptFe = feSpace_->GetFe( lpm.ptEl->elemNum );
  this->bOps_[lpm.ptEl->regionId]->CalcOpMat(bMat_, lpm, ptFe);
  coefVec = bMat_* (elemSol_);
  coefVec *= factor_;
}

template<class TYPE>
UInt CoefFunctionBOp<TYPE>::GetVecSize() const {
  return this->res_->dofNames.GetSize();

}

template<class TYPE>
std::string CoefFunctionBOp<TYPE>::ToString() const {
  std::stringstream out;
  out << "CoefFunctionBOp\n";
  out << "\tFeFunction: " << 
      SolutionTypeEnum.ToString(feFct_->GetResultInfo()->resultType);
  return out.str();

}

template class CoefFunctionBOp<Double>;
template class CoefFunctionBOp<Complex>;

// ==========================================================================
//  COEFFICIENT FUNCTION BASED ON FLUX COMPUTATION OF A BDB-FORM
// ==========================================================================
template<class TYPE, bool TRANS> CoefFunctionFlux<TYPE,TRANS>::
CoefFunctionFlux( shared_ptr<BaseFeFunction> feFct,
                  shared_ptr<ResultInfo> info,
                  TYPE factor )
                  :CoefFunctionFormBased() {
  feFct_ = dynamic_pointer_cast<FeFunction<TYPE> >(feFct);
  res_ = info;
  factor_ = factor;
  isComplex_ =  std::tr1::is_same<TYPE,Complex>::value;

  // set inherited attributes
  dimType_ = CoefFunction::VECTOR;
}

template<class TYPE, bool TRANS> CoefFunctionFlux<TYPE,TRANS>::
~CoefFunctionFlux() {


}
template<class TYPE, bool TRANS> void CoefFunctionFlux<TYPE,TRANS>::
AddIntegrator( BaseBDBInt* form,  
               RegionIdType region ) {

  // check if region has already an integrator assigned
  if( forms_.find(region) != forms_.end() ) {
    EXCEPTION( "Multiply defined region");
  }

  
  
  // check for compatibility of dimension
  PtrCoefFct coef = form->GetCoef();
  if( coef->GetDimType() == CoefFunction::TENSOR ) {
    UInt nRows, nCols;
    coef->GetTensorSize(nRows, nCols);
    UInt dMatSize = TRANS ? nCols : nRows;
    if( dMatSize != res_->dofNames.GetSize() ) {
      EXCEPTION( "All B-operators must have the same vector size");
    } 
  }
  forms_[region] = form;

}

template<class TYPE, bool TRANS> void CoefFunctionFlux<TYPE,TRANS>::
GetVector(Vector<TYPE>& coefVec,
          const LocPointMapped& lpm ) {
  
  this->feFct_->GetElemSolution( elemSol, lpm.ptEl);
  BaseBDBInt* bdb = this->forms_[lpm.ptEl->regionId];
  if( !bdb) {
    coefVec.Resize(res_->dofNames.GetSize());
    coefVec.Init();
    return;
  }
  
  // switch depending on if 
  if( TRANS ) {
    bdb->ApplydATransMat(coefVec, elemSol, lpm );
  } else {
    bdb->ApplydBMat(coefVec, elemSol, lpm );
  }
  coefVec *= factor_;
}

template<class TYPE, bool TRANS> UInt CoefFunctionFlux<TYPE,TRANS>::
GetVecSize() const{
     return res_->dofNames.GetSize();
}

template<class TYPE, bool TRANS> std::string CoefFunctionFlux<TYPE,TRANS>::
ToString() const {
  std::stringstream out;
  out << "CoefFunctionFlux\n";
  out << "ApplyTransposed: " << TRANS << std::endl;
  out << "Result: " << 
      SolutionTypeEnum.ToString(feFct_->GetResultInfo()->resultType );
  return out.str();
}

template class CoefFunctionFlux<Double,false>;
template class CoefFunctionFlux<Complex,false>;
template class CoefFunctionFlux<Double,true>;
template class CoefFunctionFlux<Complex,true>;

// ==========================================================================
//  COEFFICIENT FUNCTION BASED ON KERNEL OF BDB-INTEGRATOR
// ==========================================================================
template<class TYPE> CoefFunctionBdBKernel<TYPE>::
CoefFunctionBdBKernel(shared_ptr<BaseFeFunction> feFct,
                      TYPE factor ) 
: CoefFunctionFormBased() {
  feFct_ = dynamic_pointer_cast<FeFunction<TYPE> >(feFct);
  factor_ = factor;
  
  isComplex_ =  std::tr1::is_same<TYPE,Complex>::value;

  // set inherited attributes
  dimType_ = CoefFunction::SCALAR;
  
}

template<class TYPE> CoefFunctionBdBKernel<TYPE>::
  ~CoefFunctionBdBKernel() {
}
template<class TYPE> void CoefFunctionBdBKernel<TYPE>::
GetScalar( TYPE& coefScal,
           const LocPointMapped& lpm ) {
  
  
  // energy density is factor * elemSol^T * kernel * elemSol
  this->feFct_->GetElemSolution( elemSol_, lpm.ptEl);
  Vector<TYPE> temp;
  if( !this->forms_[lpm.ptEl->regionId]->IsComplex() ) {
    this->forms_[lpm.ptEl->regionId]->CalcKernel(kernelR_, lpm);
    temp = kernelR_ * elemSol_;
  } else {
    this->forms_[lpm.ptEl->regionId]->CalcKernel(kernel_, lpm);
    temp = kernel_ * elemSol_;
  }
  coefScal = (temp * elemSol_) * factor_;
  
}
template<class TYPE> std::string CoefFunctionBdBKernel<TYPE>::
ToString() const {
  std::stringstream out;
  out << "CoefFunctionBdBKernel\n";
  out << "Result: " << 
      SolutionTypeEnum.ToString(feFct_->GetResultInfo()->resultType );
  return out.str();
}

template class CoefFunctionBdBKernel<Double>;
template class CoefFunctionBdBKernel<Complex>;


   
} // end of namespace
