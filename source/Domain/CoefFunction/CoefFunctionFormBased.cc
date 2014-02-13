#include "CoefFunctionFormBased.hh"

#include "Forms/BiLinForms/BDBInt.hh"

namespace CoupledField  {

// ==========================================================================
//  FORM BASED COEFFICIENT FUNCTION
// ==========================================================================

CoefFunctionFormBased::CoefFunctionFormBased( ) : CoefFunction(){

  
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
  //This is not very nice. Because another operator can turn
  //this coefficient function into another type.
  //Inconsistencies are the result, messing up other CoefFunctions
  //THIS HAS TO BE CHANGED!!!!
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
void CoefFunctionBOp<TYPE>::GetTensor(Matrix<TYPE>& coefMat,
                                          const LocPointMapped& lpm ) {
  EXCEPTION("CoefFunctionBOp<TYPE>::GetTensor not implemented");
  //Vector<TYPE> tmpVec;
  //this->GetVector(tmpVec,lpm);
  ////now put into matrix can be 2x2 or 3x3
  ////either 3 or 6 vector components
  //if(tmpVec.GetSize() == 3){
  //  coefMat.Resize(2,2);
  //  coefMat.Init();
  //  for(UInt i=0;i<2;i++)
  //    coefMat[i][i] = tmpVec[i];
  //  coefMat[0][1] = tmpVec[3];
  //  coefMat[1][0] = tmpVec[3];
  //
  //}else if (tmpVec.GetSize() == 6){
  //  coefMat.Resize(3,3);
  //  coefMat.Init();
  //  for(UInt i=0;i<3;i++)
  //    coefMat[i][i] = tmpVec[i];
  //
  //  coefMat[0][1] = tmpVec[3];
  //  coefMat[0][2] = tmpVec[4];
  //  coefMat[1][2] = tmpVec[5];
  //  coefMat[1][0] = tmpVec[3];
  //  coefMat[2][0] = tmpVec[4];
  //  coefMat[2][1] = tmpVec[5];
  //}else{
  //  EXCEPTION("CoefFunctionBOp<TYPE>::GetTensor: Got an unknown Voigt vector dimension");
  //}
}

template<class TYPE>
void CoefFunctionBOp<TYPE>::GetVector(Vector<TYPE>& coefVec,
                                      const LocPointMapped& lpm ) {
  this->feFct_->GetElemSolution( elemSol_, lpm.ptEl);
  BaseFE* ptFe = feSpace_->GetFe( lpm.ptEl->elemNum );
  if(this->bOps_.find(lpm.ptEl->regionId) != this->bOps_.end()){
    this->bOps_[lpm.ptEl->regionId]->CalcOpMat(bMat_, lpm, ptFe);
    coefVec = bMat_* (elemSol_);
    coefVec *= factor_;
  } else{
    coefVec.Resize(this->res_->dofNames.GetSize());
    coefVec.Init();
  }
}

template<class TYPE>
void CoefFunctionBOp<TYPE>::GetScalar(TYPE& coefScal,
                                          const LocPointMapped& lpm ) {
  Vector<TYPE> tmpVec;
  this->GetVector(tmpVec,lpm);
  coefScal = tmpVec[0];
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

template<class TYPE, bool TRANS>
void CoefFunctionFlux<TYPE,TRANS>::GetTensor(Matrix<TYPE>& coefMat,
                                          const LocPointMapped& lpm ) {
  EXCEPTION("CoefFunctionFlux<TYPE>::GetTensor not implemented");
  //Vector<TYPE> tmpVec;
  //this->GetVector(tmpVec,lpm);
  ////now put into matrix can be 2x2 or 3x3
  ////either 3 or 6 vector components
  //if(tmpVec.GetSize() == 3){
  //  coefMat.Resize(2,2);
  //  coefMat.Init();
  //  for(UInt i=0;i<2;i++)
  //    coefMat[i][i] = tmpVec[i];
  //  coefMat[0][1] = tmpVec[3];
  //  coefMat[1][0] = tmpVec[3];
  //
  //}else if (tmpVec.GetSize() == 6){
  //  coefMat.Resize(3,3);
  //  coefMat.Init();
  //  for(UInt i=0;i<3;i++)
  //    coefMat[i][i] = tmpVec[i];
  //
  //  coefMat[0][1] = tmpVec[3];
  //  coefMat[0][2] = tmpVec[4];
  //  coefMat[1][2] = tmpVec[5];
  //  coefMat[1][0] = tmpVec[3];
  //  coefMat[2][0] = tmpVec[4];
  //  coefMat[2][1] = tmpVec[5];
  //}else{
  //  EXCEPTION("CoefFunctionBOp<TYPE>::GetTensor: Got an unknown Voigt vector dimension");
  //}
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


template<class TYPE,bool TRANS>
void CoefFunctionFlux<TYPE,TRANS>::GetScalar(TYPE& coefScal,
                                          const LocPointMapped& lpm ) {
  Vector<TYPE> tmpVec;
  this->GetVector(tmpVec,lpm);
  coefScal = tmpVec[0];
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
  this->dimType_ = CoefFunction::SCALAR;

 // std::cout << "DimType: " << CoefFunction::SCALAR << " << std::endl;

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
  coefScal = (temp * Conj(elemSol_) ) * factor_;
  
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
