#include "CoefFunctionFormBased.hh"


#include "DataInOut/Logging/LogConfigurator.hh"
#include "General/Environment.hh"
#include "PDE/StdPDE.hh"
#include "PDE/LatticeBoltzmannPDE.hh"
#include "Forms/BiLinForms/BDBInt.hh"


namespace CoupledField
{

DECLARE_LOG(cff)
DEFINE_LOG(cff, "coefFunctionForm")

// ==========================================================================
//  FORM BASED COEFFICIENT FUNCTION
// ==========================================================================

CoefFunctionFormBased::CoefFunctionFormBased( ) : CoefFunction(){

  
  // All coefficient functions are based on the solution vector, so 
  // we have a general SOLUTION dependency.
  dependType_ = CoefFunction::SOLUTION;
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
template<class TYPE> void CoefFunctionBdBKernel<TYPE>::GetScalar( TYPE& coefScal, const LocPointMapped& lpm )
{
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
  coefScal = temp.Inner(elemSol_) * factor_;
}
template<class TYPE> std::string CoefFunctionBdBKernel<TYPE>::ToString() const
{
  std::stringstream out;
  out << "CoefFunctionBdBKernel\n";
  out << "Result: " << SolutionTypeEnum.ToString(feFct_->GetResultInfo()->resultType );
  return out.str();
}


template<class TYPE> CoefFunctionDyadicStrain<TYPE>::CoefFunctionDyadicStrain(shared_ptr<BaseFeFunction> feFct) : CoefFunctionFormBased()
{
  LOG_DBG2(cff) << "CFDS:CFDS #forms=" << this->forms_.size();

  feFct_ = dynamic_pointer_cast<FeFunction<TYPE> >(feFct);

  isComplex_ =  std::tr1::is_same<TYPE,Complex>::value;

  // set inherited attributes
  this->dimType_ = CoefFunction::TENSOR;
}

template<class TYPE> CoefFunctionDyadicStrain<TYPE>::~CoefFunctionDyadicStrain() { }

template<class TYPE> void CoefFunctionDyadicStrain<TYPE>::GetTensorSize(unsigned int& numRows, unsigned int& numCols ) const
{
  LOG_DBG2(cff) << "CFDS:GTS #forms=" << this->forms_.size();
  numRows = domain->GetGrid()->GetDim() == 2 ? 3 : 6;
  numCols = numRows;
}

template<class TYPE> void CoefFunctionDyadicStrain<TYPE>::GetTensor(Matrix<TYPE>& tensor, const LocPointMapped& lpm)
{
  LOG_DBG2(cff) << "CFDS:GT #forms=" << this->forms_.size();
  BaseBDBInt* form = this->forms_[lpm.ptEl->regionId];

  // LOG_DBG2(cff) << "CFDS:GT lpm=" << lpm.ptEl->elemNum << " form=" << form->GetName() << " B=" << form->GetBOp()->GetName();

  // Matrix<TYPE> B;
  // shared_ptr<FeSpace> mySpace = domain->GetStdPDE("mechanic")->GetFeFunction(MECH_DISPLACEMENT)->GetFeSpace();
  // BaseFE* ptFe = mySpace->GetFe( lpm.ptEl->elemNum);
  // form->GetBOp()->CalcOpMat(B, lpm, ptFe);
  // LOG_DBG2(cff) << "CFDS:GT B=" << B.ToString(2);


  // the dyadic product (B^T u) (u^T B) meant for integration, so not simply strain^2
  this->feFct_->GetElemSolution(elemSol_, lpm.ptEl);
  LOG_DBG2(cff) << "CFDS:GT u=" << elemSol_.ToString(2);

  Vector<TYPE> tmp;
  form->ApplyBMat(tmp, elemSol_, lpm);
  // LOG_DBG2(cff) << "CFDS:GT B*u=" << tmp.ToString(2);

  // Vector<TYPE> manual;
  // manual = B * elemSol_;
  // LOG_DBG2(cff) << "CFDS:GT manual=" << manual.ToString(2);

  tensor.DyadicMult(tmp, tmp.Conj());
  LOG_DBG2(cff) << "CFDS:GT -> " << tensor.ToString(2);
}

template<class TYPE> std::string CoefFunctionDyadicStrain<TYPE>::ToString() const
{
  std::stringstream out;
  out << "CoefFunctionDyadicStrain\n";
  out << "Result: " << SolutionTypeEnum.ToString(feFct_->GetResultInfo()->resultType );
  return out.str();
}

// --------------------------- ddd
template<class TYPE> CoefFunctionQuadSol<TYPE>::CoefFunctionQuadSol(shared_ptr<BaseFeFunction> feFct) : CoefFunctionFormBased()
{
  feFct_ = dynamic_pointer_cast<FeFunction<TYPE> >(feFct);

  isComplex_ =  std::tr1::is_same<TYPE,Complex>::value;

  // set inherited attributes
  this->dimType_ = CoefFunction::SCALAR;
}

template<class TYPE> CoefFunctionQuadSol<TYPE>::~CoefFunctionQuadSol() { }


template<class TYPE> void CoefFunctionQuadSol<TYPE>::GetScalar(TYPE& coefScal, const LocPointMapped& lpm)
{
  // the dyadic product (B^T u) (u^T B) meant for integration, so not simply strain^2
  this->feFct_->GetElemSolution(elemSol_, lpm.ptEl);
  LOG_DBG2(cff) << "CFSQ:GS u=" << elemSol_.ToString(2);

  coefScal = elemSol_.Inner();
  LOG_DBG2(cff) << "CFSQ:GS -> " << coefScal;
}

template<class TYPE> std::string CoefFunctionQuadSol<TYPE>::ToString() const
{
  std::stringstream out;
  out << "CoefFunctionQuadSol\n";
  out << "Result: " << SolutionTypeEnum.ToString(feFct_->GetResultInfo()->resultType );
  return out.str();
}

// -----------------------

// ---------------------- Customized coeff function for optimization with lattice Boltzmann(LBM) -------------------
template<class TYPE> CoefFunctionLBM<TYPE>::CoefFunctionLBM(LatticeBoltzmannPDE* lbm, shared_ptr<BaseFeFunction> feFct,shared_ptr<ResultInfo> resInfo) : CoefFunctionFormBased()
{

//  feFct_ = dynamic_pointer_cast<FeFunction<TYPE> >(feFct);

  isComplex_ =  false;

  resType_ = resInfo->resultType;

  lbm_ = lbm;

  // set inherited attributes
  switch(resInfo->entryType ) {
    case ResultInfo::SCALAR:
      this->dimType_ = SCALAR;
      break;
    case ResultInfo::VECTOR:
      this->dimType_ = VECTOR;
      break;
    default:
      EXCEPTION("Entry type of result can only be scalar or vector!")
  }
}

template<class TYPE> CoefFunctionLBM<TYPE>::~CoefFunctionLBM() { }

template<class TYPE> void CoefFunctionLBM<TYPE>::GetScalar(TYPE& coefScal, const LocPointMapped& lpm)
{
  assert(this->dimType_ == SCALAR);
  lbm_->ExtractIntermediateSolution();
  unsigned int idx = lbm_->GetLbmId(lpm.ptEl->elemNum);
  switch (resType_)
  {
    case LBM_DENSITY:
      coefScal = lbm_->CalcLBMDensity(idx); // all Calc... functions in LbmPDE need element ids already converted to LBM world
      break;
    case LBM_PRESSURE:
      coefScal = lbm_->CalcPressure(idx);
      break;
    default:
      EXCEPTION("LBM optimization has only pressure and density as scalar solution.")
  }
}

template<class TYPE> void CoefFunctionLBM<TYPE>::GetVector(Vector<TYPE>& vec, const LocPointMapped& lpm)
{
    assert(this->dimType_ == VECTOR);
    lbm_->ExtractIntermediateSolution();
    unsigned int idx = lbm_->GetLbmId(lpm.ptEl->elemNum);
    switch (resType_)
    {
      case LBM_VELOCITY:
        vec = lbm_->CalcVelocities(idx);
        break;
      case LBM_PROBABILITY_DISTRIBUTION:
        vec = lbm_->ExtractDistribution(idx);
        break;
      default:
        EXCEPTION("LBM optimization has only velocity and distribution function values as vector solution.")
    }

}

template<class TYPE> std::string CoefFunctionLBM<TYPE>::ToString() const
{
    return "We don't need this.";
}

template<class TYPE> CoefFunctionStiffness<TYPE>::CoefFunctionStiffness(shared_ptr<BaseFeFunction> feFct, DesignMaterial::Notation notation) : CoefFunctionFormBased()
{
  LOG_DBG2(cff) << "CFS:CFS #forms=" << this->forms_.size();

  feFct_ = dynamic_pointer_cast<FeFunction<TYPE> >(feFct);

  isComplex_ =  std::tr1::is_same<TYPE,Complex>::value;

  notation_ = notation;

  // set inherited attributes
  this->dimType_ = CoefFunction::TENSOR;
}

template<class TYPE> CoefFunctionStiffness<TYPE>::~CoefFunctionStiffness() { }

template<class TYPE> unsigned int CoefFunctionStiffness<TYPE>::GetVecSize() const
{
  LOG_DBG2(cff) << "CFS:GVS #forms=" << this->forms_.size();

  unsigned int numRows = 0;
  unsigned int numCols = 0;

  this->forms_.begin()->second->GetCoef()->GetTensorSize(numRows, numCols);

  assert(numRows == (domain->GetGrid()->GetDim() == 2 ? 3 : 6) );
  assert(numCols == numRows);

  if(numRows == 3)
    return 6;
  else
    return 21;
}

template<class TYPE> void CoefFunctionStiffness<TYPE>::GetVector(Vector<TYPE>& vec, const LocPointMapped& lpm)
{
  BaseBDBInt* form = this->forms_[lpm.ptEl->regionId];

  LOG_DBG(cff) << "CFS:GTV #forms=" << this->forms_.size() << " form=" << form->GetName() << " analytic=" << form->GetCoef()->IsAnalytic() << " notation=" << notation_;

  // the tensor here is always Voigt notation as we have Voigt on the simulation side of CFS (and Hill-Mandel on the optimization side)
  Matrix<TYPE> tensor;
  form->GetCoef()->GetTensor(tensor, lpm);
  if(notation_ == DesignMaterial::HILL_MANDEL)
    tensor.VoigtToHillMandel();

  tensor.ConvertToVec_UpperTriangular(vec);
  LOG_DBG3(cff) << "CFS:GV tensor=" << tensor.ToString(2);
  LOG_DBG2(cff) << "CFS:GV -> " << vec.ToString(2);

}

template<class TYPE> std::string CoefFunctionStiffness<TYPE>::ToString() const
{
  std::stringstream out;
  out << "CoefFunctionStiffness res " << SolutionTypeEnum.ToString(feFct_->GetResultInfo()->resultType );
  return out.str();
}

template<class TYPE> void CoefFunctionStiffness<TYPE>::GetTensorSize(unsigned int& numRows, unsigned int& numCols ) const
{
  LOG_DBG2(cff) << "CFS:GTS #forms=" << this->forms_.size();

  this->forms_.begin()->second->GetCoef()->GetTensorSize(numRows, numCols);

  assert(numRows == (domain->GetGrid()->GetDim() == 2 ? 3 : 6) );
  assert(numCols == numRows);
}

template<class TYPE> void CoefFunctionStiffness<TYPE>::GetTensor(Matrix<TYPE>& tensor, const LocPointMapped& lpm)
{
  BaseBDBInt* form = this->forms_[lpm.ptEl->regionId];

  LOG_DBG(cff) << "CFS:GT #forms=" << this->forms_.size() << " form=" << form->GetName() << " analytic=" << form->GetCoef()->IsAnalytic() << " notation=" << notation_;

  // the tensor here is always Voigt notation as we have Voigt on the simulation side of CFS (and Hill-Mandel on the optimization side)
  form->GetCoef()->GetTensor(tensor, lpm);
  if(notation_ == DesignMaterial::HILL_MANDEL)
    tensor.VoigtToHillMandel();
  LOG_DBG2(cff) << "CFS:GT -> " << tensor.ToString(2);
}



template class CoefFunctionBdBKernel<double>;
template class CoefFunctionBdBKernel<Complex>;

template class CoefFunctionDyadicStrain<double>;
template class CoefFunctionDyadicStrain<Complex>;

template class CoefFunctionQuadSol<double>;
template class CoefFunctionQuadSol<Complex>;

template class CoefFunctionStiffness<double>;
template class CoefFunctionStiffness<Complex>;

template class CoefFunctionLBM<double>;

   
} // end of namespace
