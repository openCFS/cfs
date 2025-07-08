#include "CoefFunctionFormBased.hh"


#include "DataInOut/Logging/LogConfigurator.hh"
#include "General/Environment.hh"
#include "PDE/StdPDE.hh"
#include "PDE/LatticeBoltzmannPDE.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "MatVec/BLASLAPACKInterface.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/Excitation.hh"

namespace CoupledField
{
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
  // If openMP is used all the cloned Integrators must be deleted
  #ifdef USE_OPENMP
  if(omp_get_num_threads()!=1){
    std::ostringstream ostr; ostr << "Call from parallel region which is not safe!"; // Cannot except in destructor -> terminate would be generated anyways
    std::terminate();
  }
  for(UInt i=0;i<forms_.GetNumSlots();++i){
    auto myForms = forms_.Mine(i);
    for(auto form : myForms){
      delete form.second;
    }
    myForms.clear();
  }
  #endif
}

void CoefFunctionFormBased::AddIntegrator( BaseBDBInt* form,  
                                           RegionIdType region ) {
#ifdef USE_OPENMP
  if(omp_get_num_threads()!=1){
    EXCEPTION("Call from parallel region which is not safe!")
  }
#endif
  std::map<RegionIdType, BaseBDBInt* >& forms = forms_.Mine(0);

  // check if we need to skip the integrator assignment if e.g. the requested integrator name differs
  if( integratorName_.empty() || integratorName_==form->GetName() ) {
    // check if region has already an integrator assigned
    if( forms.find(region) != forms.end() ) {
      EXCEPTION( "Multiply defined region");
    }
    //here we check if we are from a single threaded region.
    //if so, we fill all slots. If not, we just fill the current thread
#ifdef USE_OPENMP
    for(UInt i=0;i<forms_.GetNumSlots();++i){
      forms_.Mine(i)[region] = form->Clone();
    }
#else
    forms_.Mine()[region] = form;
#endif
  } else {
    WARN("Skipped adding integrator " << form->GetName() << " to " << this->GetName());
  }
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
  
  isComplex_ =  std::is_same<TYPE,Complex>::value;
  
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
  AddBOperator( form->GetBOp(), region, form->GetName() );
}


template<class TYPE>  void CoefFunctionBOp<TYPE>::
AddBOperator( BaseBOperator* bOp,
              RegionIdType region,
              std::string integratorName ) {

  // check if we need to skip the integrator assignment if e.g. the requested integrator name differs
  if( integratorName_.empty() || integratorName_==integratorName ) {
    // check for compatibility of dimension
    if( bOp->GetDimDMat() != res_->dofNames.GetSize() ) {
      EXCEPTION( "Implementation error: All B-operators must have the same vector size");
    }
    
    // check if region has already an integrator assigned
    if( bOps_.find(region) != bOps_.end() ) {
      EXCEPTION( "Implementation error: Multiply defined region. Only one integrator can be assigned to a region.");
    }
    
    bOps_[region] = bOp;
  } else {
    WARN("Implementation warning: Skipped adding bOp " << bOp->GetName() << " from integrator " << integratorName << " to " << this->GetName());
  }

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
  Vector<TYPE>& elemSol = elemSol_.Mine();
  Matrix<TYPE>& bMat = bMat_.Mine();

  this->feFct_->GetElemSolution( elemSol, lpm.ptEl);
  BaseFE* ptFe = feSpace_->GetFe( lpm.ptEl->elemNum );
  if(this->bOps_.find(lpm.ptEl->regionId) != this->bOps_.end()){
    this->bOps_[lpm.ptEl->regionId]->CalcOpMat(bMat, lpm, ptFe);
    coefVec = bMat* (elemSol);
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
  out << "Complex: " << IsComplex() << "\n";
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
  isComplex_ =  std::is_same<TYPE,Complex>::value;

  // set inherited attributes
  dimType_ = CoefFunction::VECTOR;

}

template<class TYPE, bool TRANS> CoefFunctionFlux<TYPE,TRANS>::
~CoefFunctionFlux() {


}
template<class TYPE, bool TRANS> void CoefFunctionFlux<TYPE,TRANS>::
AddIntegrator( BaseBDBInt* form,  
               RegionIdType region ) {
#ifdef USE_OPENMP
  if(omp_get_num_threads()!=1){
    EXCEPTION("Call from parallel region which is not safe.")
  }
#endif
  std::map<RegionIdType, BaseBDBInt* >& forms = forms_.Mine(0);

  // check if we need to skip the integrator assignment if e.g. the requested integrator name differs
  if( integratorName_.empty() || integratorName_==form->GetName() ) {
    // check if region has already an integrator assigned
    if( forms.find(region) != forms.end() ) {
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
    //now we clone each integrator and we already checked if we are in single thread region
#ifdef USE_OPENMP
    for(UInt i=0;i<forms_.GetNumSlots();++i){
      forms_.Mine(i)[region] = form->Clone();
    }
#else
    forms_.Mine()[region] = form;
#endif
  } else {
    WARN("Skipped adding integrator " << form->GetName() << " to " << this->GetName());
  }
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

  std::map<RegionIdType, BaseBDBInt* >& forms = this->forms_.Mine();
  Vector<TYPE>& elemSol = elemSol_.Mine();

  this->feFct_->GetElemSolution( elemSol, lpm.ptEl);
  BaseBDBInt* bdb = forms[lpm.ptEl->regionId];
  if( !bdb) {
    coefVec.Resize(res_->dofNames.GetSize());
    coefVec.Init();
    return;
  }
  
  // switch depending on if 
  if( TRANS ) {
    bdb->ApplydATransMat(coefVec, elemSol, lpm );
  } else {
    // this applies physical design in dData_->GetTensor(dMat_,lpm) if we do optimization
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
  out << "Complex: " << IsComplex() << "\n";
  out << "Result: " << 
      SolutionTypeEnum.ToString(feFct_->GetResultInfo()->resultType );
  return out.str();
}

template class CoefFunctionFlux<Double,false>;
template class CoefFunctionFlux<Complex,false>;
template class CoefFunctionFlux<Double,true>;
template class CoefFunctionFlux<Complex,true>;

// ==========================================================================
//  COEFFICIENT FUNCTION FOR EIGENSTRESSES
// ==========================================================================
// Example for using CoefFunctionEigen: See MechPDE::DefinePostProcResults()
// Principal Stresses and Strains
CoefFunctionEigen::CoefFunctionEigen( shared_ptr<BaseFeFunction> feFct,
                  shared_ptr<ResultInfo> info,
                  PtrCoefFct stressCoef,
                  Double factor )
                  :CoefFunctionFormBased() {
  feFct_ = dynamic_pointer_cast<FeFunction<Double> >(feFct);
  res_ = info;
  factor_ = factor;
  isComplex_ = false; //this eigenvalue solver can only handle real values!

  // set inherited attributes
  dimType_ = CoefFunction::VECTOR;
  stressCoef_ = stressCoef;
}

CoefFunctionEigen::~CoefFunctionEigen() {

}

void CoefFunctionEigen::GetVector(Vector<Double>& coefVec, const LocPointMapped& lpm) {

  //Gets the stress/strain values from the CoefFunction class which is passed through stressCoef
  stressCoef_->GetVector(coefVec, lpm);

  GetEigenFromCoefVec(coefVec);
}

UInt CoefFunctionEigen::GetVecSize() const{
     return this->res_->dofNames.GetSize();
}

void CoefFunctionEigen::GetEigenFromCoefVec(Vector<Double> &solVec)
  {
  //This function calculates the principal stresses and principal strains from a given
  //stress or strain state in coefVec and return it to solVec.
  //Generic result, so it must be "sliced" in MechPDE with MathParser
  //
  //Returns principal stresses/strains in the following order:
  //-minimum principal stress vector
  //-medium principal stress vector (for 3D)
  //-maximum principal stress vector
  //-minimim principal stress scalar
  //-medium principal stress scalar (for 3D)
  //-maximum principal stress scalar

  //Step 1: Conversion of this (coefVec) to a std::vector a in the correct order

  std::vector<double> lp_a, lp_w, lp_work;

  // throwing otherwise "Not a valid Voigt stress/strain tensor for eigenvalue calculation"
  // is expensive as exception throwing has a slight overhead
  assert(solVec.GetSize() == 3 || solVec.GetSize() == 6);

  int lp_n = solVec.GetSize() == 3 ? 2 : 3;

  char lp_uplo = 'L';
  char lp_jobz = 'V'; //Documentation in LAPACK: dsyev
  int lp_lda = lp_n;
  int lp_info;
  int lp_lwork;

  assert(lp_n == 2 || lp_n == 3);

  if (lp_n == 2){                             //2D-case
    lp_lwork = 68;                            //68 is the optimum value for a 2x2 matrix
    lp_a.push_back(solVec.GetDoubleEntry(0)); //Reordering from Voigt notation to column-major matrix
    lp_a.push_back(solVec.GetDoubleEntry(2)); //style in FORTRAN
    lp_a.push_back(solVec.GetDoubleEntry(2));
    lp_a.push_back(solVec.GetDoubleEntry(1));
  }
  else {                                      //3D-case
    lp_lwork = 102;                           //102 is the optimum value for a 3x3 matrix
    lp_a.push_back(solVec.GetDoubleEntry(0)); //Reordering from Voigt notation to column-major matrix
    lp_a.push_back(solVec.GetDoubleEntry(5)); //style in FORTRAN
    lp_a.push_back(solVec.GetDoubleEntry(4));
    lp_a.push_back(solVec.GetDoubleEntry(5));
    lp_a.push_back(solVec.GetDoubleEntry(1));
    lp_a.push_back(solVec.GetDoubleEntry(3));
    lp_a.push_back(solVec.GetDoubleEntry(4));
    lp_a.push_back(solVec.GetDoubleEntry(5));
    lp_a.push_back(solVec.GetDoubleEntry(2));
  }

  lp_w.resize(lp_n);
  lp_work.resize(lp_lwork);

  //Step 2: Eigensolver

  solVec.Resize(lp_n * (lp_n+1));
  dsyev(&lp_jobz, &lp_uplo, &lp_n, &*lp_a.begin(), &lp_lda, &*lp_w.begin(), &*lp_work.begin(), &lp_lwork, &lp_info);

  //Step 3: Back-conversion of a to solVec
  for(int i = 0; i < lp_n; ++i) {
    for (int j = 0; j < lp_n; ++j){
      solVec[i*lp_n + j] = lp_w[i] * lp_a[i*lp_n + j]; //Stacking the eigenvectors in a 4x1 (2D) or 9x1(3D) vector
    }
  }
  for(int i = 0; i < lp_n; ++i) {
    solVec[lp_n * lp_n + i] = lp_w[i];
  }

  //lp_w: vector of eigenvalues
  //lp_a: vector of normalized eigenvectors
  //solVec before dsyev: "stacked" stress or strain tensor
  //solVec after dsyev: eigenvectors and eigenvalues in order as mentioned above.

}

std::string CoefFunctionEigen::ToString() const {
  std::stringstream out;
  out << "CoefFunctionEigen\n";
  out << "Complex: " << IsComplex() << "\n";
  out << "ApplyTransposed: false. No such implementation in CoefFunctionEigen\n";
  out << "Result: " <<
      SolutionTypeEnum.ToString(feFct_->GetResultInfo()->resultType );
  return out.str();
}

// ==========================================================================
//  COEFFICIENT FUNCTION BASED ON KERNEL OF BDB-INTEGRATOR
// ==========================================================================
template<class TYPE> CoefFunctionBdBKernel<TYPE>::
CoefFunctionBdBKernel(shared_ptr<BaseFeFunction> feFct,
                      TYPE factor ) 
: CoefFunctionFormBased() {
  feFct_ = dynamic_pointer_cast<FeFunction<TYPE> >(feFct);
  factor_ = factor;
  
  isComplex_ =  std::is_same<TYPE,Complex>::value;

  // set inherited attributes
  this->dimType_ = CoefFunction::SCALAR;
 }

template<class TYPE> CoefFunctionBdBKernel<TYPE>::
  ~CoefFunctionBdBKernel() {
}
template<class TYPE> void CoefFunctionBdBKernel<TYPE>::GetScalar( TYPE& coefScal, const LocPointMapped& lpm )
{
  // energy density is factor * elemSol^T * kernel * elemSol
  Vector<TYPE>& elemSol = elemSol_.Mine();
  std::map<RegionIdType, BaseBDBInt* >& forms = this->forms_.Mine();
  this->feFct_->GetElemSolution( elemSol, lpm.ptEl);
  
  Vector<TYPE> temp;
  if( !forms[lpm.ptEl->regionId]->IsComplex() ) {
    Matrix<Double>& kernelR = kernelR_.Mine();
    forms[lpm.ptEl->regionId]->CalcKernel(kernelR, lpm);
    temp = kernelR * elemSol;
  } else {
    Matrix<TYPE>& kernel = kernel_.Mine();
    forms[lpm.ptEl->regionId]->CalcKernel(kernel, lpm);
    temp = kernel * elemSol;
  }
  coefScal = temp.Inner(elemSol) * factor_;
}
template<class TYPE> std::string CoefFunctionBdBKernel<TYPE>::ToString() const
{
  std::stringstream out;
  out << "CoefFunctionBdBKernel\n";
  out << "Complex: " << IsComplex() << "\n";
  out << "Result: " << SolutionTypeEnum.ToString(feFct_->GetResultInfo()->resultType );
  return out.str();
}


template<class TYPE> CoefFunctionDyadicStrain<TYPE>::CoefFunctionDyadicStrain(shared_ptr<BaseFeFunction> feFct) : CoefFunctionFormBased()
{
  std::map<RegionIdType, BaseBDBInt* >& forms = this->forms_.Mine();
  LOG_DBG2(cff) << "CFDS:CFDS #forms=" << forms.size();

  feFct_ = dynamic_pointer_cast<FeFunction<TYPE> >(feFct);

  isComplex_ =  std::is_same<TYPE,Complex>::value;

  // set inherited attributes
  this->dimType_ = CoefFunction::TENSOR;
}

template<class TYPE> CoefFunctionDyadicStrain<TYPE>::~CoefFunctionDyadicStrain() { }

template<class TYPE> void CoefFunctionDyadicStrain<TYPE>::GetTensorSize(unsigned int& numRows, unsigned int& numCols ) const
{
  const std::map<RegionIdType, BaseBDBInt* > & forms = this->forms_.ConstMine();
  LOG_DBG2(cff) << "CFDS:GTS #forms=" << forms.size();
  numRows = domain->GetGrid()->GetDim() == 2 ? 3 : 6;
  numCols = numRows;
}

template<class TYPE> void CoefFunctionDyadicStrain<TYPE>::GetTensor(Matrix<TYPE>& tensor, const LocPointMapped& lpm)
{
  std::map<RegionIdType, BaseBDBInt* >& forms = this->forms_.Mine();

  LOG_DBG2(cff) << "CFDS:GT #forms=" << forms.size();
  BaseBDBInt* form = forms[lpm.ptEl->regionId];

  // LOG_DBG2(cff) << "CFDS:GT lpm=" << lpm.ptEl->elemNum << " form=" << form->GetName() << " B=" << form->GetBOp()->GetName();

  // Matrix<TYPE> B;
  // shared_ptr<FeSpace> mySpace = domain->GetStdPDE("mechanic")->GetFeFunction(MECH_DISPLACEMENT)->GetFeSpace();
  // BaseFE* ptFe = mySpace->GetFe( lpm.ptEl->elemNum);
  // form->GetBOp()->CalcOpMat(B, lpm, ptFe);
  // LOG_DBG2(cff) << "CFDS:GT B=" << B.ToString();


  // the dyadic product (B^T u) (u^T B) meant for integration, so not simply strain^2
  Vector<TYPE>& elemSol = elemSol_.Mine();
  this->feFct_->GetElemSolution(elemSol, lpm.ptEl);
  LOG_DBG2(cff) << "CFDS:GT u=" << elemSol.ToString();

  Vector<TYPE> tmp;
  form->ApplyBMat(tmp, elemSol, lpm);
  // LOG_DBG2(cff) << "CFDS:GT B*u=" << tmp.ToString();

  // Vector<TYPE> manual;
  // manual = B * ;
  // LOG_DBG2(cff) << "CFDS:GT manual=" << manual.ToString();

  tensor.DyadicMult(tmp, tmp.Conj());
  LOG_DBG2(cff) << "CFDS:GT -> " << tensor.ToString();
}

template<class TYPE> std::string CoefFunctionDyadicStrain<TYPE>::ToString() const
{
  std::stringstream out;
  out << "CoefFunctionDyadicStrain\n";
  out << "Complex: " << IsComplex() << "\n";
  out << "Result: " << SolutionTypeEnum.ToString(feFct_->GetResultInfo()->resultType );
  return out.str();
}

// --------------------------- ddd
template<class TYPE> CoefFunctionQuadSol<TYPE>::CoefFunctionQuadSol(shared_ptr<BaseFeFunction> feFct) : CoefFunctionFormBased()
{
  feFct_ = dynamic_pointer_cast<FeFunction<TYPE> >(feFct);

  isComplex_ =  std::is_same<TYPE,Complex>::value;

  // set inherited attributes
  this->dimType_ = CoefFunction::SCALAR;
}

template<class TYPE> CoefFunctionQuadSol<TYPE>::~CoefFunctionQuadSol() { }


template<class TYPE> void CoefFunctionQuadSol<TYPE>::GetScalar(TYPE& coefScal, const LocPointMapped& lpm)
{
  Vector<TYPE>& elemSol = elemSol_.Mine();
  // the dyadic product (B^T u) (u^T B) meant for integration, so not simply strain^2
  this->feFct_->GetElemSolution(elemSol, lpm.ptEl);
  LOG_DBG2(cff) << "CFSQ:GS u=" << elemSol.ToString();

  coefScal = elemSol.Inner();
  LOG_DBG2(cff) << "CFSQ:GS -> " << coefScal;
}

template<class TYPE> std::string CoefFunctionQuadSol<TYPE>::ToString() const
{
  std::stringstream out;
  out << "CoefFunctionQuadSol\n";
  out << "Complex: " << IsComplex() << "\n";
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
  std::stringstream out;
  out << "CoefFunctionLBM\n";
  out << "Complex: " << IsComplex() << "\n";
  return out.str();
}

//----------- Material TENSOR -----------

template<class TYPE, App::Type APP> CoefFunctionHomogenization<TYPE,APP>::CoefFunctionHomogenization(shared_ptr<BaseFeFunction> feFct, MaterialTensorNotation notation) : CoefFunctionFormBased()
{
  std::map<RegionIdType, BaseBDBInt* >& forms = this->forms_.Mine();

  LOG_DBG2(cff) << "CFH:CFH #forms=" << forms.size();

  feFct_ = dynamic_pointer_cast<FeFunction<TYPE> >(feFct);

  isComplex_ =  std::is_same<TYPE,Complex>::value;

  if (APP == App::MECH)
    assert( notation != NO_NOTATION);

  notation_ = notation;

  // set inherited attributes
  dimType_ = APP == App::ELEC ? CoefFunction::SCALAR : CoefFunction::TENSOR;

}

template<class TYPE, App::Type APP> CoefFunctionHomogenization<TYPE,APP>::~CoefFunctionHomogenization() { }

template<class TYPE, App::Type APP> unsigned int CoefFunctionHomogenization<TYPE,APP>::GetVecSize() const
{
  const std::map<RegionIdType, BaseBDBInt* >& forms = this->forms_.ConstMine();
  LOG_DBG2(cff) << "CFH:GVS #forms=" << forms.size();

  unsigned int numRows = 0;
  unsigned int numCols = 0;

  forms.begin()->second->GetCoef()->GetTensorSize(numRows, numCols);
  assert(numRows == domain->GetOptimization()->GetMultipleExcitation()->GetNumberHomogenization(APP));
  assert(numCols == numRows);

  if(numRows == 3)
    return 6;
  else
    return 21;
}

template<class TYPE, App::Type APP> void CoefFunctionHomogenization<TYPE,APP>::GetVector(Vector<TYPE>& vec, const LocPointMapped& lpm)
{
  std::map<RegionIdType, BaseBDBInt* >& forms = this->forms_.Mine();
  BaseBDBInt* form = forms[lpm.ptEl->regionId];

  LOG_DBG(cff) << "CFH:GTV #forms=" << forms.size() << " form=" << form->GetName() << " analytic=" << form->GetCoef()->IsAnalytic() << " notation=" << notation_;

  // the tensor here is always Voigt notation as we have Voigt on the simulation side of CFS (and Hill-Mandel on the optimization side)
  Matrix<TYPE> tensor;
  form->GetCoef()->GetTensor(tensor, lpm);

  if(notation_ == HILL_MANDEL) {
    MaterialTensor<TYPE> mt(VOIGT, &tensor, false);
    mt.ToHillMandel();
    tensor = mt.GetMatrix(HILL_MANDEL);
  }

  tensor.ConvertToVec_UpperTriangular(vec);
  LOG_DBG3(cff) << "CFH:GV tensor=" << tensor.ToString();
  LOG_DBG2(cff) << "CFH:GV -> " << vec.ToString();

}

template<class TYPE, App::Type APP> std::string CoefFunctionHomogenization<TYPE,APP>::ToString() const
{
  std::stringstream out;
  out << "CoefFunctionStiffness\n";
  out << "Complex: " << IsComplex() << "\n";
  out << "res " << SolutionTypeEnum.ToString(feFct_->GetResultInfo()->resultType );
  return out.str();
}

template<class TYPE, App::Type APP> void CoefFunctionHomogenization<TYPE,APP>::GetTensorSize(unsigned int& numRows, unsigned int& numCols ) const
{
  const std::map<RegionIdType, BaseBDBInt* >& forms = this->forms_.ConstMine();
  LOG_DBG2(cff) << "CFH:GTS #forms=" << forms.size();

  forms.begin()->second->GetCoef()->GetTensorSize(numRows, numCols);

  assert(numRows == domain->GetOptimization()->GetMultipleExcitation()->GetNumberHomogenization(APP));
  assert(numCols == numRows);
}

template<class TYPE, App::Type APP> void CoefFunctionHomogenization<TYPE,APP>::GetTensor(Matrix<TYPE>& tensor, const LocPointMapped& lpm)
{
  std::map<RegionIdType, BaseBDBInt* >& forms = this->forms_.Mine();
  BaseBDBInt* form = forms[lpm.ptEl->regionId];

  LOG_DBG(cff) << "CFH:GT #forms=" << forms.size() << " form=" << form->GetName() << " analytic=" << form->GetCoef()->IsAnalytic() << " notation=" << notation_;

  // the tensor here is always Voigt notation as we have Voigt on the simulation side of CFS (and Hill-Mandel on the optimization side)
  form->GetCoef()->GetTensor(tensor, lpm);

  if(notation_ == HILL_MANDEL && APP == App::MECH) {
    MaterialTensor<TYPE> mt(VOIGT, &tensor, false);
    mt.ToHillMandel();
    tensor = mt.GetMatrix(HILL_MANDEL);
  }

  LOG_DBG2(cff) << "CFH:GT -> " << tensor.ToString();
}

template<class TYPE, App::Type APP> void CoefFunctionHomogenization<TYPE,APP>::GetScalar(TYPE& coefScal, const LocPointMapped& lpm)
{
  std::map<RegionIdType, BaseBDBInt* >& forms = this->forms_.Mine();
  BaseBDBInt* form = forms[lpm.ptEl->regionId];

  form->GetCoef()->GetScalar(coefScal, lpm);
}

template class CoefFunctionBdBKernel<double>;
template class CoefFunctionBdBKernel<Complex>;

template class CoefFunctionDyadicStrain<double>;
template class CoefFunctionDyadicStrain<Complex>;

template class CoefFunctionQuadSol<double>;
template class CoefFunctionQuadSol<Complex>;

template class CoefFunctionHomogenization<double, App::MECH>;
template class CoefFunctionHomogenization<Complex, App::MECH>;
template class CoefFunctionHomogenization<double, App::HEAT>;

template class CoefFunctionHomogenization<Complex, App::ELEC>;
template class CoefFunctionHomogenization<double, App::ELEC>;
template class CoefFunctionHomogenization<double, App::MAG>;
template class CoefFunctionLBM<double>;

   
} // end of namespace
