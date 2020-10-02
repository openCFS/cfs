#include "BiLinWrappedLinForm.hh"
#include "Forms/LinForms/LinearForm.hh"

namespace CoupledField {


BiLinWrappedLinForm::BiLinWrappedLinForm(LinearForm* linForm,
                                         bool assembleTranposed)
: BiLinearForm(linForm->IsCoordUpdate()) {
  linForm_ = linForm;
  assembleTransposed_ = assembleTranposed;
  type_ = BILIN_WRAPPED_LIN;
  name_ = linForm->GetName();
}

BiLinWrappedLinForm::BiLinWrappedLinForm(const BiLinWrappedLinForm& right)
: BiLinearForm(right){
    this->linForm_ = right.linForm_->Clone();
    this->assembleTransposed_ = right.assembleTransposed_;
}

BiLinWrappedLinForm::~BiLinWrappedLinForm() {

  // we assume responsibility of the encapsulated LinearForm 
  // and delete it.
  delete linForm_;
}

void BiLinWrappedLinForm::CalcElementMatrix( Matrix<Double>& elemMat,
                                             EntityIterator& ent1,
                                             EntityIterator& ent2) {
  Vector<Double> elemVec;
  if( assembleTransposed_ ) {
    linForm_->CalcElemVector(elemVec, ent2);
    elemMat.Resize(1, elemVec.GetSize());
    for( UInt i = 0; i < elemVec.GetSize(); i++) {
      elemMat[0][i] = elemVec[i];
    }
  } else {
    linForm_->CalcElemVector(elemVec, ent1);
    elemMat.Resize(elemVec.GetSize(), 1);
    for( UInt i = 0; i < elemVec.GetSize(); i++) {
      elemMat[i][0] = elemVec[i];
    }
  }

}

void BiLinWrappedLinForm::CalcElementMatrix( Matrix<Complex>& elemMat,
                                             EntityIterator& ent1,
                                             EntityIterator& ent2) {
  Vector<Complex> elemVec;
  if( assembleTransposed_ ) {
    linForm_->CalcElemVector(elemVec, ent2);
    elemMat.Resize(1, elemVec.GetSize());
    for( UInt i = 0; i < elemVec.GetSize(); i++) {
      elemMat[0][i] = elemVec[i];
    }
  } else {
    linForm_->CalcElemVector(elemVec, ent1);
    elemMat.Resize(elemVec.GetSize(), 1);
    for( UInt i = 0; i < elemVec.GetSize(); i++) {
      elemMat[i][0] = elemVec[i];
    }
  }

}

bool BiLinWrappedLinForm::IsSolDependent(){

  return linForm_->IsSolDependent();

}

bool BiLinWrappedLinForm::IsComplex() const {

  return linForm_->IsComplex();

}

void BiLinWrappedLinForm::SetFeSpace( shared_ptr<FeSpace> feSpace ){

  linForm_->SetFeSpace(feSpace);

}

void BiLinWrappedLinForm::SetFeSpace( shared_ptr<FeSpace> feSpace1, shared_ptr<FeSpace> feSpace2 ){

  if( assembleTransposed_ ){
    linForm_->SetFeSpace(feSpace2);
  } else {
    linForm_->SetFeSpace(feSpace1);
  }

}

} // Namespace

