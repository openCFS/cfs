#include "BiLinWrappredLinForm.hh"
#include "Forms/LinForms/LinearForm.hh"

namespace CoupledField {


BiLinWrappedLinForm::BiLinWrappedLinForm(LinearForm* linForm,
                                         bool assembleTranposed)
: BiLinearForm(linForm->IsCoordUpdate()) {
  linForm_ = linForm;
  assembleTransposed_ = assembleTranposed;
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
    linForm_->CalcElemVector(elemVec, ent1);
    elemMat.Resize(1, elemVec.GetSize());
    for( UInt i = 0; i < elemVec.GetSize(); i++) {
      elemMat[0][i] = elemVec[i];  
    }
  } else {
    linForm_->CalcElemVector(elemVec, ent2);
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
    linForm_->CalcElemVector(elemVec, ent1);
    elemMat.Resize(1, elemVec.GetSize());
    for( UInt i = 0; i < elemVec.GetSize(); i++) {
      elemMat[0][i] = elemVec[i];  
    }
  } else {
    linForm_->CalcElemVector(elemVec, ent2);
    elemMat.Resize(elemVec.GetSize(), 1);
    for( UInt i = 0; i < elemVec.GetSize(); i++) {
      elemMat[i][0] = elemVec[i];  
    }
  }
}


} // Namespace

