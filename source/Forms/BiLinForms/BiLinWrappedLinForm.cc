#include "BiLinWrappedLinForm.hh"
#include "Forms/LinForms/LinearForm.hh"

namespace CoupledField {


BiLinWrappedLinForm::BiLinWrappedLinForm(LinearForm* linForm,
                                         bool assembleTranposed)
: BiLinearForm(linForm->IsCoordUpdate()) {
  linForm_ = linForm;
  assembleTransposed_ = assembleTranposed;
  name_ = linForm->GetName();
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
  /*linForm_->CalcElemVector(elemVec, ent1);
  if( assembleTransposed_ ){
    elemMat.Resize(1, elemVec.GetSize());
    for( UInt i = 0; i < elemVec.GetSize(); i++) {
      elemMat[0][i] = elemVec[i];
    }
  } else {
    elemMat.Resize(elemVec.GetSize(), 1);
    for( UInt i = 0; i < elemVec.GetSize(); i++) {
      elemMat[i][0] = elemVec[i];
    }
  }*/
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
  /*linForm_->CalcElemVector(elemVec, ent1);
  if( assembleTransposed_ ){
    elemMat.Resize(1, elemVec.GetSize());
    for( UInt i = 0; i < elemVec.GetSize(); i++) {
      elemMat[0][i] = elemVec[i];
    }
  } else {
    elemMat.Resize(elemVec.GetSize(), 1);
    for( UInt i = 0; i < elemVec.GetSize(); i++) {
      elemMat[i][0] = elemVec[i];
    }
  }*/
}

bool BiLinWrappedLinForm::IsSolDependent(){

  return linForm_->IsSolDependent();

}

bool BiLinWrappedLinForm::IsComplex(){

  return linForm_->IsComplex();

}

void BiLinWrappedLinForm::SetFeSpace( shared_ptr<FeSpace> feSpace ){

  //this->ptFeSpace1_ = feSpace;
  //this->ptFeSpace2_ = feSpace;
  linForm_->SetFeSpace(feSpace);
  //this->intScheme_ = feSpace->GetIntScheme();

}

void BiLinWrappedLinForm::SetFeSpace( shared_ptr<FeSpace> feSpace1, shared_ptr<FeSpace> feSpace2 ){

  //this->ptFeSpace1_ = feSpace1;
  //this->ptFeSpace2_ = feSpace2;
  if( assembleTransposed_ ){
    //this->intScheme_ = feSpace2->GetIntScheme();
    linForm_->SetFeSpace(feSpace2);
  } else {
    //this->intScheme_ = feSpace1->GetIntScheme();
    linForm_->SetFeSpace(feSpace1);
  }

}

} // Namespace

