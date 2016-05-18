namespace CoupledField {

template< class VEC_DATA_TYPE >
KXIntegrator<VEC_DATA_TYPE>::KXIntegrator( BiLinearForm * biLinForm,
                                           VEC_DATA_TYPE factor,
                                           shared_ptr<BaseFeFunction> feFct) {
  this->name_ = "RhsBUIntegrator";
  this->factor_ = factor;
  assert(biLinForm != NULL);
  this->form_ = biLinForm;
  this->feFct_ = dynamic_pointer_cast<FeFunction<VEC_DATA_TYPE> > (feFct );

}

//! Copy constructor
template<class VEC_DATA_TYPE >
KXIntegrator<VEC_DATA_TYPE>::KXIntegrator(const KXIntegrator& right ){
  this->form_ = right.form_->Clone();
  this->factor_ = right.factor_;
  this->feFct_ = right.feFct_;
}

template<class VEC_DATA_TYPE >
void KXIntegrator<VEC_DATA_TYPE>::CalcElemVector( Vector<VEC_DATA_TYPE> & elemVec,
                                                  EntityIterator& ent){

  Matrix<VEC_DATA_TYPE> elemMat;
  Vector<VEC_DATA_TYPE> elemSol;

  // calculate element matrix
  this->form_->CalcElementMatrix(elemMat, ent, ent );

  // obtain element solution
  feFct_->GetEntitySolution( elemSol, ent );
  elemVec = (elemMat * elemSol) * factor_;
}

//template<class VEC_DATA_TYPE >
//void KXIntegrator<VEC_DATA_TYPE>::CalcRhsVector( Vector<VEC_DATA_TYPE> & elemVec,
//                                              shared_ptr<CoefFunction > rhsCoefs,
//                                                  EntityIterator& ent){
//
//  Matrix<VEC_DATA_TYPE> elemMat;
//  Vector<VEC_DATA_TYPE> elemSol;
//  Vector<VEC_DATA_TYPE> cVec;
//  LocPointMapped lp;
//
//  // calculate element matrix
//  this->form_->CalcElementMatrix(elemMat, ent, ent );
//
//  // obtain element solution
//  //feFct_->GetEntitySolution( elemSol, ent );
//
//    if( rhsCoefs->GetDimType() == CoefFunction::SCALAR ) {
//      elemSol.Resize(1);
//      rhsCoefs->GetScalar(elemSol[0],lp);
//    } else {
//      rhsCoefs->GetVector(elemSol,lp);
//    }
//
//  elemVec = (elemMat * elemSol) * factor_;
//}

} // end of namespace
