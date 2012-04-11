#include "ResultFunctor.hh"

namespace CoupledField {

template<class TYPE> ResultFunctorIntegrate<TYPE>::
ResultFunctorIntegrate( PtrCoefFct coef,
                        shared_ptr<BaseFeFunction> feFct,
                        shared_ptr<ResultInfo> inf ) :
                        ResultFunctor( inf) {
  derivType_ = INTEGRATED;
  coef_ = coef;
  feFct_ = feFct;
}
    
template<class TYPE> ResultFunctorIntegrate<TYPE>::
~ResultFunctorIntegrate() {
  
}
  

template<class TYPE> void ResultFunctorIntegrate<TYPE>:: 
 EvalResult(shared_ptr<BaseResult> res ) {
  Result<TYPE>& actSol = static_cast<Result<TYPE>& >(*res);
  EntityIterator regionIt = actSol.GetEntityList()->GetIterator();
  shared_ptr<FeSpace> feSpace = feFct_->GetFeSpace();
  shared_ptr<IntScheme> intScheme = feSpace->GetIntScheme();
  Vector<TYPE>& vec = actSol.GetVector();
  vec.Resize( regionIt.GetSize() );

  // Loop over regions
  for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {
    ElemList actSDList(ptGrid_);
    actSDList.SetRegion( regionIt.GetRegion() );
    EntityIterator elemIt = actSDList.GetIterator();

    TYPE tempVal = 0.0;
    // loop over elements
    for ( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
      const Elem * el = elemIt.GetElem();


      // Obtain FE element from feSpace and integration scheme
      IntegOrder order;
      IntScheme::IntegMethod method;

      feSpace->GetFe( elemIt, method, order );
      // Get shape map from grid
      shared_ptr<ElemShapeMap> esm = 
          domain->GetGrid()->GetElemShapeMap( el, true );

      // Get integration points
      StdVector<LocPoint> intPoints;
      StdVector<Double> weights;
      intScheme->GetIntPoints( Elem::GetShapeType(el->type), method, order, 
                               intPoints, weights );

      // Loop over all integration points
      tempVal = 0.0;
      LocPointMapped lpm;
      TYPE elemVal = 0.0;
      for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

        // Calculate for each integration point the LocPointMapped
        lpm.Set( intPoints[i], esm );
        coef_->GetScalar(tempVal, lpm );
        elemVal += tempVal * lpm.jacDet * weights[i];
      } // loop integration points
      vec[regionIt.GetPos()] += elemVal;
    } // loop elements
  } // loop regions

}

template class ResultFunctorIntegrate<Complex>;
template class ResultFunctorIntegrate<Double>;
} // end of namespace
