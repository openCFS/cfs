// =====================================================================================
// 
//       Filename:  buInt.cc
// 
//    Description:  Implementation file for BUIntegrators
//                  TAKE care:
//                  This file is just for inclusion in the header file!
// 
//        Version:  1.0
//        Created:  11/03/2011 08:08:21 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================



namespace CoupledField{

template< class VEC_DATA_TYPE, bool SURFACE >
BUIntegrator<VEC_DATA_TYPE,SURFACE>::
BUIntegrator(BaseBOperator * bOp,
             VEC_DATA_TYPE factor,
             shared_ptr<CoefFunction > rhsCoef, 
             bool coordUpdate,
             bool fullEvaluation )
             : LinearForm( coordUpdate ),
               fullEvaluation_(fullEvaluation) {
  factor_ = factor;
  this->name_ = "RhsBUIntegrator";
  this->bOperator_= bOp;
//  assert(rhsCoef->GetDimType() == CoefFunction::VECTOR);
//#ifndef NDEBUG
//  if(rhsCoef->GetDimType() != CoefFunction::VECTOR){
//    Exception("BDB integrator expects the coefficient function to be vectorial!\n \
//                     For scalar valued Things, create a vectorial function with one component");
//  }
//#endif
  this->rhsCoefs_ = rhsCoef;

}

template< class VEC_DATA_TYPE, bool SURFACE >
BUIntegrator<VEC_DATA_TYPE,SURFACE>::
BUIntegrator(BaseBOperator * bOp,
             VEC_DATA_TYPE factor,
             shared_ptr<CoefFunction > rhsCoef,
             const std::set<RegionIdType>& volRegions,
             bool coordUpdate,
             bool fullEvaluation )
             : LinearForm( coordUpdate ), 
               fullEvaluation_(fullEvaluation) 
               {
  factor_ = factor;
  this->name_ = "RhsBUIntegrator";
  this->bOperator_= bOp;

  assert(rhsCoef->GetDimType() == CoefFunction::VECTOR ||
         rhsCoef->GetDimType() == CoefFunction::SCALAR);
#ifndef NDEBUG
  if(rhsCoef->GetDimType()  != CoefFunction::VECTOR &&
      rhsCoef->GetDimType() != CoefFunction::SCALAR){
    Exception("BDB integrator expects the coefficient function to be vectorial or scalar!");
  }
#endif
  this->rhsCoefs_ = rhsCoef;
  volRegions_ = volRegions;
}

  template< class VEC_DATA_TYPE, bool SURFACE >
  void BUIntegrator<VEC_DATA_TYPE,SURFACE>::
  CalcElemVector( Vector<VEC_DATA_TYPE> & elemVec,
                  EntityIterator& ent){

	 assert(rhsCoefs_->GetDimType() != CoefFunction::NO_DIM);

	  // Declare necessary variables
     const Elem* ptElem = ent.GetElem();
     Matrix<Double> bMat;
     Vector<VEC_DATA_TYPE> cVec;
     StdVector<LocPoint> intPoints;
     StdVector<Double> weights;
     UInt nrFncs = 0;
     VEC_DATA_TYPE fac;

     // Obtain FE element from feSpace and integration scheme
     IntegOrder order;
     IntScheme::IntegMethod method;
     BaseFE* ptFe = ptFeSpace_->GetFe( ent, method, order );

     nrFncs = ptFe->GetNumFncs();

     // Get shape map from grid
     shared_ptr<ElemShapeMap> esm =
         ent.GetGrid()->GetElemShapeMap( ptElem, this->coordUpdate_ );

     // Get integration points
     intScheme_->GetIntPoints( Elem::GetShapeType(ptElem->type), method, order, 
                               intPoints, weights );

     LocPointMapped lp;
     elemVec.Resize( nrFncs * Bdim_);
     elemVec.Init();
     
     // Pre-evaluate coefficient function in case of reduced accuracy
     if(! fullEvaluation_ ) {
       const ElemShape sh = Elem::shapes[ptElem->type];
       lp.Set( sh.midPointCoord, esm, sh.volume );
       if( rhsCoefs_->GetDimType() == CoefFunction::SCALAR ) {
         cVec.Resize(1);
         rhsCoefs_->GetScalar(cVec[0],lp);
       } else {
         rhsCoefs_->GetVector(cVec,lp);
       }
     }

     // Loop over all integration points
     for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

       // Calculate for each integration point the LocPointMapped
       if (SURFACE) {
         lp.Set( intPoints[i], esm, volRegions_, weights[i] );
       } else {
         lp.Set( intPoints[i], esm, weights[i] );
       }

       //calc factor
       fac = VEC_DATA_TYPE(lp.jacDet * weights[i]);
       fac *= factor_;
       // Call the CalcBMat()-method
       bOperator_->CalcOpMatTransposed( bMat, lp, ptFe);

       // Evaluate coefficient function in integration point
       // ( in case of full order)
       if( fullEvaluation_ ) {
         if( rhsCoefs_->GetDimType() == CoefFunction::SCALAR ) {
           cVec.Resize(1);
           rhsCoefs_->GetScalar(cVec[0],lp);
         } else {
           rhsCoefs_->GetVector(cVec,lp);
         }
       }
       elemVec += bMat * cVec * fac;
     }

  }
}
