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

template< class B_OP, class VEC_DATA_TYPE, bool SURFACE >
BDUIntegrator<B_OP,VEC_DATA_TYPE,SURFACE>::
BDUIntegrator(VEC_DATA_TYPE factor,
             shared_ptr<CoefFunction > rhsCoef, 
             shared_ptr<CoefFunction > dCoef,
             bool coordUpdate )
             :LinearForm( coordUpdate ){
  factor_ = factor;
  this->name_ = "RhsBDUIntegrator";

//  assert(rhsCoef->GetDimType() == CoefFunction::VECTOR);
//#ifndef NDEBUG
//  if(rhsCoef->GetDimType() != CoefFunction::VECTOR){
//    Exception("BDB integrator expects the coefficient function to be vectorial!\n \
//                     For scalar valued Things, create a vectorial function with one component");
//  }
//#endif
  this->rhsCoefs_ = rhsCoef;
  this->dCoef_ = dCoef;

}

template< class B_OP, class VEC_DATA_TYPE, bool SURFACE >
BDUIntegrator<B_OP,VEC_DATA_TYPE,SURFACE>::
BDUIntegrator(VEC_DATA_TYPE factor,
             shared_ptr<CoefFunction > rhsCoef,
             shared_ptr<CoefFunction > dCoef,
             const std::set<RegionIdType>& volRegions,
             bool coordUpdate )
             :LinearForm( coordUpdate ){
  factor_ = factor;
  this->name_ = "RhsBDUIntegrator";

  assert(rhsCoef->GetDimType() == CoefFunction::VECTOR);
#ifndef NDEBUG
  if(rhsCoef->GetDimType() != CoefFunction::VECTOR){
    Exception("BDB integrator expects the coefficient function to be vectorial!\n \
                     For scalar valued Things, create a vectorial function with one component");
  }
#endif
  this->rhsCoefs_ = rhsCoef;
  this->dCoef_ = dCoef;
  volRegions_ = volRegions;
}

  template<class B_OP, class VEC_DATA_TYPE, bool SURFACE >
  void BDUIntegrator<B_OP,VEC_DATA_TYPE,SURFACE>::
  CalcElemVector( Vector<VEC_DATA_TYPE> & elemVec,
                  EntityIterator& ent){
    // Declare necessary variables
     const Elem* ptElem = ent.GetElem();
     Matrix<VEC_DATA_TYPE> dMat, bdMat, bMat;
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

     elemVec.Resize( nrFncs * Bdim_);
     elemVec.Init();

     // Loop over all integration points
     LocPointMapped lp;
     for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

       // Calculate for each integration point the LocPointMapped
       if (SURFACE) {
         lp.Set( intPoints[i], esm, volRegions_, weights[i] );
       } else {
         lp.Set( intPoints[i], esm, weights[i] );
       }

       // obtain d matrix
       this->dCoef_->GetTensor( dMat, lp );
       
       //calc factor
       fac = VEC_DATA_TYPE(lp.jacDet * weights[i]);
       fac *= factor_;
       
       // Call the CalcBMat()-method
       operator_.CalcOpMatTransposed( bMat, lp, ptFe);
       
       bdMat.Resize(nrFncs * B_OP::DIM_DOF, dMat.GetNumRows());
       // Calculate BdMat
       bMat.Mult_Blas(dMat, bdMat,false,false,fac,0.0);

       rhsCoefs_->GetVector(cVec,lp);  
       //elemVec += bMat * cVec * fac;
       elemVec += bdMat * cVec;
     }

  }
}
