// =====================================================================================
// 
// 
//    Description:  Implementation file for BBIntegrators
//                  TAKE care:
//                  This file is just for inclusion in the header file!
// 
//        Version:  1.0
//        Created:  11/03/2011 08:04:43 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================
#include "ABInt.hh"

namespace CoupledField{


   template< class COEF_DATA_TYPE, class B_DATA_TYPE>
   ABInt<COEF_DATA_TYPE, B_DATA_TYPE>::
   ABInt(BaseBOperator * aOp, BaseBOperator * bOp, 
         PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
         bool coordUpdate )
     : BBInt<COEF_DATA_TYPE, B_DATA_TYPE>(bOp, scalCoef, factor, coordUpdate ) {
     this->name_ = "ABInt";
     this->aOperator_ = aOp;
     
     // Note: In general the AB-Integrator is not symmetric, as is should 
     // get only used, if the A- and B differential operators are distinct.
     this->isSymmetric_ = false;
   }

   template< class COEF_DATA_TYPE, class B_DATA_TYPE>
   void ABInt<COEF_DATA_TYPE, B_DATA_TYPE>::
   CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                      EntityIterator& ent1,
                      EntityIterator& ent2 ) {
     // Extract physical element
     const Elem* ptElem = ent1.GetElem();

     MAT_DATA_TYPE fac = 0.0;

     // Obtain FE element from feSpace and integration scheme
     IntegOrder order1, order2;
     IntScheme::IntegMethod method1, method2;
     BaseFE* ptFeA = this->ptFeSpace1_->GetFe( ent1, method1, order1 );
     BaseFE* ptFeB = this->ptFeSpace2_->GetFe( ent1, method2, order2 );

     const UInt nrFncsA = ptFeA->GetNumFncs();
     const UInt nrFncsB = ptFeB->GetNumFncs();

     // Get shape map from grid
     shared_ptr<ElemShapeMap> esm = 
         domain->GetGrid()->GetElemShapeMap( ptElem, this->coordUpdate_ );

     // Get integration points
     StdVector<LocPoint> intPoints;
     StdVector<Double> weights;
     this->intScheme_->GetIntPoints( Elem::GetShapeType(ptElem->type), 
                                     method1, order1, method2, order2,
                                     intPoints, weights );

     elemMat.Resize( nrFncsA * aOperator_->GetDimDof(), 
                     nrFncsB * this->bOperator_->GetDimDof() );
     elemMat.Init();

#define USE_BLAS_VERSION
     // Loop over all integration points
     LocPointMapped lp;
     const UInt numIntPts = intPoints.GetSize();
     for( UInt i = 0; i < numIntPts; ++i  ) {

       // Calculate for each integration point the LocPointMapped
       lp.Set( intPoints[i], esm );

       // Calculate A-matrix (first differential operator)
       this->aOperator_->CalcOpMat( aMat_, lp, ptFeA );
       
       // Calculate B-matrix (second differential operator)
       this->bOperator_->CalcOpMat( this->bMat_, lp, ptFeB );

       // Calculate scalar factor
       this->coefScalar_->GetScalar(fac, lp);
       fac *= MAT_DATA_TYPE(lp.jacDet * weights[i]); 


#ifdef USE_BLAS_VERSION
       aMat_.Mult_Blas(this->bMat_,elemMat,true,false,this->factor_*fac,1.0);
#else
       elemMat += Transpose(aMat_) * this->bMat * this->factor_;
#endif

     }

   }

   template< class COEF_DATA_TYPE, class B_DATA_TYPE>
   SurfaceABInt<COEF_DATA_TYPE, B_DATA_TYPE>::
   SurfaceABInt( BaseBOperator * aOp, BaseBOperator * bOp,
                 PtrCoefFct scalCoef, MAT_DATA_TYPE factor,
                 const std::set<RegionIdType>& volRegions,
                 bool coordUpdate )
            : ABInt<COEF_DATA_TYPE,B_DATA_TYPE>(aOp, bOp, scalCoef, factor, coordUpdate ){
     this->name_ = "SurfaceABInt";
     volRegions_ = volRegions;
     this->isSymmetric_ = false;
   }

   template< class COEF_DATA_TYPE, class B_DATA_TYPE>
   SurfaceABInt<COEF_DATA_TYPE, B_DATA_TYPE>::
   SurfaceABInt( BaseBOperator * aOp, BaseBOperator * bOp,
                 const std::map< RegionIdType, PtrCoefFct >& regionCoefs,
                 MAT_DATA_TYPE factor,
                 const std::set<RegionIdType>& volRegions,
                 bool coordUpdate)
     : ABInt<COEF_DATA_TYPE,B_DATA_TYPE>(aOp, bOp,regionCoefs.begin()->second, 
                                         factor, coordUpdate )
   {
     this->name_ = "SurfaceABInt";
     volRegions_ = volRegions;
     this->isSymmetric_ = false;
     regionCoefs_ = regionCoefs;
   }

   template< class COEF_DATA_TYPE, class B_DATA_TYPE>
   void SurfaceABInt<COEF_DATA_TYPE, B_DATA_TYPE>::
   CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                      EntityIterator& ent1,
                      EntityIterator& ent2 ) {
     // Extract physical element
     const Elem* ptElem1 = ent1.GetElem();
     const Elem* ptElem2 = ent2.GetElem();

     Matrix<MAT_DATA_TYPE> aMat, bMat;
     MAT_DATA_TYPE fac = 0.0;

     // Obtain FE element from feSpace and integration scheme
     IntegOrder order1, order2;
     IntScheme::IntegMethod method1, method2;
     BaseFE* ptFeA = this->ptFeSpace1_->GetFe( ent1, method1, order1 );
     BaseFE* ptFeB = this->ptFeSpace2_->GetFe( ent2, method2, order2 );


     const UInt nrFncsA = ptFeA->GetNumFncs();
     const UInt nrFncsB = ptFeB->GetNumFncs();

     // Get shape map from grid
     shared_ptr<ElemShapeMap> esm1 =
         domain->GetGrid()->GetElemShapeMap( ptElem1, this->coordUpdate_ );
     shared_ptr<ElemShapeMap> esm2 =
         domain->GetGrid()->GetElemShapeMap( ptElem2, this->coordUpdate_ );

     // Get integration points
     StdVector<LocPoint> intPoints;
     StdVector<Double> weights;
     this->intScheme_->GetIntPoints( Elem::GetShapeType(ptElem1->type), 
                                     method1, order1, method2, order2,
                                     intPoints, weights );

     elemMat.Resize( nrFncsA * this->aOperator_->GetDimDof(),
                     nrFncsB * this->bOperator_->GetDimDof() );
     elemMat.Init();

#define USE_BLAS_VERSION
     // Loop over all integration points
     LocPointMapped lp1,lp2;
     const UInt numIntPts = intPoints.GetSize();
     for( UInt i = 0; i < numIntPts; ++i ) {

       // Calculate for each integration point the LocPointMapped
       lp1.Set( intPoints[i], esm1, volRegions_ );
       lp2.Set( intPoints[i], esm2, volRegions_ );

       // Calculate A-matrix (first differential operator)
       this->aOperator_->CalcOpMat( this->aMat_, lp1, ptFeA );

       // Calculate B-matrix (second differential operator)
       this->bOperator_->CalcOpMat( this->bMat_, lp2, ptFeB );

       // Calculate scalar factor
       if(!regionCoefs_.empty()) {
         regionCoefs_[lp1.lpmVol->ptEl->regionId]->GetScalar(fac, lp1);
       } else {
         this->coefScalar_->GetScalar(fac, lp1);
       }
       
       fac *= MAT_DATA_TYPE(lp1.jacDet * weights[i]);

#ifdef USE_BLAS_VERSION
       this->aMat_.Mult_Blas(this->bMat_, elemMat, true, false,
                             this->factor_*fac, 1.0);
#else
       elemMat += Transpose(this->aMat_) * thisbMat_ * this->factor_;
#endif

     }

   }
 // Explicit template instantiation
 template class ABInt<Double,Double>;
 template class ABInt<Double,Complex>;
 template class ABInt<Complex,Double>;
 template class ABInt<Complex,Complex>;
 template class SurfaceABInt<Double,Double>;
 template class SurfaceABInt<Double,Complex>;
 template class SurfaceABInt<Complex,Double>;
 template class SurfaceABInt<Complex,Complex>;
  
} // name
