// =====================================================================================
// 
//       Filename:  bbInt.cc
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

namespace CoupledField{


   template< class B_OP, class MAT_DATA_TYPE> 
   BBInt<B_OP,MAT_DATA_TYPE>::
   BBInt(shared_ptr<CoefFunction> scalCoef, MAT_DATA_TYPE factor) {
     this->name_ = "BBInt";
     isSymmetric_ = true;
     this->scalCoef_ = scalCoef;
     factor_ = factor;
   }

   template<class B_OP, class MAT_DATA_TYPE>
   void BBInt<B_OP,MAT_DATA_TYPE>::CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                                                      EntityIterator& ent1,
                                                      EntityIterator& ent2) {
     // Extract physical element
     const Elem* ptElem = ent1.GetElem();

     Matrix<MAT_DATA_TYPE> bMat;
     MAT_DATA_TYPE fac = 0.0;

     // Obtain FE element from feSpace and integration scheme
     shared_ptr<IntScheme> intScheme;
     BaseFE* ptFe = ptFeSpace1_->GetFe( ent1, intScheme );
     UInt nrFncs = ptFe->GetNumFncs();

     // Get shape map from grid
     shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap( ptElem );

     // Get integration points
     StdVector<LocPoint> intPoints;
     StdVector<Double> weights;
     intScheme->GetIntPoints( Elem::GetShapeType(ptElem->type), intPoints, weights );

     elemMat.Resize( nrFncs * B_OP::DIM_DOF );
     elemMat.Init();

#define USE_BLAS_VERSION
     // Loop over all integration points
     LocPointMapped lp;
     for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

       // Calculate for each integration point the LocPointMapped
       lp.Set( intPoints[i], esm );

       // Call the CalcBMat()-method
       bOperator_.CalcOpMat( bMat, lp, ptFe);

       // Calculate scalar factor
       scalCoef_->GetScalar(fac, lp);
       fac *= MAT_DATA_TYPE(lp.jacDet * weights[i]); 

       bOperator_.TransformJacDet(fac,lp,ptFe);

#ifdef USE_BLAS_VERSION
       bMat.Mult_Blas(bMat, elemMat, true, false, this->factor_ * fac, 1.0);
#else
       elemMat += Transpose(bMat) * bMat * this->factor_ * fac;
#endif

     }

   }

   template< class B_OP, class MAT_DATA_TYPE>
   void BBIntMassEdge<B_OP,MAT_DATA_TYPE>::
   CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                      EntityIterator& ent1,
                      EntityIterator& ent2 ) {
     // Extract physical element
     const Elem* ptElem = ent1.GetElem();

     Matrix<MAT_DATA_TYPE> bMat,bbMat;
     MAT_DATA_TYPE fac = 0.0;

     // Obtain FE element from feSpace and integration scheme
     shared_ptr<IntScheme> intScheme;
     //TODO: we have to find another solution for this. The only reason for this class is within
     //these two lines. This functionaliy has to be moved to a pre-Integration function
     FeHCurlHi* ptFe = static_cast<FeHCurlHi*>(this->ptFeSpace1_->GetFe( ent1, intScheme ));
     // Special: Only use lower order functions
     ptFe->SetOnlyLowestOrder(true);

     UInt nrFncs = ptFe->BaseFE::GetNumFncs();

     // Get shape map from grid
     shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap( ptElem );

     // Get integration points
     StdVector<LocPoint> intPoints;
     StdVector<Double> weights;
     intScheme->GetIntPoints( Elem::GetShapeType(ptElem->type), intPoints, weights );

     elemMat.Resize( nrFncs * B_OP::DIM_DOF );
     elemMat.Init();

     // Loop over all integration points
     LocPointMapped lp;
     for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

       // Calculate for each integration point the LocPointMapped
       lp.Set( intPoints[i], esm );

       // Call the CalcBMat()-method
       this->bOperator_.CalcOpMat( bMat, lp, ptFe);

       // Calculate scalar factor
       this->scalCoef_->GetScalar(fac, lp);
       fac *= MAT_DATA_TYPE(lp.jacDet * weights[i]); 
       this->bOperator_.TransformJacDet(fac,lp,ptFe);

#ifdef USE_BLAS_VERSION
       bMat.Mult_Blas(bMat, elemMat, true, false, this->factor_*fac, 1.0);
#else
       elemMat += Transpose(bMat) * bMat * this->factor_ * fac;
#endif

     }
     ptFe->SetOnlyLowestOrder(false);
   }
}
