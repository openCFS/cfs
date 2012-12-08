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


   template< class B_OP, class MAT_DATA_TYPE, class COEF_DATA_TYPE> 
   BBInt<B_OP,MAT_DATA_TYPE,COEF_DATA_TYPE>::
   BBInt(PtrCoefFct scalCoef, MAT_DATA_TYPE factor, 
         bool coordUpdate) :
   BaseBDBInt(coordUpdate) {
     this->name_ = "BBInt";
     this->isSymmetric_ = true;
     this->coefScalar_ = scalCoef;
     this->factor_ = factor;
     
     // Ensure, that the coefficient set is a scalar valued one
     if( this->coefScalar_->GetDimType() != CoefFunction::SCALAR ) {
       EXCEPTION( "The BBInt-class only works with scalar-valued "
           << "coefficient functions!");
     }
   }

   template< class B_OP, class MAT_DATA_TYPE, class COEF_DATA_TYPE> 
   void BBInt<B_OP,MAT_DATA_TYPE,COEF_DATA_TYPE>::
   CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                      EntityIterator& ent1,
                      EntityIterator& ent2) {
     // Extract physical element
     const Elem* ptElem = ent1.GetElem();

     Matrix<MAT_DATA_TYPE> bMat;
     MAT_DATA_TYPE fac = 0.0;

     // Obtain FE element from feSpace and integration scheme
     IntegOrder order;
      IntScheme::IntegMethod method;
      BaseFE* ptFe = ptFeSpace1_->GetFe( ent1, method, order );
     UInt nrFncs = ptFe->GetNumFncs();

     // Get shape map from grid
     shared_ptr<ElemShapeMap> esm = 
         domain->GetGrid()->GetElemShapeMap( ptElem, this->coordUpdate_ );

     // Get integration points
     StdVector<LocPoint> intPoints;
     StdVector<Double> weights;
     intScheme_->GetIntPoints( Elem::GetShapeType(ptElem->type), method, order,
                               intPoints, weights );

     elemMat.Resize( nrFncs * B_OP::DIM_DOF );
     elemMat.Init();

#define USE_BLAS_VERSION
     // Loop over all integration points
     LocPointMapped lp;
     for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

       // Calculate for each integration point the LocPointMapped
       lp.Set( intPoints[i], esm );

       // Call the CalcBMat()-method
       this->bOperator_.CalcOpMat( bMat, lp, ptFe);

       // Calculate scalar factor
       this->coefScalar_->GetScalar(fac, lp);
       fac *= MAT_DATA_TYPE(lp.jacDet * weights[i]); 

       this->bOperator_.TransformJacDet(fac,lp,ptFe);

#ifdef USE_BLAS_VERSION
       bMat.Mult_Blas(bMat, elemMat, true, false, this->factor_ * fac, 1.0);
#else
       elemMat += Transpose(bMat) * bMat * this->factor_ * fac;
#endif

     }
   }
   
   //! Multiply element matrix with vector
   template< class B_OP, class MAT_DATA_TYPE, class COEF_DATA_TYPE> 
   void BBInt<B_OP,MAT_DATA_TYPE,COEF_DATA_TYPE>::
   ApplyElemMat( Vector<MAT_DATA_TYPE>&ret, const Vector<Double>& sol,
                 EntityIterator& ent1,
                 EntityIterator& ent2 ) {
     Matrix<MAT_DATA_TYPE> elemMat;
     CalcElementMatrix(elemMat, ent1, ent2);
     ret = elemMat * sol;
   }

   //! Apply B-operator on vector
   template< class B_OP, class MAT_DATA_TYPE, class COEF_DATA_TYPE> 
   void BBInt<B_OP,MAT_DATA_TYPE,COEF_DATA_TYPE>::
   ApplyBMat( Vector<MAT_DATA_TYPE>&ret, 
              const Vector<MAT_DATA_TYPE>& sol,
              const LocPointMapped& lpm ) {
     Matrix<MAT_DATA_TYPE> bOp;
     BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
     bOperator_.CalcOpMat(bOp, lpm, ptFe);
     ret = bOp * sol;
   }

   //! Apply dB-operator on vector
   template< class B_OP, class MAT_DATA_TYPE, class COEF_DATA_TYPE> 
   void BBInt<B_OP,MAT_DATA_TYPE,COEF_DATA_TYPE>::
   ApplydBMat( Vector<MAT_DATA_TYPE>&ret, 
               const Vector<MAT_DATA_TYPE>& sol,
               const LocPointMapped& lpm ) {
     Matrix<MAT_DATA_TYPE> bMat;
     COEF_DATA_TYPE fac;
     BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
     bOperator_.CalcOpMat(bMat, lpm, ptFe);
     this->coefScalar_->GetScalar(fac, lpm);

     ret = (bMat * sol) * fac;
   }

   template< class B_OP, class MAT_DATA_TYPE, class COEF_DATA_TYPE> 
   void BBInt<B_OP,MAT_DATA_TYPE,COEF_DATA_TYPE>::
   CalcKernel( Matrix<MAT_DATA_TYPE>& kernel, 
               const LocPointMapped& lpm ) {

     Matrix<MAT_DATA_TYPE> bMat;
     MAT_DATA_TYPE fac = 0.0;

     // Obtain FE element from feSpace and integration scheme
     BaseFE* ptFe = this->ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
     UInt nrFncs = ptFe->GetNumFncs();

     kernel.Resize( nrFncs * B_OP::DIM_DOF );
     kernel.Init();

#define USE_BLAS_VERSION

     // Call the CalcBMat()-method
     this->bOperator_.CalcOpMat( bMat, lpm, ptFe);

     // Calculate scalar factor
     this->coefScalar_->GetScalar(fac, lpm);
     this->bOperator_.TransformJacDet(fac, lpm, ptFe);

#ifdef USE_BLAS_VERSION
     bMat.Mult_Blas(bMat, kernel, true, false, this->factor_ * fac, 1.0);
#else
     elemMat += Transpose(bMat) * bMat * this->factor_ * fac;
#endif
   }
   
   
   // =======================================================================
   
   
   template<class B_OP, class MAT_DATA_TYPE>
     void BBIntMassEdge<B_OP,MAT_DATA_TYPE>::
   CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                      EntityIterator& ent1,
                      EntityIterator& ent2 ) {
     // Extract physical element
     const Elem* ptElem = ent1.GetElem();

     Matrix<MAT_DATA_TYPE> bMat,bbMat;
     MAT_DATA_TYPE fac = 0.0;

     // Obtain FE element from feSpace and integration scheme
     IntegOrder order;
     IntScheme::IntegMethod method;
     //TODO: we have to find another solution for this. The only reason for this class is within
     //these two lines. This functionaliy has to be moved to a pre-Integration function
     FeHCurlHi* ptFe = 
         static_cast<FeHCurlHi*>(this->ptFeSpace1_->GetFe( ent1, method, order ));
     // Special: Only use lower order functions
     ptFe->SetOnlyLowestOrder(true);

     UInt nrFncs = ptFe->BaseFE::GetNumFncs();

     // Get shape map from grid
     shared_ptr<ElemShapeMap> esm = 
         domain->GetGrid()->GetElemShapeMap( ptElem, this->coordUpdate_ );

     // Get integration points
     StdVector<LocPoint> intPoints;
     StdVector<Double> weights;
     this->intScheme_->GetIntPoints( Elem::GetShapeType(ptElem->type), method, order,
                                     intPoints, weights );

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
       this->coefScalar_->GetScalar(fac, lp);
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

   template<class B_OP, class MAT_DATA_TYPE,class COEF_DATA_TYPE>
     void SurfaceBBInt<B_OP,MAT_DATA_TYPE,COEF_DATA_TYPE>::
   CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                      EntityIterator& ent1,
                      EntityIterator& ent2 ) {
     // Extract physical element
     const Elem* ptElem1 = ent1.GetElem();
     const Elem* ptElem2 = ent2.GetElem();

     Matrix<MAT_DATA_TYPE> bMat,bMatT;
     MAT_DATA_TYPE fac = 0.0;

     // Obtain FE element from feSpace and integration scheme
     IntegOrder order1, order2;
     IntScheme::IntegMethod method1, method2;
     BaseFE* ptFe1 = this->ptFeSpace1_->GetFe( ent1, method1, order1 );
     BaseFE* ptFe2 = this->ptFeSpace1_->GetFe( ent2, method2, order2 );
     UInt nrFncs1 = ptFe1->GetNumFncs();


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

     elemMat.Resize( nrFncs1 * B_OP::DIM_DOF );
     elemMat.Init();

#define USE_BLAS_VERSION
     // Loop over all integration points
     LocPointMapped lp1,lp2;
     for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

       // Calculate for each integration point the LocPointMapped
       lp1.Set( intPoints[i], esm1, volRegions_ );
       lp2.Set( intPoints[i], esm2, volRegions_ );

       // Call the CalcBMat()-method
       this->bOperator_.CalcOpMatTransposed( bMatT, lp1, ptFe1);
       this->bOperator_.CalcOpMat( bMat, lp2, ptFe2);

       // Calculate scalar factor
       //TODO: Which point to take? lp1 or lp2?
       // porposal: let PDE decide via constructor parameter
       this->coefScalar_->GetScalar(fac, lp1);
       fac *= MAT_DATA_TYPE(lp1.jacDet * weights[i]);

       this->bOperator_.TransformJacDet(fac,lp1,ptFe1);

#ifdef USE_BLAS_VERSION
       bMatT.Mult_Blas(bMat, elemMat, false, false, this->factor_ * fac, 1.0);
#else
       elemMat += Transpose(bMat) * bMat * this->factor_ * fac;
#endif

     }
   }


}
