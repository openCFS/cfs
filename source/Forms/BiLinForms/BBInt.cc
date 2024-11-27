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

#include "BBInt.hh"
#include "Domain/Domain.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

DEFINE_LOG(bbint, "bbint")

namespace CoupledField{


   template< class COEF_DATA_TYPE, class B_DATA_TYPE> 
   BBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
   BBInt(BaseBOperator * bOp,
         PtrCoefFct scalCoef, MAT_DATA_TYPE factor, 
         bool coordUpdate) :
   BaseBDBInt(coordUpdate) {
     this->type_ = BB_INT;
     this->name_ = type.ToString(type_);
     this->isSymmetric_ = true;
     this->bOperator_ = bOp;
     this->coefScalar_ = scalCoef;
     this->factor_ = factor;

     // Ensure, that the coefficient set is a scalar valued one
     if( this->coefScalar_->GetDimType() != CoefFunction::SCALAR ) {
       EXCEPTION( "The BBInt-class only works with scalar-valued "
           << "coefficient functions!");
     }
   }

   template< class COEF_DATA_TYPE, class B_DATA_TYPE> 
   void BBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
   CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                      EntityIterator& ent1,
                      EntityIterator& ent2) {
     IntegOrder order;
     IntScheme::IntegMethod method;
     StdVector<LocPoint> intPoints;
     StdVector<Double> weights;
     Matrix<MAT_DATA_TYPE> bMat;
     MAT_DATA_TYPE fac = 0.0;

     // Extract physical element
     const Elem* ptElem = ent1.GetElem();
     BaseFE* ptFe = ptFeSpace1_->GetFe( ent1, method, order );
     const UInt nrFncs = ptFe->GetNumFncs();

     // Get shape map from grid
    //  shared_ptr<ElemShapeMap> esm = ent1.GetGrid()->GetElemShapeMap( ptElem, this->coordUpdate_ );
    // shared_ptr<ElemShapeMap> esm(ptElem->ptrShapeMap);
    shared_ptr<ElemShapeMap> esm = (ptElem)->GetElemShapeMap(ent1.GetGrid(), this->coordUpdate_);

     // Get integration points
     intScheme_->GetIntPoints( Elem::GetShapeType(ptElem->type), method, order,
                               intPoints, weights );

     elemMat.Resize( nrFncs * bOperator_->GetDimDof() );
     elemMat.Init();

     // Loop over all integration points
     LocPointMapped lp;
     const UInt numIntPts = intPoints.GetSize();
     for( UInt i = 0; i < numIntPts; ++i ) {

       // Calculate for each integration point the LocPointMapped
       lp.Set( intPoints[i], esm, weights[i] );

       // Call the CalcBMat()-method
       this->bOperator_->CalcOpMat( bMat, lp, ptFe);
       LOG_DBG3(bbint) << "e= " << ptElem->elemNum << " bMat= " << bMat;

       // Calculate scalar factor
       this->coefScalar_->GetScalar(fac, lp);
       LOG_DBG3(bbint) << "e= " << ptElem->elemNum << " nu= " << fac;

       fac *= MAT_DATA_TYPE(lp.jacDet * weights[i]); 

#ifdef NDEBUG
       bMat.Mult_Blas(bMat, elemMat, true, false, this->factor_ * fac, 1.0);
#else
       elemMat += Transpose(bMat) * bMat * this->factor_ * fac;
       LOG_DBG3(bbint) << "CEM e=" << ptElem->elemNum << " ip=" << i << " fac=" << fac << " factor_=" << fac << " bmat=" << bMat.ToString();
       LOG_DBG3(bbint) << "CEM e=" << ptElem->elemNum << " -> K_" << i << "=" << elemMat.ToString();
#endif
     }
   }
   
   // ===============
   //  Apply ElemMat
   // ===============
   // 1) General template for double-case -> not implemented
   template< class COEF_DATA_TYPE, class B_DATA_TYPE> 
   void BBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
   ApplyElemMat( Vector<Double>&ret, const Vector<Double>& sol,
                 EntityIterator& ent1,
                 EntityIterator& ent2 ) {
     EXCEPTION("Not implemented");
   }
   
   // 2) Special double-only case
   template<> 
   void BBInt<Double, Double>::
   ApplyElemMat( Vector<Double>&ret, 
                 const Vector<Double>& sol,
                 EntityIterator& ent1,
                 EntityIterator& ent2 ) {
     Matrix<MAT_DATA_TYPE> elemMat;
     CalcElementMatrix(elemMat, ent1, ent2);
     ret = elemMat * sol;
   }

   // 3) Complex case is always possible
   template< class COEF_DATA_TYPE, class B_DATA_TYPE> 
   void BBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
   ApplyElemMat( Vector<Complex>&ret, 
                 const Vector<Complex>& sol,
                 EntityIterator& ent1,
                 EntityIterator& ent2 ) {
     Matrix<MAT_DATA_TYPE> elemMat;
     CalcElementMatrix(elemMat, ent1, ent2);
     ret = elemMat * sol;
   }
   
   // ===============
   //  Apply BMatrix
   // ===============
   // 1) General template for double-case -> not implemented
   template< class COEF_DATA_TYPE, class B_DATA_TYPE> 
   void BBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
   ApplyBMat( Vector<Double>&ret, 
              const Vector<Double>& sol,
              const LocPointMapped& lpm ) {
     EXCEPTION( "Not implemented" );
   }

   // 2) Special double-only case
   template<> 
   void BBInt<Double, Double>::
   ApplyBMat( Vector<Double>&ret, 
              const Vector<Double>& sol,
              const LocPointMapped& lpm ) {
     BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
     bOperator_->CalcOpMat(bMat_, lpm, ptFe);
     ret = bMat_ * sol;
   }
   
   // 3) Complex case is always possible
   template< class COEF_DATA_TYPE, class B_DATA_TYPE> 
   void BBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
   ApplyBMat( Vector<Complex>&ret, 
              const Vector<Complex>& sol,
              const LocPointMapped& lpm ) {
     BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
     bOperator_->CalcOpMat(bMat_, lpm, ptFe);
     ret = bMat_ * sol;
   }

   // ================
   //  Apply dBMatrix
   // ================   
   // 1) General template for double-case -> not implemented
    template< class COEF_DATA_TYPE, class B_DATA_TYPE> 
    void BBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
    ApplydBMat( Vector<Double>&ret, 
               const Vector<Double>& sol,
               const LocPointMapped& lpm ) {
      Exception( "Not implemented" );
   }
    
    template<> 
    void BBInt<Double, Double>::
    ApplydBMat( Vector<Double>&ret, 
                const Vector<Double>& sol,
                const LocPointMapped& lpm ) {
      Double fac;
      BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
      bOperator_->CalcOpMat(bMat_, lpm, ptFe);
      this->coefScalar_->GetScalar(fac, lpm);

      ret = (bMat_ * sol) * fac;
    }

    template< class COEF_DATA_TYPE, class B_DATA_TYPE> 
    void BBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
    ApplydBMat( Vector<Complex>&ret, 
                const Vector<Complex>& sol,
                const LocPointMapped& lpm ) {
      COEF_DATA_TYPE fac;
      BaseFE* ptFe = ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
      bOperator_->CalcOpMat(bMat_, lpm, ptFe);
      this->coefScalar_->GetScalar(fac, lpm);

      ret = (bMat_ * sol) * Complex(fac);
    }

    
   template< class COEF_DATA_TYPE, class B_DATA_TYPE> 
   void BBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
   CalcKernel( Matrix<MAT_DATA_TYPE>& kernel, 
               const LocPointMapped& lpm ) {

     MAT_DATA_TYPE fac = 0.0;

     // Obtain FE element from feSpace and integration scheme
     BaseFE* ptFe = this->ptFeSpace1_->GetFe( lpm.ptEl->elemNum );
     const UInt nrFncs = ptFe->GetNumFncs();

     kernel.Resize( nrFncs * bOperator_->GetDimDof() );
     kernel.Init();


     // Call the CalcBMat()-method
     this->bOperator_->CalcOpMat( bMat_, lpm, ptFe);

     // Calculate scalar factor
     this->coefScalar_->GetScalar(fac, lpm);

#ifdef NDEBUG
     bMat_.Mult_Blas(bMat_, kernel, true, false, this->factor_ * fac, 1.0);
#else
     kernel += Transpose(bMat_) * bMat_ * this->factor_ * fac;
#endif
   }
   
   
   // =======================================================================
   
   
   template< class COEF_DATA_TYPE, class B_DATA_TYPE> 
     void BBIntMassEdge<COEF_DATA_TYPE, B_DATA_TYPE>::
   CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                      EntityIterator& ent1,
                      EntityIterator& ent2 ) {
     // Extract physical element
     const Elem* ptElem = ent1.GetElem();

     MAT_DATA_TYPE fac = 0.0;

     // Obtain FE element from feSpace and integration scheme
     IntegOrder order;
     IntScheme::IntegMethod method;
     //TODO: we have to find another solution for this. The only reason for this class is within
     //these two lines. This functionaliy has to be moved to a pre-Integration function
     FeHCurlHi* ptFe = 
         static_cast<FeHCurlHi*>(this->ptFeSpace1_->GetFe( ent1, method, order ));
     // Special: Only use lower order functions
     //ptFe->SetOnlyLowestOrder(true);

     const UInt nrFncs = ptFe->GetNumFncs();

     // Get shape map from grid
    //  shared_ptr<ElemShapeMap> esm =
    //      ent1.GetGrid()->GetElemShapeMap( ptElem, this->coordUpdate_ );
    // shared_ptr<ElemShapeMap> esm(ptElem->ptrShapeMap);
    shared_ptr<ElemShapeMap> esm = (ptElem)->GetElemShapeMap(ent1.GetGrid(), this->coordUpdate_);

     // Get integration points
     StdVector<LocPoint> intPoints;
     StdVector<Double> weights;
     this->intScheme_->GetIntPoints( Elem::GetShapeType(ptElem->type), method, order,
                                     intPoints, weights );

     elemMat.Resize( nrFncs * this->bOperator_->GetDimDof() );
     elemMat.Init();

     // Loop over all integration points
     LocPointMapped lp;
     const UInt numIntPts = intPoints.GetSize();
     for( UInt i = 0; i < numIntPts; ++i  ) {

       // Calculate for each integration point the LocPointMapped
       lp.Set( intPoints[i], esm, weights[i] );

       // Call the CalcBMat()-method
       this->bOperator_->CalcOpMat( this->bMat_, lp, ptFe);

       // Calculate scalar factor
       this->coefScalar_->GetScalar(fac, lp);

       fac *= MAT_DATA_TYPE(lp.jacDet * weights[i]); 

#ifdef NDEBUG
       this->bMat_.Mult_Blas(this->bMat_, elemMat, true, false, this->factor_*fac, 1.0);
#else
       elemMat += Transpose(this->bMat_) * this->bMat_ * this->factor_ * fac;
#endif
     }
     //ptFe->SetOnlyLowestOrder(false);
   }

   template< class COEF_DATA_TYPE, class B_DATA_TYPE> 
   void SurfaceBBInt<COEF_DATA_TYPE, B_DATA_TYPE>::
   CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                      EntityIterator& ent1,
                      EntityIterator& ent2 ) {
     // Extract physical element
     const Elem* ptElem1 = ent1.GetElem();
     const Elem* ptElem2 = ent2.GetElem();

     MAT_DATA_TYPE fac = 0.0;

     // Obtain FE element from feSpace and integration scheme
     IntegOrder order1, order2;
     IntScheme::IntegMethod method1, method2;
     BaseFE* ptFe1 = this->ptFeSpace1_->GetFe( ent1, method1, order1 );
     BaseFE* ptFe2 = this->ptFeSpace1_->GetFe( ent2, method2, order2 );
     const UInt nrFncs1 = ptFe1->GetNumFncs();


     // Get shape map from grid
    //  shared_ptr<ElemShapeMap> esm1 =
    //      ent1.GetGrid()->GetElemShapeMap( ptElem1, this->coordUpdate_ );
    //  shared_ptr<ElemShapeMap> esm2 =
    //      ent2.GetGrid()->GetElemShapeMap( ptElem2, this->coordUpdate_ );
    // shared_ptr<ElemShapeMap> esm1(ptElem1->ptrShapeMap);
    // shared_ptr<ElemShapeMap> esm2(ptElem2->ptrShapeMap);
    shared_ptr<ElemShapeMap> esm1 = (ptElem1)->GetElemShapeMap(ent1.GetGrid(), this->coordUpdate_);
    shared_ptr<ElemShapeMap> esm2 = (ptElem2)->GetElemShapeMap(ent2.GetGrid(), this->coordUpdate_);

     // Get integration points
     StdVector<LocPoint> intPoints;
     StdVector<Double> weights;
     this->intScheme_->GetIntPoints( Elem::GetShapeType(ptElem1->type), 
                                     method1, order1, method2, order2,
                                     intPoints, weights );

     elemMat.Resize( nrFncs1 * this->bOperator_->GetDimDof() );
     elemMat.Init();

     // Loop over all integration points
     LocPointMapped lp1,lp2;
     const UInt numIntPts = intPoints.GetSize();
     for( UInt i = 0; i < numIntPts; ++i  ) {

       // Calculate for each integration point the LocPointMapped
       lp1.SetWithSurface( intPoints[i], esm1, volRegions_, weights[i] );
       lp2.SetWithSurface( intPoints[i], esm2, volRegions_, weights[i] );

       // Call the CalcBMat()-method
       this->bOperator_->CalcOpMatTransposed( bMatT_, lp1, ptFe1);
       this->bOperator_->CalcOpMat( this->bMat_, lp2, ptFe2);


       // Calculate scalar factor
       //TODO: Which point to take? lp1 or lp2?
       // porposal: let PDE decide via constructor parameter
       this->coefScalar_->GetScalar(fac, lp1);
       fac *= MAT_DATA_TYPE(lp1.jacDet * weights[i]);

#ifdef NDEBUG
       bMatT_.Mult_Blas(this->bMat_, elemMat, false, false, this->factor_ * fac, 1.0);
#else
       elemMat += Transpose(this->bMat_) * this->bMat_ * this->factor_ * fac;
#endif
     }
   }

   // Explicit template instantiation
  template class BBInt<Double,Double>;
  template class BBInt<Double,Complex>;
  template class BBInt<Complex,Double>;
  template class BBInt<Complex,Complex>;
  
  template class BBIntMassEdge<Double,Double>;
  template class BBIntMassEdge<Double,Complex>;
  template class BBIntMassEdge<Complex,Double>;
  template class BBIntMassEdge<Complex,Complex>;
  
  template class SurfaceBBInt<Double,Double>;
  template class SurfaceBBInt<Double,Complex>;
  template class SurfaceBBInt<Complex,Double>;
  template class SurfaceBBInt<Complex,Complex>;
}
