// =====================================================================================
// 
//       Filename:  abInt.cc
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


   template< class A_OP, class B_OP, class MAT_DATA_TYPE > 
   ABInt<A_OP, B_OP,MAT_DATA_TYPE>::
   ABInt( shared_ptr<CoefFunction> scalCoef, MAT_DATA_TYPE factor )
     : BBInt<B_OP, MAT_DATA_TYPE>(scalCoef, factor) {
     this->name_ = "ABInt";
     
     // Note: In general the AB-Integrator is not symmetric, as is should 
     // geto only used, if the A- and B differential operators are distinct.
     this->isSymmetric_ = false;
   }

   template< class A_OP, class B_OP, class MAT_DATA_TYPE > 
   void ABInt<A_OP, B_OP,MAT_DATA_TYPE>::
   CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                      EntityIterator& ent1,
                      EntityIterator& ent2 ) {
     // Extract physical element
     const Elem* ptElem = ent1.GetElem();

     Matrix<MAT_DATA_TYPE> aMat, bMat;
     MAT_DATA_TYPE fac = 0.0;

     // Obtain FE element from feSpace and integration scheme
     shared_ptr<IntScheme> intScheme;
     BaseFE* ptFeA = this->ptFeSpace1_->GetFe( ent1, intScheme );
     BaseFE* ptFeB = this->ptFeSpace2_->GetFe( ent1, intScheme );

     UInt nrFncsA = ptFeA->GetNumFncs();
     UInt nrFncsB = ptFeB->GetNumFncs();

     // Get shape map from grid
     shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap( ptElem );

     // Get integration points
     StdVector<LocPoint> intPoints;
     StdVector<Double> weights;
     intScheme->GetIntPoints( Elem::GetShapeType(ptElem->type), intPoints, weights );

     elemMat.Resize( nrFncsA * A_OP::DIM_DOF, 
                     nrFncsB * B_OP::DIM_DOF );
     elemMat.Init();

#define USE_BLAS_VERSION
     // Loop over all integration points
     LocPointMapped lp;
     for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

       // Calculate for each integration point the LocPointMapped
       lp.Set( intPoints[i], esm );

       // Calculate A-matrix (first differential operator)
       this->aOperator_.CalcOpMat( aMat, lp, ptFeA );
       
       // Calculate B-matrix (second differential operator)
       this->bOperator_.CalcOpMat( bMat, lp, ptFeB );

       // Calculate scalar factor
       this->scalCoef_->GetScalar(fac, lp);
       fac *= MAT_DATA_TYPE(lp.jacDet * weights[i]); 

       this->aOperator_.TransformJacDet(fac,lp,ptFeA);
       this->bOperator_.TransformJacDet(fac,lp,ptFeB);

#ifdef USE_BLAS_VERSION
       aMat.Mult_Blas(bMat,elemMat,true,false,this->factor_*fac,1.0);
#else
       elemMat += Transpose(aMat) * bMat * this->factor_;
#endif

     }

   }

  
} // name
