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


   template<template<class,class> class B_OP,
             class FE_TYPE,
             class MAT_DATA_TYPE>
   void BBInt<B_OP,FE_TYPE,MAT_DATA_TYPE>::CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                                   EntityIterator& ent1,
                                   EntityIterator& ent2) {
       // Extract physical element
          const Elem* ptElem = ent1.GetElem();

          operator_.SetSolDim(Bdim_);

          Matrix<MAT_DATA_TYPE> bMat,bbMat;
          MAT_DATA_TYPE fac = 0.0;

          // Obtain FE element from feSpace and integration scheme
          shared_ptr<IntScheme> intScheme;
          FE_TYPE* ptFe = dynamic_cast<FE_TYPE*>(ptFeSpace1_->GetFe( ent1, intScheme ));

          UInt nrFncs = ptFe->GetNumFncs();

          // Get shape map from grid
          shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap( ptElem );

          // Get integration points
          StdVector<LocPoint> intPoints;
          StdVector<Double> weights;
          intScheme->GetIntPoints( Elem::GetShapeType(ptElem->type), intPoints, weights );

          elemMat.Resize( nrFncs * Bdim_);
          elemMat.Init();

       #define USE_BLAS_VERSION
      #ifdef USE_BLAS_VERSION
          Matrix<MAT_DATA_TYPE> bdbMat;
          bbMat.Resize(nrFncs * Bdim_);
      #endif
          // Loop over all integration points
          LocPointMapped lp;
          for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

            // Calculate for each integration point the LocPointMapped
            lp.Set( intPoints[i], esm );

            // Call the CalcBMat()-method
            operator_.CalcOpMat( bMat, lp, ptFe);

            fac = MAT_DATA_TYPE(lp.jacDet * weights[i]);
            operator_.TransformJacDet(fac,lp,ptFe);

      #ifdef USE_BLAS_VERSION
            bMat.Mult_Blas(bMat,bbMat,true,false,1.0,0);
            elemMat += bdbMat * factor_;
      #else
            elemMat += Transpose(bMat) * bMat * factor_;
      #endif

     }
   }

   template<template<class,class> class B_OP,
               class FE_TYPE,
               class MAT_DATA_TYPE>
   void BBIntMassEdge<B_OP,FE_TYPE,MAT_DATA_TYPE>::CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                                                                               EntityIterator& ent1,
                                                                               EntityIterator& ent2 ) {
       // Extract physical element
       const Elem* ptElem = ent1.GetElem();

       this->operator_.SetSolDim(this->Bdim_);

       Matrix<MAT_DATA_TYPE> bMat,bbMat;
       MAT_DATA_TYPE fac = 0.0;

       // Obtain FE element from feSpace and integration scheme
       shared_ptr<IntScheme> intScheme;
       //TODO: we have to find another solution for this. The only reason for this class is within
       //these two lines. This functionaliy has to be moved to a pre-Integration function
       FeHCurlHi* ptFe = dynamic_cast<FeHCurlHi*>(this->ptFeSpace1_->GetFe( ent1, intScheme ));
       // Special: Only use lower order functions
       ptFe->SetOnlyLowestOrder(true);

       UInt nrFncs = ptFe->BaseFE::GetNumFncs();

       // Get shape map from grid
       shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap( ptElem );

       // Get integration points
       StdVector<LocPoint> intPoints;
       StdVector<Double> weights;
       intScheme->GetIntPoints( Elem::GetShapeType(ptElem->type), intPoints, weights );

       elemMat.Resize( nrFncs * this->Bdim_);
       elemMat.Init();

   #define USE_BLAS_VERSION
      #ifdef USE_BLAS_VERSION
        Matrix<MAT_DATA_TYPE> bdbMat;
        bbMat.Resize(nrFncs * this->Bdim_);
     #endif
       // Loop over all integration points
       LocPointMapped lp;
       for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

         // Calculate for each integration point the LocPointMapped
         lp.Set( intPoints[i], esm );

         // Call the CalcBMat()-method
         this->operator_.CalcOpMat( bMat, lp, ptFe);

         fac = MAT_DATA_TYPE(lp.jacDet * weights[i]);
         this->operator_.TransformJacDet(fac,lp,ptFe);

   #ifdef USE_BLAS_VERSION
         bMat.Mult_Blas(bMat,bbMat,true,false,1.0,0);
         elemMat += bdbMat * this->factor_;
   #else
         elemMat += Transpose(bMat) * bMat * factor_;
   #endif

       }
      ptFe->SetOnlyLowestOrder(false);
   }
}
