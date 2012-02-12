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

template< class B_OP, class VEC_DATA_TYPE >
BUIntegrator<B_OP,VEC_DATA_TYPE>::
BUIntegrator(VEC_DATA_TYPE factor,
             shared_ptr<CoefFunction > rhsCoef, 
             bool coordUpdate )
             :LinearForm( coordUpdate ){
  factor_ = factor;
  this->name_ = "RhsBUIntegrator";

         assert(rhsCoef->GetDimType() == CoefFunction::VECTOR);
 #ifndef NDEBUG
       if(rhsCoef->GetDimType() != CoefFunction::VECTOR){
         Exception("BDB integrator expects the coefficient function to be vectorial!\n \
                     For scalar valued Things, create a vectorial function with one component");
       }
 #endif
         this->rhsCoefs_ = rhsCoef;

    }

  template<class B_OP, class VEC_DATA_TYPE >
  void BUIntegrator<B_OP,VEC_DATA_TYPE>::CalcElemVector( Vector<VEC_DATA_TYPE> & elemVec,
                                                                   EntityIterator& ent){
     // Declare necessary variables
     const Elem* ptElem = ent.GetElem();
     Matrix<VEC_DATA_TYPE> bMat;
     Vector<VEC_DATA_TYPE> cVec;
     shared_ptr<IntScheme> intScheme;
     StdVector<LocPoint> intPoints;
     StdVector<Double> weights;
     UInt nrFncs = 0;
     VEC_DATA_TYPE fac;

     // Obtain FE element from feSpace and integration scheme
     BaseFE * ptFe = ptFeSpace_->GetFe( ent, intScheme );

     nrFncs = ptFe->GetNumFncs();

     // Get shape map from grid
     shared_ptr<ElemShapeMap> esm = 
         domain->GetGrid()->GetElemShapeMap( ptElem, this->coordUpdate_ );

     // Get integration points
     intScheme->GetIntPoints( Elem::GetShapeType(ptElem->type), intPoints, weights );

     elemVec.Resize( nrFncs * Bdim_);
     elemVec.Init();

     // Loop over all integration points
     LocPointMapped lp;
     for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

       // Calculate for each integration point the LocPointMapped
       lp.Set( intPoints[i], esm );

       //calc factor
       fac = VEC_DATA_TYPE(lp.jacDet * weights[i]);
       fac *= factor_;
       // Call the CalcBMat()-method
       operator_.CalcOpMatTransposed( bMat, lp, ptFe);
       
       rhsCoefs_->GetVector(cVec,lp);
       cVec *= fac;
       elemVec += bMat * cVec;
     }

  }
}
