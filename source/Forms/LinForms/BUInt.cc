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
             bool fullEvaluation,
             bool extractReal,
             const string& id)
             : LinearForm( coordUpdate ),
               fullEvaluation_(fullEvaluation) {
  factor_ = factor;
  this->name_ = "RhsBUIntegrator";
  this->bOperator_= bOp;

  this->rhsCoefs_ = rhsCoef;
  extractReal_ = extractReal;
  id_ = id;
}

template< class VEC_DATA_TYPE, bool SURFACE >
BUIntegrator<VEC_DATA_TYPE,SURFACE>::
BUIntegrator(BaseBOperator * bOp,
             VEC_DATA_TYPE factor,
             shared_ptr<CoefFunction > rhsCoef,
             const std::set<RegionIdType>& volRegions,
             bool coordUpdate,
             bool fullEvaluation,
             bool extractReal,
             const string& id)
             : LinearForm( coordUpdate ), 
               fullEvaluation_(fullEvaluation) 
               {
  factor_ = factor;
  this->name_ = "RhsBUIntegrator";
  this->bOperator_= bOp;
  extractReal_ = extractReal;
  id_ = id;

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
  VEC_DATA_TYPE fac(0.0);

  //Surface: inverse of jacobian
  Vector<VEC_DATA_TYPE> pt1;
  Vector<VEC_DATA_TYPE> pt2;
  Matrix<Double> JacT;
  Matrix<Double> TF;
  Matrix<Double> TFinv;

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
      if ( scalVec_.GetSize() > 0 ) {
        for (UInt k=0; k<scalVec_.GetSize(); k++) {
            cVec[k] *= scalVec_[k];
        }
      }
    }
  }

     // Loop over all integration points
     VEC_DATA_TYPE vol = 0.0; //for volume normalization
     for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

       // Calculate for each integration point the LocPointMapped
       if (SURFACE) {
         lp.SetWithSurface( intPoints[i], esm, volRegions_, weights[i] );
       } else {
         lp.Set( intPoints[i], esm, weights[i] );
       }

       //calc factor
       fac = VEC_DATA_TYPE(lp.jacDet * weights[i]);
       vol += fac; //computes the volume
       
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
        if ( scalVec_.GetSize() > 0 ) {
          for (UInt k=0; k<scalVec_.GetSize(); k++)
            cVec[k] *= scalVec_[k];             
        }
        if (SURFACE && (ptFeSpace_->GetSpaceType() == FeSpace::HCURL)) {
          //uxn
          pt1 = lp.normal;
          cVec.CrossProduct(pt1,pt2);
          // Jacobian of surface element
          lp.jac.Transpose(JacT);
          //Metric and its inverse
          TF = JacT * lp.jac;
          TF.Invert(TFinv);

          // Transformation of a function in curl space (see Zaglmayer Lemma 4.15)
          cVec = (TFinv * JacT) * pt2;
          fac *= lp.lpmVol->jacDet;
        }
      }
    }  
    elemVec += bMat * cVec * fac;    
  }  
  //normalize to volume
  if ( normalizeToVol_) {
    elemVec /= vol;    
  }
}

}
