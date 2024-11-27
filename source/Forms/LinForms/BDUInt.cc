// =====================================================================================
// 
//       Filename:  bduInt.cc
// 
//    Description:  Implementation file for BDUIntegrators
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
             bool coordUpdate,
			 bool extractReal)
             :LinearForm( coordUpdate ){
  factor_ = factor;
  this->name_ = "RhsBDUIntegrator";
  Bdim_ = B_OP::DIM_SPACE;

//  assert(rhsCoef->GetDimType() == CoefFunction::VECTOR);
//#ifndef NDEBUG
//  if(rhsCoef->GetDimType() != CoefFunction::VECTOR){
//    Exception("BDB integrator expects the coefficient function to be vectorial!\n \
//                     For scalar valued Things, create a vectorial function with one component");
//  }
//#endif
  this->rhsCoefs_ = rhsCoef;
  this->dCoef_ = dCoef;
  extractReal_ = extractReal;
}

template< class B_OP, class VEC_DATA_TYPE, bool SURFACE >
BDUIntegrator<B_OP,VEC_DATA_TYPE,SURFACE>::
BDUIntegrator(VEC_DATA_TYPE factor,
             shared_ptr<CoefFunction > rhsCoef,
             shared_ptr<CoefFunction > dCoef,
             const std::set<RegionIdType>& volRegions,
             bool coordUpdate,
			 bool extractReal)
             :LinearForm( coordUpdate ){
  factor_ = factor;
  this->name_ = "RhsBDUIntegrator";
  Bdim_ = B_OP::DIM_SPACE;

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
  extractReal_ = extractReal;
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
     VEC_DATA_TYPE fac(0.0);

     // Obtain FE element from feSpace and integration scheme
     IntegOrder order;
     IntScheme::IntegMethod method;
     BaseFE* ptFe = ptFeSpace_->GetFe( ent, method, order );

     nrFncs = ptFe->GetNumFncs();

     // Get shape map from grid
    //  shared_ptr<ElemShapeMap> esm =
    //      ent.GetGrid()->GetElemShapeMap( ptElem, this->coordUpdate_ );
    // shared_ptr<ElemShapeMap> esm(ptElem->ptrShapeMap);
    shared_ptr<ElemShapeMap> esm = (ptElem)->GetElemShapeMap((ent.GetGrid()), this->coordUpdate_);

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
         lp.SetWithSurface( intPoints[i], esm, volRegions_, weights[i] );
       } else {
         lp.Set( intPoints[i], esm, weights[i] );
       }

//       std::cout << "integration point: " << lp.lp.coord.ToString() << std::endl;

       // obtain d matrix
       this->dCoef_->GetTensor( dMat, lp );
       
//       std::cout << "\n elem: " << ptElem->elemNum << std::endl;
//       std::cout << "dMat: \n" << dMat.ToString(0) << std::endl;
       
       //calc factor
       fac = VEC_DATA_TYPE(lp.jacDet * weights[i]);
       fac *= factor_;
       
       // Call the CalcBMat()-method
       operator_.CalcOpMatTransposed( bMat, lp, ptFe);
       
//       Matrix<VEC_DATA_TYPE> tmp;
//       bMat.Transpose(tmp);
//       std::cout << "bMat:\n" << tmp.ToString(0) << std::endl;
       
       assert(bMat.GetNumRows() == nrFncs * B_OP::DIM_DOF);
       bdMat.Resize(bMat.GetNumRows(), dMat.GetNumCols());
       
       // Calculate BdMat
#ifdef NDEBUG
       bMat.Mult_Blas(dMat, bdMat,false,false,fac,0.0); // bdMat = fac * bMat * dMat + 0.0 * bdMat
#else
       bdMat = (bMat * dMat) * fac;
#endif

//       std::cout << "bdMat:\n" << bdMat.ToString(0) << std::endl;

       rhsCoefs_->GetVector(cVec,lp);
       //elemVec += bMat * cVec * fac;
//       std::cout << " cVec:" << cVec.ToString() << std::endl;

       elemVec += bdMat * cVec;
//       std::cout << " elemVec:" << elemVec.ToString() << std::endl;
     }

  }
}
