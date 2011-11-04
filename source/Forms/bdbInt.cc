// =====================================================================================
// 
//       Filename:  bdbInt.cc
// 
//    Description:  Implementation file for bdbIntegrators
//                  TAKE care:
//                  This file is just for inclusion in the header file!
// 
//        Version:  1.0
//        Created:  11/03/2011 07:36:28 PM
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
            class MAT_DATA_TYPE,
            class COEF_DATA_TYPE>
  BDBInt<B_OP,FE_TYPE,MAT_DATA_TYPE,COEF_DATA_TYPE>::BDBInt(shared_ptr<CoefFunction> dData, MAT_DATA_TYPE factor){
      name_ = "BDBInt";
      isSymmetric_ = true;
      assert(dData->GetType() == CoefFunction::TENSOR);
#ifndef NDEBUG
      if(dData->GetType() != CoefFunction::TENSOR){
        Exception("BDB integrator expects the coefficient function to be tensorial");
      }
#endif
      dData_ = dData;
      factor_ = factor;
  }

  //! Destructor
  template<template<class,class> class B_OP,
              class FE_TYPE,
              class MAT_DATA_TYPE,
              class COEF_DATA_TYPE>
  BDBInt<B_OP,FE_TYPE,MAT_DATA_TYPE,COEF_DATA_TYPE>::~BDBInt(){

  }

  //! Calculate the
  template<template<class,class> class B_OP,
              class FE_TYPE,
              class MAT_DATA_TYPE,
              class COEF_DATA_TYPE>
  void BDBInt<B_OP,FE_TYPE,MAT_DATA_TYPE,COEF_DATA_TYPE>::CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                                  EntityIterator& ent1,
                                  EntityIterator& ent2) {

    // Extract physical element
    const Elem* ptElem = ent1.GetElem();

    operator_.SetSolDim(Bdim_);

    Matrix<MAT_DATA_TYPE> bMat, dbMat;
    Matrix<COEF_DATA_TYPE> dMat;
    MAT_DATA_TYPE fac = 0.0;

    // Obtain FE element from feSpace and integration scheme
    shared_ptr<IntScheme> intScheme;
    FE_TYPE* ptFe = static_cast<FE_TYPE*>(ptFeSpace1_->GetFe( ent1, intScheme ));

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
    bdbMat.Resize(nrFncs * Bdim_);
#endif


    // Loop over all integration points
    LocPointMapped lp;
    for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

      // Calculate for each integration point the LocPointMapped
      lp.Set( intPoints[i], esm );

      // Call the CalcBMat()-method
      operator_.CalcOpMat( bMat, lp, ptFe);

      // Calculate D-Mat
      //calcDMat(dMat, ent1.GetElem());
      dData_->GetTensor(dMat,lp,ptElem);
      //ptMaterial_->GetTensor(dMat,ELEC_PERMITTIVITY,Global::REAL,PLANE_STRAIN);

      fac = MAT_DATA_TYPE(lp.jacDet * weights[i]);
      operator_.TransformJacDet(fac,lp,ptFe);

      dbMat.Resize(dMat.GetNumRows(),nrFncs*Bdim_);
      dbMat.Init();

#ifdef USE_BLAS_VERSION
      dMat.Mult_Blas(bMat,dbMat,false,false,1.0,0);
      dbMat = dbMat * fac;
      //Transpose(bMat).Mult(dbMat,bdbMat);
      bMat.Mult_Blas(dbMat,bdbMat,true,false,1.0,0);
      elemMat += bdbMat * factor_;

#else
      dbMat = (dMat * bMat) * fac;
      elemMat += Transpose(bMat) * dbMat * factor_;
#endif

    }
  }

}
