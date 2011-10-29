// =====================================================================================
// 
//       Filename:  bdbInt_NEW.hh
// 
//    Description:  New implementation of the BDB integrator class
//                  Takes as a template parameter the operator it should evaulate
//                  new implementation to avoid old structures
// 
//        Version:  1.0
//        Created:  10/04/2011 09:28:38 AM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#ifndef FILE_BDBINT_NEW
#define FILE_BDBINT_NEW

#include "integrator.hh"
#include "Elements/basefe.hh"
#include <boost/tr1/type_traits.hpp>

namespace CoupledField {



  //! general class for calculation of bdb forms
  template<template<class,class> class B_OP,
            class FE_TYPE,
            class MAT_DATA_TYPE=Double>
  class BDBInt_NEW : public Integrator {
    public:

      //! Constructor with pointer to BaseElem
      BDBInt_NEW(shared_ptr<CoefFunction> dData, MAT_DATA_TYPE factor);

      //! Destructor
      virtual ~BDBInt_NEW();

      //! Compute element matrix associated to BDB form
      void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                                 EntityIterator& ent1,
                                 EntityIterator& ent2 );


      bool IsComplex(){
        return std::tr1::is_same<MAT_DATA_TYPE,Complex>::value;
      }

    protected:
      B_OP<FE_TYPE,MAT_DATA_TYPE> operator_;


      //! set a constant factor for multiplication with the element matrix
      MAT_DATA_TYPE factor_;
  };

  //=========================================================================================
  //Implementation of BDB Int
  //=========================================================================================

  template<template<class,class> class B_OP,
            class FE_TYPE,
            class MAT_DATA_TYPE>
  BDBInt_NEW<B_OP,FE_TYPE,MAT_DATA_TYPE>::BDBInt_NEW(shared_ptr<CoefFunction> dData, MAT_DATA_TYPE factor){
      name_ = "BDBInt";
      isSymmetric_ = true;
      dData_ = dData;
      factor_ = factor;
  }

  //! Destructor
  template<template<class,class> class B_OP,
  class FE_TYPE,
  class MAT_DATA_TYPE>
  BDBInt_NEW<B_OP,FE_TYPE,MAT_DATA_TYPE>::~BDBInt_NEW(){

  }

  //! Calculate the
  template<template<class,class> class B_OP,
            class FE_TYPE,
            class MAT_DATA_TYPE>
  void BDBInt_NEW<B_OP,FE_TYPE,MAT_DATA_TYPE>::CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                                  EntityIterator& ent1,
                                  EntityIterator& ent2) {

    // Extract physical element
    const Elem* ptElem = ent1.GetElem();

    operator_.SetSolDim(Bdim_);

    Matrix<MAT_DATA_TYPE> dMat, bMat, dbMat;
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
      dData_->GetMatrix(dMat,lp,ptElem);
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
#endif
