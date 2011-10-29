// =====================================================================================
// 
//       Filename:  BBInt.hh
// 
//    Description:  This class implements a general symmetric integrator without any 
//                  material factors. The other opprtunity, to pass a id-Matrix to the
//                  general BDB-Integrator is not preferable due to overhead in
//                  computing an identity matrix of correct size in the coefficient
//                  function and an unnecessay matrix-matrix calculation
// 
//        Version:  1.0
//        Created:  10/29/2011 02:39:39 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#ifndef FILE_BBINT
#define FILE_BBINT

#include "integrator.hh"
#include "Elements/basefe.hh"
#include <boost/tr1/type_traits.hpp>
#include "Elements/HCurlElemsHi.hh"


namespace CoupledField {

  //! general class for calculation of bb forms
  template<template<class,class> class B_OP,
            class FE_TYPE,
            class MAT_DATA_TYPE=Double>
  class BBInt : public Integrator {
    public:

      //! Constructor with pointer to BaseElem
      BBInt(MAT_DATA_TYPE factor){
        this->name_ = "BBInt";
        isSymmetric_ = true;
        factor_ = factor;
      }

      //! Destructor
      ~BBInt(){

      }

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

  //! general class for calculation of bb forms
  template<template<class,class> class B_OP,
            class FE_TYPE,
            class MAT_DATA_TYPE=Double>
  class BBIntMassEdge : public BBInt<B_OP,FE_TYPE,MAT_DATA_TYPE> {
    public:

      //! Constructor with pointer to BaseElem
      BBIntMassEdge(MAT_DATA_TYPE factor):
        BBInt<B_OP,FE_TYPE,MAT_DATA_TYPE>(factor){
        this->name_ = "BBIntMassEdge";
      }

      //! Destructor
      ~BBIntMassEdge(){

      }

      //! Compute element matrix associated to BDB form
      void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                                 EntityIterator& ent1,
                                 EntityIterator& ent2 );
  };

  //=========================================================================================
  //Implementation of BB Int
  //=========================================================================================

  //! Calculate the
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
      //these two lines. This functionaliy has to be moved to GetFe in FeSpaceHCulHi!
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
#endif
