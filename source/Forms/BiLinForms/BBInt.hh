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

#include "BiLinearForm.hh"
#include "FeBasis/BaseFE.hh"
#include <boost/tr1/type_traits.hpp>
#include "FeBasis/HCurl/HCurlElemsHi.hh"


namespace CoupledField {

  //! general class for calculation of bb forms
  template<template<class,class> class B_OP,
            class FE_TYPE,
            class MAT_DATA_TYPE=Double>
  class BBInt : public BiLinearForm {
    public:

      //! Constructor with pointer to BaseElem
      BBInt( shared_ptr<CoefFunction> scalCoef, MAT_DATA_TYPE factor);

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

      void SetFeSpace( shared_ptr<FeSpace> feSpace ) {
        this->ptFeSpace1_ = feSpace;
        UInt opDim = feSpace->GetFeFunction()->GetResultInfo()->dofNames.GetSize();
        Bdim_ = opDim;
      }

      virtual void SetFeSpace( shared_ptr<FeSpace> feSpace1, shared_ptr<FeSpace> feSpace2) {
        this->ptFeSpace1_ = feSpace1;
        this->ptFeSpace2_ = feSpace2;
        UInt opADim = feSpace1->GetFeFunction()->GetResultInfo()->dofNames.GetSize();
        UInt opBDim = feSpace2->GetFeFunction()->GetResultInfo()->dofNames.GetSize();
        Adim_ = opADim;
        Bdim_ = opBDim;
      }

    protected:
      B_OP<FE_TYPE,MAT_DATA_TYPE> operator_;
      
      //! Pointer to scalar (!!) coefficient for BB-integrator 
      shared_ptr<CoefFunction > scalCoef_;

      //! set a constant factor for multiplication with the element matrix
      MAT_DATA_TYPE factor_;

      //! dimension of a-operator (first B-Operator)
      UInt Adim_;

      //! dimension of b-operator (second B-Operator)
      UInt Bdim_;
  };

  //! general class for calculation of bb forms
  template<template<class,class> class B_OP,
            class FE_TYPE,
            class MAT_DATA_TYPE=Double>
  class BBIntMassEdge : public BBInt<B_OP,FE_TYPE,MAT_DATA_TYPE> {
    public:

      //! Constructor with pointer to BaseElem
      BBIntMassEdge(shared_ptr<CoefFunction> scalCoef, MAT_DATA_TYPE factor):
        BBInt<B_OP,FE_TYPE,MAT_DATA_TYPE>(scalCoef, factor){
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
}

//Include template definition file
#include "BBInt.cc"
#endif
