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

#include "biLinearForm.hh"
#include "Elements/basefe.hh"
#include <boost/tr1/type_traits.hpp>

namespace CoupledField {



  //! general class for calculation of bdb forms
  template<template<class,class> class B_OP,
            class FE_TYPE,
            class MAT_DATA_TYPE=Double,
            class COEF_DATA_TYPE=Double>
  class BDBInt : public BiLinearForm {
    public:

      //! Constructor with pointer to BaseElem
      BDBInt(shared_ptr<CoefFunction> dData, MAT_DATA_TYPE factor);

      //! Destructor
      virtual ~BDBInt();

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


      //! set a constant factor for multiplication with the element matrix
      MAT_DATA_TYPE factor_;

      //! Pointer to coefficient function computing the d-matrix of the BDB Integrator
      shared_ptr<CoefFunction > dData_;

      //! dimension of a-operator (first B-Operator)
      UInt Adim_;

      //! dimension of b-operator (second B-Operator)
      UInt Bdim_;
  };

}

//Include template definition file
#include "bdbInt.cc"

#endif
