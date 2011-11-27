// =====================================================================================
// 
//       Filename:  bdbInt.hh
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

#include <boost/tr1/type_traits.hpp>

#include "BiLinearForm.hh"

#include "FeBasis/BaseFE.hh"


namespace CoupledField {


  //! Base class for BDB integrator

  //! This is the base class for all BDB-integrators. It mainly serves as
  //! a common hook to store BDB-integrators, independent of the type of
  //! finite element functions and the coefficient data, by which the derived
  //! classes are templatized.
  //! This interface mainly provides abstract methods to apply the B-operator
  //! and the dB-operator to vectors, which is used e.g. in postprocessing to
  //! evaluate flux and energy quantities.
class BaseBDBInt : public BiLinearForm {
public:

  //@{
  //! Apply / multiply element matrix to vector

  //! This method multiplies the element matrix with a vector. This can be used e.g.
  //! to calculate the energy of a field on the element level. 
  virtual void ApplyElemMat( Vector<Double>&ret, const Vector<Double>& sol,
                             EntityIterator& ent1,
                             EntityIterator& ent2 ) {
    EXCEPTION("BaseBDBInt: ApplyBMat not implemented");
  }
  virtual void ApplyElemMat( Vector<Complex>&ret, const Vector<Double>& sol,
                             EntityIterator& ent1,
                             EntityIterator& ent2 ) {
    EXCEPTION("BaseBDBInt: ApplyBMat not implemented");
  }  
  //@}

  //@{
  //! Apply B-Operator on vector
  virtual void ApplyBMat( Vector<Double>&ret, const Vector<Double>& sol,
                          const LocPointMapped& lpm ) {
    EXCEPTION("BaseBDBInt: ApplyBMat not implemented");
  }
  virtual void ApplyBMat( Vector<Complex>&ret, const Vector<Complex>& sol,
                          const LocPointMapped& lpm ) {
    EXCEPTION("BaseBDBInt: ApplyBMat not implemented");
  }
  //@}


  //@{
  //! Apply dB-Operator on vector
  virtual void ApplydBMat( Vector<Double>&ret, const Vector<Double>& sol,
                           const LocPointMapped& lpm ) {
    EXCEPTION("BaseBDBInt: ApplydBMat not implemented");
  }

  virtual void ApplydBMat( Vector<Complex>&ret, const Vector<Complex>& sol,
                           const LocPointMapped& lpm ) {
    EXCEPTION("BaseBDBInt: ApplydBMat not implemented");
  }
  //@}

  //@{
  //! Calculate integration kernel, i.e. B*d*B without integration
  virtual void CalcKernel( Matrix<Double>& kernel,
                           const LocPointMapped& lpm ) {
    EXCEPTION("BaseBDBInt: CalcKernel not implemented");
  }

  virtual void CalcKernel( Matrix<Complex>& kernel,
                           const LocPointMapped& lpm ) {
    EXCEPTION("BaseBDBInt: CalcKernel not implemented");
  }
  //@}

};


  //! general class for calculation of bdb forms
  template<template<class,class> class B_OP,
            class FE_TYPE,
            class MAT_DATA_TYPE=Double,
            class COEF_DATA_TYPE=Double>
  class BDBInt : public BaseBDBInt {
    public:

      //! Constructor with pointer to BaseElem
      BDBInt(shared_ptr<CoefFunction> dData, MAT_DATA_TYPE factor);

      //! Destructor
      virtual ~BDBInt();

      //! Compute element matrix associated to BDB form
      void CalcElementMatrix( Matrix<MAT_DATA_TYPE>& elemMat,
                              EntityIterator& ent1,
                              EntityIterator& ent2 );

      //! Multiply element matrix with vector
      template<class VEC_DATA_TYPE> 
      void ApplyElemMat( Vector<VEC_DATA_TYPE>&ret, const Vector<Double>& sol,
                         EntityIterator& ent1,
                         EntityIterator& ent2 );

      
      //! Calculate integration kernel, i.e. B*d*B without integration
      void CalcKernel( Matrix<MAT_DATA_TYPE>& kernel, 
                       const LocPointMapped& lpm );
      
      //! Apply B-Operator on vector
      //template<class VEC_DATA_TYPE>
      void ApplyBMat( Vector<MAT_DATA_TYPE>&ret, 
                      const Vector<MAT_DATA_TYPE>& sol,
                      const LocPointMapped& lpm );

      //! Apply dB-Operator on vector
      //template<class VEC_DATA_TYPE>
      void ApplydBMat( Vector<MAT_DATA_TYPE>&ret, 
                       const Vector<MAT_DATA_TYPE>& sol,
                       const LocPointMapped& lpm );

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
#include "BDBInt.cc"

#endif
