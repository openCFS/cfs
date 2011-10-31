// =====================================================================================
// 
//       Filename:  integrator.hh
// 
//    Description:  Base Class for all available integrators within CFS++
//                  Most oimportant derived classes so far:
//                  BDBInt -> for symmetric bilinear forms
//                  ADBInt -> for asymmetric bilinear forms
// 
//        Version:  1.0
//        Created:  10/04/2011 09:31:15 AM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#ifndef FILE_INTEGRATOR_
#define FILE_INTEGRATOR_

#include "coefFunction.hh"
#include "Elements/fespace.hh"
#include "Domain/domain.hh"

namespace CoupledField
{
  // forward class declaration
  class FeSpace;

  class Integrator{
    public:
      Integrator(){

      }

      ~Integrator(){

      }
      virtual void CalcElemVector(Vector<Complex> elemVec,EntityIterator& ent){
        EXCEPTION("Integrator::CalcElementVector called this may not happen!");
      }

      virtual void CalcElemVector(Vector<Double> elemVec,EntityIterator& ent){
        EXCEPTION("Integrator::CalcElementVector called this may not happen!");
      }


      virtual void CalcElementMatrix( Matrix<Double>& stiffMat,
                                          EntityIterator& ent1,
                                          EntityIterator& ent2){
        EXCEPTION("Integrator::CalcElementMatrix called this may not happen!");
      }

      virtual void CalcElementMatrix( Matrix<Complex>& stiffMat,
                                          EntityIterator& ent1,
                                          EntityIterator& ent2){
        EXCEPTION("Integrator::CalcElementMatrix called this may not happen!");
      }

      void SetFeSpace( shared_ptr<FeSpace> feSpace1, shared_ptr<FeSpace> feSpace2) {
        this->ptFeSpace1_ = feSpace1;
        this->ptFeSpace2_ = feSpace2;
        UInt opADim = feSpace1->GetFeFunction()->GetResultInfo()->dofNames.GetSize();
        UInt opBDim = feSpace2->GetFeFunction()->GetResultInfo()->dofNames.GetSize();
        Adim_ = opADim;
        Bdim_ = opBDim;
      }

      //! Set Finite Element Space
      void SetFeSpace( shared_ptr<FeSpace> feSpace ) {
        this->ptFeSpace1_ = feSpace;
        UInt opDim = feSpace->GetFeFunction()->GetResultInfo()->dofNames.GetSize();
        Bdim_ = opDim;
      }

      std::string GetName(){
        return name_;
      }

      virtual bool IsComplex()=0;

    protected:

      //! name of (bi)linearform
      std::string name_;

      //! is the (bi) linear form symmetric
      bool isSymmetric_;

      //! Pointer to coefficient function computing the d-matrix of the BDB Integrator
      shared_ptr<CoefFunction > dData_;

      //! pointer to finite element space 1
      shared_ptr<FeSpace> ptFeSpace1_;

      //! pointer to finite element space 2
      shared_ptr<FeSpace> ptFeSpace2_;

      //! dimension of a-operator
      UInt Adim_;

      //! dimension of b-operator
      UInt Bdim_;

  };
}

#endif
