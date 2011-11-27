// =====================================================================================
// 
//       Filename:  BiLinearForm.hh
// 
//    Description:  Base Class for all bilinear integrators within CFS++
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

#ifndef FILE_BILINEARFORM_
#define FILE_BILINEARFORM_

#include "Domain/CoefFunction/CoefFunction.hh"
#include "FeBasis/FeSpace.hh"
#include "Domain/Domain.hh"

namespace CoupledField
{
  // forward class declaration
  class FeSpace;

  class BiLinearForm{
    public:
      BiLinearForm(){

      }

      ~BiLinearForm(){

      }



      virtual void CalcElementMatrix( Matrix<Double>& stiffMat,
                                          EntityIterator& ent1,
                                          EntityIterator& ent2){
        EXCEPTION("BiLinearForm::CalcElementMatrix called this may not happen!");
      }

      virtual void CalcElementMatrix( Matrix<Complex>& stiffMat,
                                          EntityIterator& ent1,
                                          EntityIterator& ent2){
        EXCEPTION("Integrator::CalcElementMatrix called this may not happen!");
      }

      //! Set finite element space in cases of mixed spaces
      //! To remain generality, each subclass needs this function too
      virtual void SetFeSpace( shared_ptr<FeSpace> feSpace1, shared_ptr<FeSpace> feSpace2)=0;

      //! Set Finite Element Space
      virtual void SetFeSpace( shared_ptr<FeSpace> feSpace )=0;

      std::string GetName(){
        return name_;
      }

      virtual bool IsComplex()=0;

    protected:

      //! name of (bi)linearform
      std::string name_;

      //! is the (bi) linear form symmetric
      bool isSymmetric_;

      //! pointer to finite element space 1
      shared_ptr<FeSpace> ptFeSpace1_;

      //! pointer to finite element space 2
      shared_ptr<FeSpace> ptFeSpace2_;

  };
}

#endif
