// =====================================================================================
// 
//       Filename:  BiLinearForm.hh
// 
//    Description:  Base Class for all bilinear integrators within CFS++
//                  Most important derived classes so far:
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
#include "Utils/ThreadLocalStorage.hh"

namespace CoupledField
{
  // forward class declaration
  class FeSpace;

class BiLinearForm : public CfsCopyable{
    public:

    //for NMG integrators
    typedef enum {
      MASTER_MASTER,
      SLAVE_SLAVE,
      MASTER_SLAVE,
      SLAVE_MASTER
    } CouplingDirection;

      BiLinearForm( bool coordUpdate = false ){
        coordUpdate_ = coordUpdate;
        isSymmetric_ = false;
        isNewtonBilinearForm_ = false;
        isSymmetric_ = false;
        isSolDependent_ = false;
      }

      /** This assignment operator is only! designed for use for OMP
       *  For other purposes, its applicability is highly questionable.
       *  Furthermore, the default copy constructor is assumed to
       *  work just fine...
       *  Funny thing, the usage of this constructor is not threadsafe
       *  but object creation and access needs to be synchronized anyway
       *
       *  In general: operators are assumed to be lightweight, so
       *  we can affort a copy
       *  CoefFunctions may be not, so we need to make them thread safe...
       */
      BiLinearForm(const BiLinearForm& right){
        this->coordUpdate_ = right.coordUpdate_;
        this->isNewtonBilinearForm_ = right.isNewtonBilinearForm_;
        this->isSymmetric_ = right.isSymmetric_;

        this->name_ = right.name_;

        // we just copy the feSpace pointers and need to make sure
        // not to alter their state...
        this->ptFeSpace1_ = right.ptFeSpace1_;
        this->ptFeSpace2_ = right.ptFeSpace2_;
        this->intScheme_ = right.intScheme_;
        this->isSolDependent_ = false;
      }

      /** Create a deep copy of the current objects pointer in combination
       *  with meaningful copy constructors
       */
      virtual BiLinearForm* Clone()=0;


      virtual ~BiLinearForm(){

      }

      virtual void CalcElementMatrix( Matrix<Double>& stiffMat,
                                          EntityIterator& ent1,
                                          EntityIterator& ent2){
        EXCEPTION("BiLinearForm::CalcElementMatrix called, this must not happen!");
      }

      virtual void CalcElementMatrix( Matrix<Complex>& stiffMat,
                                          EntityIterator& ent1,
                                          EntityIterator& ent2){
        EXCEPTION("BiLinearForm::CalcElementMatrix called, this must not happen!");
      }

      //! Set finite element space in cases of mixed spaces
      //! To remain generality, each subclass needs this function too
      virtual void SetFeSpace( shared_ptr<FeSpace> feSpace1, shared_ptr<FeSpace> feSpace2)=0;

      //! Set Finite Element Space
      virtual void SetFeSpace( shared_ptr<FeSpace> feSpace )=0;

      //! Return name of bilinear form
      const std::string& GetName() const {
        return name_;
      }

      //! Set name of bilinear form
      void SetName(const std::string& name ){
        name_ = name;
      }

      virtual bool IsComplex() const =0;
      
      //! Return, if bilinear form produces symmetric matrices
      bool IsSymmetric() {return isSymmetric_;}
      
      //! Return if element matrix is solution dependent
      virtual bool IsSolDependent() = 0;
      
      //! Return if bilinearform uses updated Lagrangian formulation
      bool IsCoordUpdate() { return coordUpdate_; }
      
      /** set coordUpdate, required by shape opt */
      void SetCoordUpdate(bool val) { coordUpdate_ = val; }

      //! set bilinearform to part of Newton tangential matrix
      void SetNewtonBilinearForm() {
        isNewtonBilinearForm_ = true;
      }

      //! Return if bilinearform is part of Newton tangential matrix
      bool IsNewtonBilinearForm() {
        return isNewtonBilinearForm_;
      }

      //! Set Coefficient Function of B operator
      virtual void SetBCoefFunctionOpB(PtrCoefFct coef){
        EXCEPTION("Integrator::SetBCoefFunctionOpB not available in base class!");
      }

      //! Set Coefficient Function of A operator
      virtual void SetBCoefFunctionOpA(PtrCoefFct coef){
        EXCEPTION("Integrator::SetCoefFunctionOpA not available in base class!");
      }

      //! \copydoc BiLinearForm::IsSolDependent
      virtual void SetSolDependent(bool depend) {
    	  isSolDependent_ = depend;
      }

      IntScheme* GetIntScheme() { return intScheme_.get(); }

      FeSpace* GetFeSpace1() { return ptFeSpace1_.get(); }

    protected:

      //! name of (bi)linearform
      std::string name_;

      //! is the (bi)linear form symmetric
      bool isSymmetric_;
      
      //! flag for use of updated Lagrangian formulation
      bool coordUpdate_; 

      //! is the (bi)linear part of the Newton tangential matrix
      bool isNewtonBilinearForm_;

      //!depends on the solution
      bool isSolDependent_;

      //! pointer to finite element space 1
      shared_ptr<FeSpace> ptFeSpace1_;

      //! pointer to finite element space 2
      shared_ptr<FeSpace> ptFeSpace2_;

      //! point to integration scheme
      shared_ptr<IntScheme> intScheme_;
  };
}

#endif
