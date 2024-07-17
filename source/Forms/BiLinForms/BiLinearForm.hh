// =====================================================================================
// 
//       Filename:  BiLinearForm.hh
// 
//    Description:  Base Class for all bilinear integrators within openCFS
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
      PRIM_PRIM,
      SEC_SEC,
      PRIM_SEC,
      SEC_PRIM
    } CouplingDirection;

      BiLinearForm( bool coordUpdate = false ){
        coordUpdate_ = coordUpdate;
        isSymmetric_ = false;
        isNewtonBiLinearForm_ = false;
        isSymmetric_ = false;
        isSolDependent_ = false;

        useVolEqnA_ = false;
        useVolEqnB_ = false;
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
        this->isNewtonBiLinearForm_ = right.isNewtonBiLinearForm_;
        this->isSymmetric_ = right.isSymmetric_;

        this->name_ = right.name_;

        // we just copy the feSpace pointers and need to make sure
        // not to alter their state...
        this->ptFeSpace1_ = right.ptFeSpace1_;
        this->ptFeSpace2_ = right.ptFeSpace2_;
        this->intScheme_ = right.intScheme_;
        this->isSolDependent_ = false;

        this->useVolEqnA_ = right.useVolEqnA_;
        this->useVolEqnB_ = right.useVolEqnB_;
      }

      /** Create a deep copy of the current objects pointer in combination
       *  with meaningful copy constructors
       */
      virtual BiLinearForm* Clone()=0;

      virtual ~BiLinearForm() {}

      typedef enum { NO_BILIN_TYPE = -1, BILIN_WRAPPED_LIN, BDB_INT, BB_INT, AB_INT, ADB_INT, SINGLE_ENTRY_BILIN_INT, IC_MODES_INT} Type;
      static Enum<Type> type;

      Type GetType() const { return type_; }

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
      void SetNewtonBiLinearForm() {
        isNewtonBiLinearForm_ = true;
      }

      //! Return if bilinearform is part of Newton tangential matrix
      bool IsNewtonBiLinearForm() {
        return isNewtonBiLinearForm_;
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

      IntScheme* GetIntScheme() const { return intScheme_.get(); }
      shared_ptr<IntScheme> GetPtrIntScheme() { return intScheme_; }

      FeSpace* GetFeSpace1() const { return ptFeSpace1_.get(); }
      shared_ptr<FeSpace> GetPtrFeSpace1() { return ptFeSpace1_; }

      //! Set eqn evaluation to volume for A operator
      virtual void SetUseVolEqnA( bool useVolEqn ) {
        useVolEqnA_ = useVolEqn;
      }

      //! Get eqn evaluation for A operator
      virtual bool GetUseVolEqnA() {
        return useVolEqnA_;
      }

      //! Set eqn evaluation to volume for B operator
      virtual void SetUseVolEqnB( bool useVolEqn ) {
        useVolEqnB_ = useVolEqn;
      }

      //! Get eqn evaluation for B operator
      virtual bool GetUseVolEqnB() {
        return useVolEqnB_;
      }

    protected:

      /** name of (bi)linearform. This is in the constructor BDBInt, ... and can be overwritten by SetName() */
      std::string name_;

      Type type_ = NO_BILIN_TYPE;

      //! is the (bi)linear form symmetric
      bool isSymmetric_;
      
      //! flag for use of updated Lagrangian formulation
      bool coordUpdate_; 

      //! is the (bi)linear part of the Newton tangential matrix
      bool isNewtonBiLinearForm_;

      //!depends on the solution
      bool isSolDependent_;

      //! pointer to finite element space 1
      shared_ptr<FeSpace> ptFeSpace1_;

      //! pointer to finite element space 2
      shared_ptr<FeSpace> ptFeSpace2_;

      //! point to integration scheme
      shared_ptr<IntScheme> intScheme_;

      // Flag to indicate if the number of functions shall be aquired from the volume or surface element
      // This is needed for e.g. a gradient evaluated at a surface, since we perform the evaluation at the
      // integration points of the surface, but we need all DOFs of the element for the calculation of the gradient
      // We can set this for A- and B operator seperately
      bool useVolEqnA_;

      bool useVolEqnB_;
  };
}

#endif
