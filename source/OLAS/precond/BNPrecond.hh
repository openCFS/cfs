#ifndef OLAS_BNPRECOND_HH
#define OLAS_BNPRECOND_HH

#include "MatVec/Vector.hh"
#include "BasePrecond.hh"


namespace CoupledField {

  //! This class is the Barton-Nackman base class for
  //! standard preconditioners. It essentially casts
  //! standard Matrices and Vectors into template types
  template<class T_precond, class T_storage, typename T>
  class BNPrecond : public BasePrecond {

  public:
    
    virtual ~BNPrecond() {}

    using BasePrecond::Apply;
    using BasePrecond::Setup;

    typedef typename AssocType<T>::T_Mtype T_Mtype;
    typedef typename AssocType<T>::T_Vtype T_Vtype;
    typedef typename AssocType<T>::T_Stype T_Stype;

    typedef T_storage AssocMatrixType;
    typedef Vector<T> AssocVectorType;

    //! casts the preconditioner into its actual class
    T_precond& SubClass(){return static_cast<T_precond&>(*this);}

    //! casts the preconditioner into its actual class (const version)
    const T_precond& SubClass() const {
      return static_cast<const T_precond&>(*this);
    }

    //! cast the matrix into its full type and perform the Setup 
    //! of the preconditioner as subclass
    virtual void Setup( BaseMatrix &A, PtrParamNode analysis_id ) {
      AssocMatrixType& idA = dynamic_cast<AssocMatrixType&>(A);
      SubClass().Setup(idA, analysis_id);
      readyToUse_ = true;
    }

    //! cast the Matrix and Vectors into their full types and
    //! apply the preconditioner if Setup has been called in advance
    virtual void Apply( const BaseMatrix &A, 
                        const BaseVector &r, BaseVector &z)  {
      const AssocMatrixType& idA = dynamic_cast<const AssocMatrixType&>(A);
      const AssocVectorType& idr = dynamic_cast<const AssocVectorType&>(r);
      AssocVectorType& idz = dynamic_cast<AssocVectorType&>(z);

      if (readyToUse_) {
        SubClass().Apply(idA,idr,idz);
      }
      else {
        EXCEPTION( "Preconditioner used before Setup!");
      }
    }
  };
}

#endif
