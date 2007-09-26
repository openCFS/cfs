// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_BNPRECOND_HH
#define OLAS_BNPRECOND_HH

namespace OLAS {

  //! This class is the Barton-Nackman base class for
  //! standard preconditioners. It essentially casts
  //! standard Matrices and Vectors into template types
  template<class T_precond, class T_storage, typename T>
  class BNPrecond : public BaseStdPrecond {

  public:

    typedef typename assocType<T>::T_Mtype T_Mtype;
    typedef typename assocType<T>::T_Vtype T_Vtype;
    typedef typename assocType<T>::T_Stype T_Stype;

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
    virtual void Setup( StdMatrix &A ) {
      TRY_CAST {
	RefCast( A, AssocMatrixType, idA );
      	SubClass().Setup(idA);
	readyToUse_ = true;
      } CATCH_CAST;
    }

    //! cast the Matrix and Vectors into their full types and
    //! apply the preconditioner if Setup has been called in advance
    virtual void Apply( const StdMatrix &A, 
			const SparseVector &r, SparseVector &z) const {
    	
      TRY_CAST {
	ConstRefCast(A,AssocMatrixType,idA);
	ConstRefCast(r,AssocVectorType,idr);
	RefCast(z,AssocVectorType,idz);

	if (readyToUse_) {
	  SubClass().Apply(idA,idr,idz);
	}
	else {
	  Error( "Preconditioner used before Setup!", __FILE__, __LINE__ );
	}
      } CATCH_CAST;
    }

  };

}

#endif
