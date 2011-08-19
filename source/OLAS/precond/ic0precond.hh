// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_IC0PRECOND_HH
#define OLAS_IC0PRECOND_HH

#include <def_expl_templ_inst.hh>

#include "baseprecond.hh"
#include "bnprecond.hh"

namespace CoupledField {

  template<typename> class SCRS_Matrix;

  //! IC0 preconditioner

  //! This class implements iterative preconditioning by application of
  //! Incomplete Cholesky with zero fill in
  template <typename T>
  class IC0Precond : public BNPrecond<IC0Precond<T>,SCRS_Matrix<T>,T> {

  public:

    using BNPrecond<IC0Precond<T>,SCRS_Matrix<T>,T>::Apply;
    using BNPrecond<IC0Precond<T>,SCRS_Matrix<T>,T>::Setup;

    //! Constructor

    //! This constructor takes as input a system matrix from which the problem
    //! size, matrix and entry types are derived and two pointers to the
    //! communication objects. This is the constructor required by the
    //! GeneratePrecondObject function.
    IC0Precond( const StdMatrix &mat, PtrParamNode solverNode, PtrParamNode olasInfo);

    //! Default Destructor

    //! The default destructor is deep
    //! and frees internally allocated memory
    ~IC0Precond();

    //! Apply IC0 preconditioner

    //! This method applies the IC0 preconditioner and does
    //! the forward and backward substitution
    void Apply( const SCRS_Matrix<T> &sysmat,
		const Vector<T> &r, Vector<T> &z ) const;

    //! Triggers setup of the IC0 Preconditioner

    //! The setup phase generates the incomplete cholesky factorization
    //! with zero fill in
    void Setup( SCRS_Matrix<T> &sysMat );


    //! Query type of preconditioner object

    //! When called this method returns the type of the preconditioner object.
    //! In the case of an object of this class the return value is JACOBI.
    BasePrecond::PrecondType GetPrecondType() const{
      return BasePrecond::IC0;
    };


  private:

    //! Default constructor

    //! The default constructor is not allowed, since we need size information
    //! and pointers to communication objects for corrected initialisation.
    IC0Precond(){
      EXCEPTION( "Default constructor of IC0Precond should never be called!" );
    };

    //! Checks if a double value is larger a prescribed value
    UInt CheckPositivity(Double val);

     //! Checks if a complex value is larger a prescribed value
    UInt CheckPositivity(Complex val);
 
    //! Number of unknowns in linear system

    //! We store the number of unknowns in linear system since this is also
    //! the number of diagonal entries in the system matrix and thus the
    //! length of the internal data array storing their inverses.
    UInt size_;

    //@{ \name Administration of matrix factors
    //! The matrix factor \f$U=L^T\f$ is stored in CRS storage format with
    //! the diagonal entries, which we know are all one. Within each row, the
    //! non-zero entries are stored in lexicographic ordering with respect to
    //! their column indices starting from the diagonal.

    //! Values of non-zero entries of factor \f$U=L^T\f$
    T * dataU_;

    //! Column indices of non-zero entries factor \f$U=L^T\f$
    UInt * cidxU_;

    //! Pointers to row starts in contiguous storage of factor \f$U=L^T\f$
    UInt * rptrU_;

    //! The Setup routine will set this boolean to true, once it has computed
    //! a factorisation. This gives us a confortable way of keeping track of
    //! the status of the solver and the internal memory allocation.
    bool amFactorised_;
  };

}//namespace

#ifndef EXPLICIT_TEMPLATE_INSTANTIATION
//#include "ic0precond.cc"
#endif

#endif // OLAS_IC0PRECOND_HH
