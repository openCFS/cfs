// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_SSORPRECOND_HH
#define OLAS_SSORPRECOND_HH

#include <def_expl_templ_inst.hh>

#include "baseprecond.hh"
#include "bnprecond.hh"

namespace CoupledField {

  template<typename> class Vector;
  template<typename> class CRS_Matrix;
  class StdMatrix;
  

  //! SSOR preconditioner

  //! This class implements iterative preconditioning by application of
  //! Symmetric Successive Over-Relaxation (SSOR). Applying the preconditioner
  //! as \f$M^{-1} r = z\f$ means that we perform one iteration sweep of SSOR
  //! to obtain an approximate solution for \f$Mz=r\f$ starting from a zero
  //! initial guess.
  template <typename T>
  class SSORPrecond : public BNPrecond<SSORPrecond<T>,CRS_Matrix<T>,T> {

  public:

    using BNPrecond<SSORPrecond<T>,CRS_Matrix<T>,T>::Apply;
    using BNPrecond<SSORPrecond<T>,CRS_Matrix<T>,T>::Setup;

    //! Typename of matrix entries (=T)
    typedef typename AssocType<T>::T_Mtype T_Mtype;

    //! Tiny vector of the same dimension as matrix block
    typedef typename AssocType<T>::T_Vtype T_Vtype;

    //! Scalar of the same primitive data type as matrix
    typedef typename AssocType<T>::T_Stype T_Stype;

    //! Matrix class associated with this preconditioner
    typedef CRS_Matrix<T> MyMatrixClass;

    //! Vector class associated with this preconditioner
    typedef Vector<T> MyVectorClass;

    // \param size: number of matrix rows/cols
    // No source code yet! Thus, this will counter instantiation.
    // SSORPrecond( Integer size );

    //! Constructor

    //! This constructor takes as input a system matrix from which the problem
    //! size, matrix and entry types are derived and two pointers to the
    //! communication objects. This is the constructor required by the
    //! GeneratePrecondObject function.
    SSORPrecond( const StdMatrix &mat, PtrParamNode myParams, 
                 PtrParamNode olasInfo);

    //! Default Destructor

    //! The default destructor is deep
    //! and frees internally allocated memory
    ~SSORPrecond();

    //! Apply SSOR preconditioner

    //! This method applies the SSOR preconditioner. It acts as a driver
    //! routine and calls the appropriate private method which implements
    //! SOR for the specific storage type of the matrix.
    void Apply( const CRS_Matrix<T> &sysmat,
		const Vector<T> &r, Vector<T> &z ) ;

    //! Triggers setup of the SSOR Preconditioner

    //! The setup phase generates a vector containing the inverses of the
    //! diagonal entries of the system matrix
    void Setup( CRS_Matrix<T> &sysMat, PtrParamNode analysis_id );

    //! Query type of preconditioner object

    //! When called this method returns the type of the preconditioner object.
    //! In the case of an object of this class the return value is JACOBI.
    BasePrecond::PrecondType GetPrecondType() const{
      return BasePrecond::SSOR;
    };


  private:

    //! Default constructor

    //! The default constructor is not allowed, since we need size information
    //! and pointers to communication objects for corrected initialisation.
    SSORPrecond(){
      EXCEPTION( "Default constructor of SSORPrecond should never be called!" );
    };
 
    //! Array containing inverses of diagonal entries of system matrix
    T_Mtype *diagInv_;

    //! Number of unknowns in linear system

    //! We store the number of unknowns in linear system since this is also
    //! the number of diagonal entries in the system matrix and thus the
    //! length of the internal data array storing their inverses.
    UInt size_;

  };

}//namespace

#ifndef EXPLICIT_TEMPLATE_INSTANTIATION
//#include "ssorprecond.cc"
#endif

#endif // OLAS_SSORPRECOND_HH
