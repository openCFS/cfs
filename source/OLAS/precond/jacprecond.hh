// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_JACPRECOND_HH
#define OLAS_JACPRECOND_HH

#include <string>
#include <iostream>
#include <fstream>


#include "baseprecond.hh"
#include "bnprecond.hh"
#include "MatVec/matrix.hh"

namespace CoupledField {

  // ========================================================================
  //  S C A L A R   J A C O B I   P R E C O N D I T I O N E R
  // ========================================================================

  //! Jacobi Preconditioner

  //! This class implements Jacobi preconditioning. Applying the Jacobi
  //! preconditioner gives us \f$z = M^{-1} r\f$ where
  //! \f$M=\mbox{diag}(a_{11},a_{22}, \ldots, a_{nn})\f$ is a diagonal
  //! matrix consisting of the diagonal entries of the problem matrix \f$A\f$.
  //! Thus, it is evident, that Jacobi preconditioning requires \f$A\f$ to
  //! not have any zero entries on the diagonal.
  template <class T_storage,typename T>
  class JacPrecond : public BNPrecond<JacPrecond<T_storage,T>,T_storage,T> {

  public:

    using BNPrecond<JacPrecond<T_storage,T>,T_storage,T>::Apply;
    using BNPrecond<JacPrecond<T_storage,T>,T_storage,T>::Setup;

    ///
    typedef typename AssocType<T>::T_Mtype T_Mtype;
    typedef typename AssocType<T>::T_Vtype T_Vtype;
    typedef typename AssocType<T>::T_Stype T_Stype;

    //! Constructor

    //! This constructor takes as input a system matrix from which the problem
    //! size, matrix and entry types are derived and two pointers to the
    //! communication objects. This is the constructor required by the
    //! GeneratePrecondObject function.
    JacPrecond( const StdMatrix &mat, PtrParamNode solverNode, 
                PtrParamNode olasInfo );

    //! Default Destructor
    ~JacPrecond();

    //! Scales the residual with the inverse diagonal of system matrix
    void Apply( const T_storage &sysmat, const Vector<T> &r,
		Vector<T> &z );

    //! Triggers setup of the Jacobi Preconditioner

    //! The setup phase generates a vector containing the inverses of the
    //! diagonal entries of the system matrix
    void Setup( T_storage &sysmat,
                PtrParamNode analysis_id);

    //! Query type of preconditioner object

    //! When called this method returns the type of the preconditioner object.
    //! In the case of an object of this class the return value is JACOBI.
    BasePrecond::PrecondType GetPrecondType() const {
      return BasePrecond::JACOBI;
    };

  private:
 
    //! Default constructor

    //! The default constructor is not allowed, since we need size information
    //! and pointers to communication objects for corrected initialisation.
    JacPrecond(){
      EXCEPTION( "Default constructor of JacPrecond should never be called!");
    };

    //! Array containing inverses of diagonal entries of system matrix
    T_Mtype *diagInv_;

    //! Dimension of system matrix and thus length of diagInv_
    UInt size_;

  };

  // ================================================================s========
  //   B L O C K   J A C O B I   P R E C O N D I T I O N E R
  // ========================================================================
  
  //! Block Jacobi Preconditioner

  //! This class implements Block-Jacobi preconditioning. Applying the Jacobi
  //! preconditioner gives us \f$z = M^{-1} r\f$ where
  //! \f$M=\mbox{diag}(a_{11},a_{22}, \ldots, a_{nn})\f$ is a diagonal
  //! matrix consisting of the diagonal entries of the problem matrix \f$A\f$.
  //! Thus, it is evident, that Jacobi preconditioning requires \f$A\f$ to
  //! not have any zero entries on the diagonal.
  template <class T_storage,typename T>
  class BlockJacPrecond : public BNPrecond<BlockJacPrecond<T_storage,T>,T_storage,T> {

  public:

    using BNPrecond<BlockJacPrecond<T_storage,T>,T_storage,T>::Apply;
    using BNPrecond<BlockJacPrecond<T_storage,T>,T_storage,T>::Setup;

    ///
    typedef typename AssocType<T>::T_Mtype T_Mtype;
    typedef typename AssocType<T>::T_Vtype T_Vtype;
    typedef typename AssocType<T>::T_Stype T_Stype;

    //! Constructor

    //! This constructor takes as input a system matrix from which the problem
    //! size, matrix and entry types are derived and two pointers to the
    //! communication objects. This is the constructor required by the
    //! GeneratePrecondObject function.
    BlockJacPrecond( const StdMatrix &mat, PtrParamNode solverNode, 
                     PtrParamNode olasInfo );

    //! Default Destructor
    ~BlockJacPrecond();

    //! Scales the residual with the inverse diagonal of system matrix
    void Apply( const T_storage &sysmat, const Vector<T> &r,
                Vector<T> &z );

    //! Triggers setup of the Jacobi Preconditioner

    //! The setup phase generates a vector containing the inverses of the
    //! diagonal entries of the system matrix
    void Setup( T_storage &sysmat,
                PtrParamNode analysis_id);

    //! Query type of preconditioner object

    //! When called this method returns the type of the preconditioner object.
    //! In the case of an object of this class the return value is BLOCK_JACOBI.
    BasePrecond::PrecondType GetPrecondType() const {
      return BasePrecond::BLOCK_JACOBI;
    };

  private:

    //! Default constructor

    //! The default constructor is not allowed, since we need size information
    //! and pointers to communication objects for corrected initialisation.
    BlockJacPrecond(){
      EXCEPTION( "Default constructor of JacPrecond should never be called!");
    };
    
    //! Number of block rows 
    UInt numRows_;
    
    //! Vector containing matrices with factorization
    StdVector<Matrix<T_Stype> > factors_;

  };

  
}//namespace

#endif // OLAS_JACPRECOND_HH
