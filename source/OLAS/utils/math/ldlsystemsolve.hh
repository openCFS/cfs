// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef LDL_SYSTEM_SOLVE_HH
#define LDL_SYSTEM__STEP_HH

#include <vector>

#include "General/defs.hh"

// Add debug macro for derived class' macros here in order to get
// automatic debugging of the Solve step of the solver / preconditioner
// when turning on debugging of the solver / preconditioner
#ifdef DEBUG_LDLSOLVER_SOLVE
#define DEBUG_LDLSYSTEMSOLVE
#endif


namespace CoupledField {

  template <typename T> class Vector;
  
  //! Solution of a linear system of the form \f$LDL^T x = y\f$

  //! This class implements the solution of a linear system of the form
  //! \f$LDL^T x = y\f$ for the vector \f$x\f$ of unknowns. This functionality
  //! is used in the LDLSolver class and the corresponding preconditioner
  //! classes, which are derived from this one.
  template <typename T>
  class LDLSystemSolve {

  public:

    //! Default constructor
    LDLSystemSolve() {
    };

    //! Destructor
    virtual ~LDLSystemSolve() {
    };

    //! Solution of \f$LDL^Tx=y\f$ by forward backward substitutions

    //! \param cidx      array containing column indices of U
    //! \param rptr      array containing row pointers of U
    //! \param entry     array containing non-zero entries of U
    //! \param invDiag   array containing inverses of entries of diagonal
    //!                  matrix D
    //! \param x         vector containing the solution of the linear system
    //! \param y         vector containing right hand side of linear system
    //! \param dim       dimension of problem matrix
    void SolveLDLSystem( const UInt *cidx, const UInt *rptr, const T *entry,
                         const T *invDiag, Vector<T> &x, const Vector<T> &y,
                         UInt dim ) const {



      // We solve the linear system L D L^T x = y by performing the following
      // three steps:
      //
      // (1) L z = y       <--  forward substitution
      //
      // (2) w = D^{-1} z  <--  scaling
      //
      // (3) L^T x = w     <-- backward substitution
      //
      // The intermediate vectors are omitted in the implemenation, since
      // the steps can be done in place.


      UInt i, k;


      // ============================================
      //  Forward substitution with L = U^T (step 1)
      // ============================================

#ifdef DEBUG_LDLSYSTEMSOLVE
      (*debug) << "Starting forward substitution" << std::endl;
#endif

      // U is stored in CRS format, so L = U^T is "stored" in CCS format.
      // The diagonal entries are omitted, since they are all equal to one.
      // To solve L z = y we make use of the representation
      //
      // y = l1 * z1 + l2 * z2 + ... + ln  * zn
      //
      // where zk is the k-th component of z and lk the k-th column vector
      // of L. So after zk was computed we replace y by y - lk * zk.
      for ( k = 0; k < dim; k++ ) {
        x[k] = y[k];
      }

      for ( k = 0; k < dim; k++ ) {
        for ( i = rptr[k]; i < rptr[k+1]; i++ ) {
          x[ cidx[i] ] -= entry[i] * x[k];
        }
      }


      // ============================================
      //  Backward substitution with U (steps 2 + 3)
      // ============================================

#ifdef DEBUG_LDLSYSTEMSOLVE
      (*debug) << "Starting backward substitution (incl. scaling)"
               << std::endl;
#endif

      // U is stored in CRS format, where the entries of a row are in
      // ascending lexicographical ordering w.r.t. their column indices.
      // The diagonal entries are omitted, since they are all equal to one.
      k = dim;
      while ( k > 0 ) {
        k--;
        x[k] = invDiag[k] * x[k];
        for ( i = rptr[k]; i < rptr[k+1]; i++ ) {
          x[k] = x[k] - entry[i] * x[ cidx[i] ];
        }
      }
    }
  };

}

#endif
