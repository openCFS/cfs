// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

#include "multigrid/jacobi.hh"

/**********************************************************/
#ifdef DEBUG_TO_CERR
#ifndef DEBUG_JACOBI
#define DEBUG_JACOBI
#endif // DEBUG_JACOBI
#define  debug  &std::cerr
#endif // DEBUG_TO_CERR
/**********************************************************/

namespace OLAS {
/**********************************************************/

template <typename T>
Jacobi<T>::Jacobi()
    : DiagonalInverse_( NULL ),
      Size_( 0 ),
      Omega_( 0.8 ),
      PenaltyFlags_( NULL )
{
    ENTER_FCN("Jacobi::Jacobi");
}

/**********************************************************/

template <typename T>
Jacobi<T>::~Jacobi()
{
    ENTER_FCN("Jacobi::~Jacobi");
    
    Reset();
}

/**********************************************************/

template <typename T>
bool Jacobi<T>::Setup( const CRS_Matrix<T>& matrix )
{
    ENTER_FCN("Jacobi::Setup");

#ifdef  DEBUG_JACOBI
    if( matrix.GetNrows() <= 0 ) {
        Warning( "Jacobi::Setup called with an empty "
                 "matrix -> object reseted", __FILE__, __LINE__ );
        Reset();
        return false;
    }
    for( Integer i = 1; i < matrix.GetNrows(); i++ ) {
        if( matrix.GetColPointer()[matrix.GetRowPointer()[i]] != i ) {
            Error( "Jacobi::Setup: non-diagonal at leading position",
                   __FILE__, __LINE__ );
        }
        if( matrix.GetDataPointer()[matrix.GetRowPointer()[i]] == (T)0 ) {
            Error( "Jacobi::Setup: zero diagonal entry",
                   __FILE__, __LINE__ );
        }
    }
#endif

    // Create a new array for the diagonal inverses, only if the
    // old one cannot be reused. So first check, whether the old
    // one is present and has appropriate size.
    // Note: (Size_ == 0) <==> (DiagonalInverse_ == NULL)
    // So it is sufficient to check (Size_ != matrix.GetNrows()),
    // except of calls with an empty matrix.
    if( Size_ != matrix.GetNrows() ) {
        DeleteArray( DiagonalInverse_ );
        DiagonalInverse_ = NULL;
    }
    // create a new array for the diagonal inverses
    if( DiagonalInverse_ == NULL ) {
        Size_ = matrix.GetNrows();
        NewArray( DiagonalInverse_, T_Mtype, Size_ );
    }

    const Integer *const pRow = matrix.GetRowPointer();
    const T_Mtype *const pDat = matrix.GetDataPointer();
    // fill the array with the inverses of the diagonal entries
    for( Integer i = 1; i <= Size_; i++ ) {
        DiagonalInverse_[i] = opType<T_Mtype>::invert(pDat[pRow[i]]);
    }

    // create auxiliary vector
    auxVec_.Resize( matrix.GetNrows() );

    // set flag for the prepared state
    return this->prepared_ = true;
}

/**********************************************************/

template <typename T>
bool Jacobi<T>::Setup( const CRS_Matrix<T>& matrix,
                       const bool *const    penalty_flags )
{
    ENTER_FCN("Jacobi::Setup");


    if( Setup(matrix) ) {
        PenaltyFlags_ = penalty_flags;
        return true;
    } else {
        return false;
    }
}

/**********************************************************/

template <typename T>
void Jacobi<T>::
Step( const CRS_Matrix<T>&                  matrix,
      const Vector<T>&                      rhs,
            Vector<T>&                      sol,
      const typename Smoother<T>::Direction direction,
      const bool                            force_setup )
{
    ENTER_FCN("Jacobi::Step");

    // call Setup, if object is not prepared or
    // preparation forced by force_setup == true
    if( !(this->IsPrepared()) || force_setup ) {
        if( false == Setup(matrix) )  return;
    }

#ifdef  DEBUG_JACOBI
    // check, wheter the setup is suppressed explicitely, although
    // the matrix size and the size of the Jacobi preparation do
    // not match
    if( matrix.GetNrows() != Size_ ) {
        Warning( "Jacobi::Step: non-matching dimensions",
                 __FILE__, __LINE__ );
        return;
    }
#endif

    // apply one Jacobi step
    T_Vtype accumulator;
    const Integer *const RowPtr  = matrix.GetRowPointer();
    const Integer *const ColPtr  = matrix.GetColPointer();
    const T       *const DataPtr = matrix.GetDataPointer();

    //////////////////////////////////////////////
    // If there are penalty flags, use them.
    //////////////////////////////////////////////

    if( PenaltyFlags_ ) {
        for( Integer i = 1; i <= Size_; i++ ) {
            if( PenaltyFlags_[i] ) {
                auxVec_[i] = rhs[i] * DiagonalInverse_[i];
            } else {
                accumulator = rhs[i];
                // skip the diagonal entry at index RowPtr[i] and
                // start with the first offdiagonal entry at RowPtr[i]+1
                for( Integer ij = RowPtr[i]+1; ij < RowPtr[i+1]; ij++ ) {
                    accumulator -= DataPtr[ij] * sol[ColPtr[ij]];
                }
                // update sol[i]; I guess, it is OK to call this decition
                // (Omega_ == 1.0) Size_ times (and although this is a floating
                // point comparison
                if( Omega_ == 1.0 ) {
                    auxVec_[i] = DiagonalInverse_[i] * accumulator;
                } else {
                    auxVec_[i] =  (1.0 - Omega_) * sol[i]
                                 + Omega_ * DiagonalInverse_[i] * accumulator;
                }
            }
        }

    //////////////////////////////////////////////
    // normal Jacobi
    //////////////////////////////////////////////

    } else {
        for( Integer i = 1; i <= Size_; i++ ) {
            accumulator = rhs[i];
            // skip the diagonal entry at index RowPtr[i] and
            // start with the first offdiagonal entry at RowPtr[i]+1
            for( Integer ij = RowPtr[i]+1; ij < RowPtr[i+1]; ij++ ) {
                accumulator -= DataPtr[ij] * sol[ColPtr[ij]];
            }
            // update sol[i]; I guess, it is OK to call this decition
            // (Omega_ == 1.0) Size_ times (and although this is a floating
            // point comparison
            if( Omega_ == 1.0 ) {
                auxVec_[i] = DiagonalInverse_[i] * accumulator;
            } else {
                auxVec_[i] =  (1.0 - Omega_) * sol[i]
                             + Omega_ * DiagonalInverse_[i] * accumulator;
            }
        }
    }

    // copy new solution on passed vector
    for( Integer i = 1; i <= Size_; i++ )  sol[i] = auxVec_[i];
}

/**********************************************************/

template <typename T>
void Jacobi<T>::Reset()
{
    ENTER_FCN("Jacobi::Reset");
    
    DeleteArray( DiagonalInverse_ ); // delete diagonal inverse
    DiagonalInverse_ = NULL;
    auxVec_.Resize( 0 );
    Size_  =   0; // reset the size of the LES
    Omega_ = 1.0; // reset damping factor

    DeleteArray( PenaltyFlags_ );
    PenaltyFlags_ = NULL;

    // call Reset() of base class
    Smoother<T>::Reset();
}

/**********************************************************/
} // namespace OLAS

/**********************************************************/
#ifdef DEBUG_TO_CERR
#undef debug
#endif // DEBUG_TO_CERR
/**********************************************************/
