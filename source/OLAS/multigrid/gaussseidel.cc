// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

#include "multigrid/gaussseidel.hh"

/**********************************************************/
#ifdef DEBUG_TO_CERR
#ifndef DEBUG_GAUSSSEIDEL
#define DEBUG_GAUSSSEIDEL
#endif // DEBUG_GAUSSSEIDEL
#define  debug  &std::cerr
#endif // DEBUG_TO_CERR
/**********************************************************/

namespace OLAS {
/**********************************************************/

template <typename T>
GaussSeidel<T>::GaussSeidel()
    : DiagonalInverse_( NULL ),
      Size_( 0 ),
      Omega_( 1.0 ),
      PenaltyFlags_( NULL )
{
}

/**********************************************************/

template <typename T>
GaussSeidel<T>::~GaussSeidel()
{
    
    Reset();
}

/**********************************************************/

template <typename T>
bool GaussSeidel<T>::Setup( const CRS_Matrix<T>& matrix )
{

#ifdef  DEBUG_GAUSSSEIDEL
    if( matrix.GetNrows() <= 0 ) {
        Warning( "GaussSeidel::Setup called with an empty "
                 "matrix -> object reseted", __FILE__, __LINE__ );
        Reset();
        return false;
    }
    for( Integer i = 1; i < matrix.GetNrows(); i++ ) {
        if( matrix.GetColPointer()[matrix.GetRowPointer()[i]] != i ) {
            Error( "GaussSeidel::Setup: non-diagonal at leading position",
                   __FILE__, __LINE__ );
        }
        if( matrix.GetDataPointer()[matrix.GetRowPointer()[i]] == (T)0 ) {
            Error( "GaussSeidel::Setup: zero diagonal entry",
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

    // set flag for the prepared state
    return this->prepared_ = true;
}

/**********************************************************/

template <typename T>
bool GaussSeidel<T>::Setup( const CRS_Matrix<T>& matrix,
                            const bool *const    penalty_flags )
{


    if( Setup(matrix) ) {
        PenaltyFlags_ = penalty_flags;
        return true;
    } else {
        return false;
    }
}

/**********************************************************/

template <typename T>
void GaussSeidel<T>::
Step( const CRS_Matrix<T>&                  matrix,
      const Vector<T>&                      rhs,
            Vector<T>&                      sol,
      const typename Smoother<T>::Direction direction,
      const bool                            force_setup )
{

    // call Setup, if object is not prepared or
    // preparation forced by force_setup == true
    if( !(this->IsPrepared()) || force_setup ) {
        if( false == Setup(matrix) )  return;
    }

#ifdef  DEBUG_GAUSSSEIDEL
    // check, wheter the setup is suppressed explicitely, although
    // the matrix size and the size of the GaussSeidel preparation
    // do not match
    if( matrix.GetNrows() != Size_ ) {
        Warning( "suppressed setup in GaussSeidel::Step, but the "
                 "dimensions do not match", __FILE__, __LINE__ );
        return;
    }
#endif

    // apply one Gauss-Seidel step
    T_Vtype accumulator;
    const Integer *const RowPtr   = matrix.GetRowPointer();
    const Integer *const ColPtr   = matrix.GetColPointer();
    const T       *const DataPtr  = matrix.GetDataPointer();

    //////////////////////////////////////////////
    // If there are penalty flags, use them.
    //////////////////////////////////////////////
    if( PenaltyFlags_ ) {
        if( direction == Smoother<T>::FORWARD ) {
            // Gauss Seidel FORWARD step
            for( Integer i = 1; i <= Size_; i++ ) {
                if( PenaltyFlags_[i] ) {
                    sol[i] = rhs[i] * DiagonalInverse_[i];
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
                        sol[i] = DiagonalInverse_[i] * accumulator;
                    } else {
                        sol[i] += Omega_ *
                                  ( (DiagonalInverse_[i] * accumulator)
                                   - sol[i]);
                    }
                }
            }
        } else {
            // Gauss Seidel BACKWARD step
            for( Integer i = Size_; i > 0; i-- ) {
                if( PenaltyFlags_[i] ) {
                    sol[i] = rhs[i] * DiagonalInverse_[i];
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
                        sol[i] = DiagonalInverse_[i] * accumulator;
                    } else {
                        sol[i] += Omega_ *
                                  ( (DiagonalInverse_[i] * accumulator)
                                   - sol[i]);
                    }
                }
            }
        }
    //////////////////////////////////////////////
    // normal Gauss-Seidel
    //////////////////////////////////////////////
    } else {
        if( direction == Smoother<T>::FORWARD ) {
            // Gauss Seidel FORWARD step
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
                    sol[i] = DiagonalInverse_[i] * accumulator;
                } else {
                    sol[i] += Omega_ *
                              ((DiagonalInverse_[i] * accumulator) - sol[i]);
                }
            }
        } else {
            // Gauss Seidel BACKWARD step
            for( Integer i = Size_; i > 0; i-- ) {
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
                    sol[i] = DiagonalInverse_[i] * accumulator;
                } else {
                    sol[i] += Omega_ *
                              ((DiagonalInverse_[i] * accumulator) - sol[i]);
                }
            }
        }
    }
}

/**********************************************************/

template <typename T>
void GaussSeidel<T>::Reset()
{
    
    DeleteArray( DiagonalInverse_ ); // delete diagonal inverse
    DiagonalInverse_ = NULL;
    Size_            =    0;         // reset the size of the LES
    Omega_           =  1.0;         // reset damping factor

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
