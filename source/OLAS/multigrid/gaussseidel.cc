// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

#include "OLAS/multigrid/gaussseidel.hh"




namespace CoupledField {
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

    // Create a new array for the diagonal inverses, only if the
    // old one cannot be reused. So first check, whether the old
    // one is present and has appropriate size.
    // Note: (Size_ == 0) <==> (DiagonalInverse_ == NULL)
    // So it is sufficient to check (Size_ != matrix.GetNumRows()),
    // except of calls with an empty matrix.
    if( Size_ != (Integer)matrix.GetNumRows() ) {
        delete [] ( DiagonalInverse_ );  DiagonalInverse_  = NULL;
        DiagonalInverse_ = NULL;
    }
    // create a new array for the diagonal inverses
    if( DiagonalInverse_ == NULL ) {
        Size_ = matrix.GetNumRows();
        NEWARRAY( DiagonalInverse_, T_Mtype, Size_ );
    }

    const UInt *const pRow = matrix.GetRowPointer();
    const T_Mtype *const pDat = matrix.GetDataPointer();
    // fill the array with the inverses of the diagonal entries
    for( Integer i = 0; i < Size_; i++ ) {
        DiagonalInverse_[i] = OpType<T_Mtype>::invert(pDat[pRow[i]]);
    }

    // set flag for the prepared state
    return this->prepared_ = true;
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

    // apply one Gauss-Seidel step
    T_Vtype accumulator;
    const UInt *const RowPtr   = matrix.GetRowPointer();
    const UInt *const ColPtr   = matrix.GetColPointer();
    const T       *const DataPtr  = matrix.GetDataPointer();


    //////////////////////////////////////////////
    // normal Gauss-Seidel
    //////////////////////////////////////////////
    // Gauss Seidel FORWARD step
    for( Integer i = 0; i < Size_; i++ ) {
        accumulator = rhs[i];
        // skip the diagonal entry at index RowPtr[i] and
        // start with the first offdiagonal entry at RowPtr[i]+1
        for( UInt ij = RowPtr[i]+1; ij < RowPtr[i+1]; ij++ ) {
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

/**********************************************************/

template <typename T>
void GaussSeidel<T>::Reset()
{
    
    delete [] ( DiagonalInverse_ ); // delete diagonal inverse
    DiagonalInverse_ = NULL;
    Size_            =    0;         // reset the size of the LES
    Omega_           =  1.0;         // reset damping factor

    delete [] ( PenaltyFlags_ );
    PenaltyFlags_ = NULL;

    // call Reset() of base class
    Smoother<T>::Reset();
}

/**********************************************************/
} // namespace CoupledField

/**********************************************************/

