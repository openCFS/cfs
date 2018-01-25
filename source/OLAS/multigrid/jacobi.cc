// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

#include "OLAS/multigrid/jacobi.hh"


namespace CoupledField {
/**********************************************************/

template <typename T>
Jacobi<T>::Jacobi()
    : DiagonalInverse_( NULL ),
      Size_( 0 ),
      Omega_( 2.0/3.0 )
{
}

/**********************************************************/

template <typename T>
Jacobi<T>::~Jacobi()
{
    Reset();
}

/**********************************************************/

template <typename T>
bool Jacobi<T>::Setup( const CRS_Matrix<T>& matrix )
{
    // Create a new array for the diagonal inverses, only if the
    // old one cannot be reused. So first check, whether the old
    // one is present and has appropriate size.
    // Note: (Size_ == 0) <==> (DiagonalInverse_ == NULL)
    // So it is sufficient to check (Size_ != matrix.GetNumRows()),
    // except of calls with an empty matrix.
    if( Size_ != matrix.GetNumRows() ) {
        delete [] ( DiagonalInverse_ );  DiagonalInverse_  = NULL;
        DiagonalInverse_ = NULL;
    }
    // create a new array for the diaFindDiagonalEntriesgonal inverses
    if( DiagonalInverse_ == NULL ) {
        Size_ = matrix.GetNumRows();
        NEWARRAY( DiagonalInverse_, T, Size_ );
    }

    const UInt * const diagP = matrix.GetDiagPointer();
    const T * const dataP = matrix.GetDataPointer();

    // fill the array with the inverses of the diagonal entries
    for( Integer i = 0; i < (Integer)Size_; i++ ) {
        DiagonalInverse_[i] = T(1.0) / dataP[ diagP[i] ];
    }

    // create auxiliary vector
    auxVec_.Resize( matrix.GetNumRows() );

    // set flag for the prepared state
    return this->prepared_ = true;
}

/**********************************************************/

template <typename T>
bool Jacobi<T>::Setup( const CRS_Matrix<T>& matrix,
                       const bool *const    penalty_flags )
{

    if( Setup(matrix) ) {
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

    // call Setup, if object is not prepared or
    // preparation forced by force_setup == true
    if( !(this->IsPrepared()) || force_setup ) {
        if( false == Setup(matrix) )  return;
    }

    //////////////////////////////////////////////
    // normal Jacobi
    //////////////////////////////////////////////
    auxVec_.Init(T(0.0));

    Vector<T> r_h = rhs;
    matrix.MultSub(sol, r_h );

    for (UInt i = 0; i < Size_; ++i) auxVec_[i] = sol[i] + Omega_ * r_h[i] * DiagonalInverse_[i];

    // copy new solution on passed vector
    for( Integer i = 0; i < (Integer)Size_; i++ ) sol[i] = auxVec_[i];

}

/**********************************************************/

template <typename T>
void Jacobi<T>::Reset()
{
    
    delete [] ( DiagonalInverse_ );  DiagonalInverse_  = NULL; // delete diagonal inverse
    auxVec_.Clear();
    Size_  =   0; // reset the size of the LES
    Omega_ = 1.0; // reset damping factor


    // call Reset() of base class
    Smoother<T>::Reset();
}

/**********************************************************/
} // namespace CoupledField
