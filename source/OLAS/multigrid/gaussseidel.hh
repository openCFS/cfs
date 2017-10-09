// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

#ifndef OLAS_GAUSSSEIDEL_HH
#define OLAS_GAUSSSEIDEL_HH

/**********************************************************/

#include "OLAS/multigrid/smoother.hh"

namespace CoupledField {
/**********************************************************/

//! Gauss-Seidel smoother for AMG

/*! Gauss-Seidel smoother for AMG.
 *  This class implements both a pure Gauss-Seidel iterative
 *  solver/smoother and a damped Gauss-Seidel for the LES \f$Ax=b\f$,
 *  i.e. Successive Overrelaxation (SOR), with parameter \f$ \omega \f$.
 *  <br>
 *  \f$ \displaystyle u_{k+1} =
 *      (1-\omega) u^k - \omega D^{-1} \left( b - Lu^{k+1} - Uu^k \right)
 *  \f$
 *  <br>
 *  with a splitting of \f$A\f$ in its upper triangular, diagonal
 *  and lower triangular part: \f$ A = U + D + L \f$.
 */
template <typename T>
class GaussSeidel : public Smoother<T>
{
    public:

        GaussSeidel();
        ~GaussSeidel();

        //@{
        //! get value for the SOR damping parameter Omega
        Double& GetOmega() { return Omega_; }
        const Double& GetOmega() const { return Omega_; }

        //! set value for the SOR damping parameter Omega
        void SetOmega( const Double& omega ) {
            Omega_ = omega;
        }
        //@}
        
        //! setup of the Gauss-Seidel smoother

        /*! This function prepares the Gauss-Seidel smoother.
         *  That means that an array is created, containing the
         *  inverses of the diagonal entries of the matrix "matrix".
         *  (reimplementation of the virtual function of the base
         *  class Smoother
         *  \param matrix system matrix
         *  \return true, if the setup was successful
         */
        bool Setup( const CRS_Matrix<T>& matrix );


        //! one Gauss-Seidel/SOR step

        /*! One Gauss-Seidel/SOR step
         *  \param matrix The problem matrix A
         *  \param rhs Vector with the right hand side b
         *  \param sol Vector, where the solution x should be stored
         *         into
         *  \param direction \c Smoother::FORWARD: the unknowns are
         *         updated in ascending order starting at 1.
         *         \c Smoother::BACKWARD: vice versa. For a general
         *         description of this parameter see Smoother::Step.
         *  \param force_setup flag for the control of the setup
         *         (in class GaussSeidel the creation of array containing
         *         inverse values of the matrix' diagonal entries, for
         *         further information of this parameter see
         *         Smoother::Setup)
         */
        void Step( const CRS_Matrix<T>& matrix,
                   const Vector<T>&     rhs,
                         Vector<T>&     sol,
                   const typename Smoother<T>::Direction direction
                                        = Smoother<T>::FORWARD,
                   const bool           force_setup = true );
        
        //! resets the object to the state after its creation
        void Reset();


    protected:

        T    *DiagonalInverse_; //!< array with inverse diagonal entries
        Integer     Size_;            //!< size of the LES
        Double     Omega_;           //!< damping factor for a SOR
};

/**********************************************************/
} // namespace CoupledField

#endif // OLAS_GAUSSSEIDEL_HH
